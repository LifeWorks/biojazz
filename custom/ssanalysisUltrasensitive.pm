#-#####################################################################################
#- File:     Ultrasensitive.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Ultrasensitive;
use Class::Std::Storable;
use base qw(Scoring);
{
    use Carp;
    use Utils;
    use Globals qw($verbosity $TAG);

    use Stimulus;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
#    my %A_of :ATTR(get => 'A', set => 'A', init_arg => 'A');  # constructor must supply initialization
#    my %B_of :ATTR(get => 'B', set => 'B', default => 'yyy'); # default value is yyy

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: score_genome
    # Synopsys: Score for ultrasensitivity.
    #--------------------------------------------------------------------------------------
    sub score_genome {
        my $self = shift;
        my $genome_model_ref = shift;

        confess "ERROR: internal error, $genome_model_ref not a GenomeModel" if !$genome_model_ref->isa('GenomeModel');

        my $config_ref = $self->get_config_ref();
        my $genome_name = $genome_model_ref->get_name();
        my $work_dir = $self->get_work_dir();
        my $local_dir = $self->get_local_dir();
        my $matlab_work = $self->get_matlab_work();

        my $stats_ref = $genome_model_ref->get_stats_ref();
        if (!defined $stats_ref) {
            printn "WARNING: stats_ref is not defined for $genome_name";
            $stats_ref = {};
            $genome_model_ref->set_stats_ref($stats_ref);
        }
        my $history_ref = $genome_model_ref->get_history_ref();

        printn "Ultrasensitive::score_genome scoring genome $genome_name";

        #---------------------------------------------------------
        # INIT SCORING
        #---------------------------------------------------------
        my $elite_flag = $genome_model_ref->get_elite_flag();
        if ($elite_flag) {
            printn "Ultrasensitive::score_genome elite individual already scored, previous score=$stats_ref->{score}";
            return if ($config_ref->{rescore_elite} == 0);
            printn "Ultrasensitive::score_genome re-scoring elite individual" if $verbosity > 1;

            # Should clear stats in the evolution to prevent stats loss when rescore genomes.
            $genome_model_ref->clear_stats();
            $stats_ref->{score} = 0;
        } else {
            printn "Ultrasensitive::score_genome scoring non-elite individual..." if $verbosity > 1;
 
            # Should clear stats in the evolution to prevent stats loss when rescore genomes.
            $genome_model_ref->clear_stats();
            $stats_ref->{score} = 0;
        }

        #---------------------------------------------------------
        # CREATE I/O GENES
        #---------------------------------------------------------
        my $lg_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
                START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                regulated_concentration => $config_ref->{regulated_concentration_min}, # all-zeroes
                UNUSED => "0000",
                domains => [
                    {
                        allosteric_flag => 0,
                        RT_transition_rate => $config_ref->{RT_transition_rate_min}, # all-zeroes
                        TR_transition_rate => $config_ref->{TR_transition_rate_min}, # all-zeroes
                        RT_phi => $config_ref->{RT_phi_min}, # all-zeroes
                        protodomains => [
                            {
                                type => "bsite",
                                substrate_polarity => 0,
                                binding_profile => $config_ref->{lg_binding_profile},
                                kf_profile => "0",
                                kb_profile => "0",
                                kp_profile => "0",
                                Keq_ratio => $config_ref->{Keq_ratio_min},
                                kf_polarity_mask => "0",
                                kb_polarity_mask => "0",
                                kf_conformation_mask => "0",
                                kb_conformation_mask => "0",
                                kp_conformation_mask => "0",
                                UNUSED => "0",
                            },
                        ],
                        UNUSED => "0",
                    },
                ],
            });
        printn "lg_sequence=".$lg_sequence_ref->get_sequence() if $verbosity > 1;
        my $tg_sequence_ref = $genome_model_ref->get_gene_parser_ref()->create_sequence({
                START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                regulated_concentration => $config_ref->{TG_init}, # all-zeroes
                UNUSED => "0000",
                domains => [
                    {
                        allosteric_flag => 0,
                        RT_transition_rate => $config_ref->{RT_transition_rate_min}, # all-zeroes
                        TR_transition_rate => $config_ref->{TR_transition_rate_min}, # all-zeroes
                        RT_phi => $config_ref->{RT_phi_min}, # all-zeroes
                        protodomains => [
                            {
                                type => "msite",
                                substrate_polarity => 0,
                                binding_profile => $config_ref->{tg_binding_profile},
                                kf_profile => "0",
                                kb_profile => "0",
                                kp_profile => "0",
                                Keq_ratio => $config_ref->{Keq_ratio_min},
                                kf_polarity_mask => "0",
                                kb_polarity_mask => "0",
                                kf_conformation_mask => "0",
                                kb_conformation_mask => "0",
                                kp_conformation_mask => "0",
                                UNUSED => "0",
                            },
                        ],
                        UNUSED => "0",
                    },
                ],
            });
        printn "tg_sequence=".$tg_sequence_ref->get_sequence() if $verbosity > 1;

        #---------------------------------------------------------
        # STIMULUS/SAMPLING EQUATIONS
        #---------------------------------------------------------
        # if random delay, delay from 1/4 period to 1/2 period inclusive
        my $stimulus_sub_ref = \&{$config_ref->{stimulus}};
        my $stimulus_ref = undef;
        if ($config_ref->{stimulus} eq "ss_ramp_equation") {
            $stimulus_ref = &$stimulus_sub_ref(
                NODE => "LG0000",
                DELAY => $config_ref->{LG_delay},	
                RANGE => $config_ref->{LG_range},
                STRENGTH => $config_ref->{LG_strength},
                RAMP_TIME => $config_ref->{LG_ramp_time},
                STEPS => $config_ref->{LG_steps},
            );
        } elsif ($config_ref->{stimulus} eq "staircase_equation") {
            $stimulus_ref = &$stimulus_sub_ref(
                NODE => "LG0000",
                PERIOD => $config_ref->{LG_period},
                DELAY => $config_ref->{LG_delay},
                STRENGTH => $config_ref->{LG_strength},
                CONCENTRATION => $config_ref->{LG_concentration},
                DUTY => $config_ref->{LG_duty},
                RFTIME => $config_ref->{LG_rftime},
                STEPS => $config_ref->{LG_steps},
            );
        } else {
            confess "ERROR: unknown stimulus subroutine";
        }
        my ($lg_source_eqn, $lg_sink_eqn) = @{$stimulus_ref->{equations}};
        my @stimulus_event_times = @{$stimulus_ref->{events}};
        my @stimulus_values_list = @{$stimulus_ref->{values}};
        if ($verbosity > 1) {
            printn "Stimulus:";
            printn $lg_source_eqn;
            printn $lg_sink_eqn;
            printn join " ", @stimulus_event_times;
            printn join " ", @stimulus_values_list;
        }

        my @input_vector = @stimulus_values_list;

        #---------------------------------------------------------
        # PARSE/TRANSLATE GENOME AND I/O GENES
        #---------------------------------------------------------
        my $genome_iref = $genome_model_ref->parse(
            [
                sequence_ref => $lg_sequence_ref,
                prefix => "L",
            ],
            [
                sequence_ref => $tg_sequence_ref,
                prefix => "T",
            ],
        );
        my $parse_successful = $stats_ref->{parse_successful} = $genome_model_ref->check();

        my $history = $genome_model_ref->sprint_history(10);
        printn $history if $verbosity > 1 || $config_ref->{sprint_history};

        #############################################################################
        my $network_connectivity = $stats_ref->{network_connectivity} = 0;
        if ($parse_successful) {
            my $transcript = $genome_iref->sprint(colour_flag => 0);
            printn $transcript if $verbosity > 2 || $config_ref->{sprint_transcript};
            burp_file("$matlab_work/$genome_name.tsc", "$history\n$transcript") if $config_ref->{save_transcript};
            $genome_model_ref->translate();

            #---------------------------------------------------------
            # BUILD/PRUNE NETWORK
            #---------------------------------------------------------
            # BUILD NETWORK
            my $genome_ref = $genome_model_ref->get_parser_ref();
            $genome_ref->build_network();

            # REPORT PROTODOMAIN CONNECTIVITY
            printn "Protodomains: ".join(",", map {$_->get_name()} @{$genome_ref->get_adjacency_matrix_node_refs()->{protodomains}}) if $verbosity > 1;
            printn $genome_ref->get_adjacency_matrix_ref()->{protodomains}->[0]->sprint_matrix() if $verbosity > 2;
            printn $genome_ref->get_connectivity_matrix_ref()->{protodomains}->sprint_matrix() if $verbosity > 2;
            # REPORT GENE CONNECTIVITY
            printn "Genes: ".join(",", map {$_->get_name()} @{$genome_ref->get_adjacency_matrix_node_refs()->{genes}}) if $verbosity > 1;
            printn $genome_ref->get_adjacency_matrix_ref()->{genes}->[0]->sprint_matrix() if $verbosity > 2;
            printn $genome_ref->get_connectivity_matrix_ref()->{genes}->sprint_matrix() if $verbosity > 2;

            # PRUNE GENES
            $stats_ref->{num_pruned_genes} = scalar $genome_ref->prune_isolated_genes();

            my $gene_ref = $genome_model_ref->get_gene_parser_ref();
            my $lg_gene_ref = $gene_ref->lookup_object_instance_by_name("LG0000");
            my $tg_gene_ref = $gene_ref->lookup_object_instance_by_name("TG0000");
            my $protodomain_ref = $genome_model_ref->get_protodomain_parser_ref();
            my ($lg_protodomain_ref)= $protodomain_ref->grep_instances_by_name("LPD");
            my ($tg_protodomain_ref)= $protodomain_ref->grep_instances_by_name("TPD");

            #---------------------------------------------------------
            # SCORING: 2 pts -- LG/TG connected to anything
            #---------------------------------------------------------
            if ($lg_gene_ref->get_export_flag()) {
                $network_connectivity++;
                printn "LG is connected" if $verbosity > 1;
            }
            if ($tg_gene_ref->get_export_flag()) {
                $network_connectivity++;
                printn "TG is connected" if $verbosity > 1;
            }

            #---------------------------------------------------------
            # SCORING: 90 + 400 pts -- LG/TG subnets
            #---------------------------------------------------------
            my (@lg_subnet, @tg_subnet);   # gene subnets
            my (@tg_adjacent_kinases, @tg_adjacent_phosphatases, @lg_adjacent_protodomains);
            if ($network_connectivity == 2) { # LG/TF connected
                my (@lg_pd_subnet, @tg0_pd_subnet, @tg1_pd_subnet);   # protodomain subnets
                @lg_pd_subnet = $genome_ref->get_connected(key => "protodomains", ref => $lg_protodomain_ref);
                @tg0_pd_subnet = $genome_ref->get_connected(key => "protodomains", ref => $tg_protodomain_ref, state => 0);
                @tg1_pd_subnet = $genome_ref->get_connected(key => "protodomains", ref => $tg_protodomain_ref, state => 1);

                printn "LG protodomain connects to ".join ",", (map {$_->[2]} @lg_pd_subnet) if $verbosity > 1;
                printn "TG/0 protodomain connects to ".join ",", (map {$_->[2]} @tg0_pd_subnet) if $verbosity > 1;
                printn "TG/1 protodomain connects to ".join ",", (map {$_->[2]} @tg1_pd_subnet) if $verbosity > 1;

                # max 90 points for subnet size
                my $lg_pd_subnet_size = (@lg_pd_subnet > 30) ? 30 : @lg_pd_subnet;
                my $tg0_pd_subnet_size = (@tg0_pd_subnet > 30) ? 30 : @tg0_pd_subnet;
                my $tg1_pd_subnet_size = (@tg1_pd_subnet > 30) ? 30 : @tg1_pd_subnet;
                $network_connectivity += ($lg_pd_subnet_size + $tg0_pd_subnet_size + $tg1_pd_subnet_size);

                @tg_adjacent_kinases = $genome_ref->find_adjacent_csites($tg_protodomain_ref, 0);
                @tg_adjacent_phosphatases = $genome_ref->find_adjacent_csites($tg_protodomain_ref, 1);
                $stats_ref->{num_adjacent_kinases} = scalar(@tg_adjacent_kinases);
                $stats_ref->{num_adjacent_phosphatases} = scalar(@tg_adjacent_phosphatases);
                printn "Found ".@tg_adjacent_kinases." adjacent kinases";
                printn "Found ".@tg_adjacent_phosphatases." adjacent phosphatases";

                @lg_adjacent_protodomains = union(
                    [map {$_->[0]} $genome_ref->get_adjacent(key => "protodomains", ref => $lg_protodomain_ref)],
                );
                @lg_adjacent_protodomains = simple_difference(
                    \@lg_adjacent_protodomains,
                    [$lg_protodomain_ref]
                );
                $stats_ref->{num_receptive_protodomains} = scalar (@lg_adjacent_protodomains);

                # now use the gene subnet to determine the connectivity
                @lg_subnet = map {$_->[0]} $genome_ref->get_connected(key => "genes", ref => $lg_gene_ref);
                @tg_subnet = map {$_->[0]} $genome_ref->get_connected(key => "genes", ref => $tg_gene_ref);
                printn "LG protein connects to ".join ",", (map {$_->get_name} @lg_subnet) if $verbosity > 1;
                printn "TG protein connects to ".join ",", (map {$_->get_name} @tg_subnet) if $verbosity > 1;

                #########################################################################################
                # score -- LG/TG connected to each other
                if (grep /LG/, (map {$_->get_name()} @tg_subnet)) {
                    printn "TG0000 fans out to LG0000" if $verbosity > 1;
                    $network_connectivity += 100;
                }
                if (grep /TG/, (map {$_->get_name()} @lg_subnet)) {
                    printn "LG0000 fans out to TG0000" if $verbosity > 1;
                    $network_connectivity += 100;
                }
                if (scalar(@tg_adjacent_kinases) > 0) {
                    $network_connectivity += 100;
                }
                if (scalar(@tg_adjacent_phosphatases) > 0) {
                    $network_connectivity += 100;
                }

                $network_connectivity += 100 if $network_connectivity >= 400;  # max LG/TG connectivity score
            }
            
            if ($network_connectivity >= 500) {
                $stats_ref->{network_connected_flag} = 1;
                # exclude proteins not in LG/TG subnet from export
                my @proteins_not_in_subnet = simple_difference([$genome_model_ref->get_genes()], [union(\@lg_subnet, \@tg_subnet)]);
                map {$_->set_export_flag(0)} @proteins_not_in_subnet;
                $stats_ref->{num_protein_out_subnet} = scalar @proteins_not_in_subnet;

                #---------------------------------------------------------
                # GENERATE ANC/FACILE MODEL
                #---------------------------------------------------------
                my $anc_model = $genome_model_ref->get_genome_parser_ref()->export_anc(
                    max_external_iterations => $config_ref->{max_external_iterations},
                    max_internal_iterations => $config_ref->{max_internal_iterations},
                    max_complex_size => $config_ref->{max_complex_size},
                    max_species => $config_ref->{max_species},
                    max_csite_bound_to_msite_number => $config_ref->{max_csite_bound_to_msite_number},
                    default_max_count => $config_ref->{default_max_count},
                    default_steric_factor => $config_ref->{default_steric_factor},
                    export_graphviz => ref $config_ref->{export_graphviz} ? (join ",",@{$config_ref->{export_graphviz}}) : $config_ref->{export_graphviz},
                    equations => [$lg_source_eqn, $lg_sink_eqn],
                    matlab_ode_solver => $config_ref->{solver},
                    matlab_solver_options => ('matlab_solver_options{InitialStep} = ' . "$config_ref->{InitialStep};\n" .
                        'matlab_solver_options{AbsTol} = ' . "$config_ref->{AbsTol};\n" . 
                        'matlab_solver_options{RelTol} = ' . "$config_ref->{RelTol};\n" . 
                        'matlab_solver_options{MaxStep} = ' . "$config_ref->{MaxStep}"),
                    t_final => $config_ref->{LG_timeout},
                    t_vector =>"[t0:$config_ref->{sampling_interval}:tf]",
                    ode_event_times => (join " ", @stimulus_event_times),
                    SS_timescale => $config_ref->{SS_timescale},
                );
                burp_file("$matlab_work/$genome_name.mod", $anc_model);
                system("$ENV{ANC_HOME}/anc.pl --report=species $matlab_work/$genome_name.mod");
                my @facile_model = slurp_file("$matlab_work/$genome_name.eqn");

                $self->anc_process_species_report("$matlab_work/$genome_name.species.rpt");
                my @anc_species = $self->anc_get_species();
                $stats_ref->{num_anc_species} = @anc_species;
                printn "ANC NUM SPECIES: ".scalar(@anc_species) if $verbosity > 1;
                printn "ANC SPECIES: @anc_species" if $verbosity > 2;

                #---------------------------------------------------------
                # RUN FACILE
                #---------------------------------------------------------
                my $sampling_interval = $config_ref->{sampling_interval};
                $self->facile_run(
                    EQN_FILE => "$matlab_work/$genome_name.eqn",
                    SIM_TYPE => "matlab",
                );

                ###############################################################################
                #---------------------------------------------------------
                # SCORE COMPLEXITY
                # Basically compute the number of genes, domains, protodomains, rules
                # and put those values in account as how complex is the network
                #---------------------------------------------------------
                my $num_protodomains = @{[$anc_model =~ /ReactionSite :/g]};
                my $num_domains = @{[$anc_model =~ /AllostericStructure :/g]};
                my $num_proteins = @{[$anc_model =~ /\sStructure :/g]};
                my $num_rules = @{[$anc_model =~ /CanBindRule :/g]};
                $stats_ref->{num_rules} = $num_rules;
                printn "ANC model complexity: $num_protodomains + $num_domains + $num_proteins + $num_rules" if $verbosity >= 1;
                $stats_ref->{complexity} = $num_protodomains + $num_domains + $num_proteins + $num_rules;
                my $complexity_threshold = defined $config_ref->{complexity_threshold} ? $config_ref->{complexity_threshold} : 250;
                $stats_ref->{complexity_score} = n_hill($stats_ref->{complexity}, $complexity_threshold, 1);
                #---------------------------------------------------------
                # CHECK ANC/FACILE MODEL
                #---------------------------------------------------------
                # check that TF_0 and TF_1 are products of at least 1 reaction each
                my $num_reactions_tg_1 = grep (/(->|<-).* TG00001/, @facile_model);
                my $num_reactions_tg_0 = grep (/(->|<-).* TG00000/, @facile_model);
                $stats_ref->{num_reactions_tg_0} = $num_reactions_tg_0;
                $stats_ref->{num_reactions_tg_1} = $num_reactions_tg_1;
                $network_connectivity += 100 * ($num_reactions_tg_0 > 1 ? 1 : $num_reactions_tg_0);
                $network_connectivity += 100 * ($num_reactions_tg_1 > 1 ? 1 : $num_reactions_tg_1);

                $stats_ref->{ANC_ok_flag} = ($network_connectivity >= 700) ? 1 : 0;

                # check that number of species is less than maximum
                if (!defined $config_ref->{max_species} || $config_ref->{max_species} < 0 || $stats_ref->{num_anc_species} < $config_ref->{max_species}) {
                    $network_connectivity += 100;
                }
            }

            #########################################################################################
            #---------------------------------------------------------
            # NETWORK SIMULATION
            #---------------------------------------------------------
            # sim_flag indicates that network was successfully simulated
            # and that calculated results are valid
            my $ANC_ok_flag = $stats_ref->{ANC_ok_flag};
            $stats_ref->{sim_flag} = 0;
            if ($ANC_ok_flag) {
                $stats_ref->{sim_flag} = 1;
                #---------------------------------------------------------
                # RUN MATLAB SIM
                #---------------------------------------------------------
                printn "Ultrasensitive::score_genome: running matlab driver..." if $verbosity > 1;
                my $matlab_ref = $self->get_matlab_ref();
                $matlab_ref->cmd("clear all; ${genome_name}Driver");
                $matlab_ref->wait_on("Facile.*done");

                my @event_times = $self->matlab_get_variable(name => "event_times");
                my @event_flags = $self->matlab_get_variable(name => "event_flags");
                my $timeout_flag = $stats_ref->{timeout_flag} = (
                    ($event_times[$#event_times] == $config_ref->{LG_timeout}) ||
                    (grep {$_ != 1.0} @event_flags) != 0
                ) ? 1 : 0;

                #---------------------------------------------------------
                # PLOT RESULTS
                #---------------------------------------------------------
                if (defined $config_ref->{plot_input} && $config_ref->{plot_input}) {
                    $self->matlab_plot_complex(figure => 900,
                        complex => "LG0000",
                        title_prefix => "$genome_name",
                        filename => "$genome_name" . "_input",
                    );
                }
                if (defined $config_ref->{plot_output} && $config_ref->{plot_output}) {
                    $self->matlab_plot_complex(figure => 901,
                        complex => "TG00000",
                        title_prefix => "$genome_name",
                        filename => "$genome_name" . "_output0",
                    );
                    $self->matlab_plot_complex(figure => 902,
                        complex => "TG00001",
                        title_prefix => "$genome_name",
                        filename => "$genome_name" . "_output1",
                    );
                }
                if (defined $config_ref->{plot_species} && $config_ref->{plot_species}) {
                    $self->matlab_plot_all_complexes(
                        file_dir => "species_plot",
                    );
                }
                $self->matlab_cmd("disp('Done plotting')\n");
                $self->matlab_wait_on("Done plotting");

                #---------------------------------------------------------
                # PHASE PLOT
                #---------------------------------------------------------
                if (defined $config_ref->{plot_phase} && $config_ref->{plot_phase}) {
                    my $output_node = "TG00001";
                    $self->matlab_plot_phase(
                        figure => 904,
                        X_complex => "LG0000",
                        Y_complex => $output_node,
                        title_prefix => "$genome_name",
                        axis_ref => [0, $config_ref->{LG_range},
                            0, $config_ref->{TG_init}],
                        filename => "$genome_name" . "_phase",
                    );
                }

                $matlab_ref->cmd("save(\'$genome_name\')")
            }
        }  # if $parse_successful
        #---------------------------------------------------------
        # REMOVE FILES
        #---------------------------------------------------------
        if (defined $local_dir) {
            #`echo $local_dir/matlab/${genome_name}*     | xargs rm -f`;
 
            my $file_glob = "$matlab_work/${genome_name}*";
            my @files = glob($file_glob);
            if (@files) {
                printn "Moving @files to $work_dir/matlab" if $verbosity > 1;
                system("mv @files $work_dir/matlab");
            }
        }
    }
}



# Package BEGIN must return true value
return 1;

