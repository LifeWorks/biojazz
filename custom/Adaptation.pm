#===============================================================================
#
#         FILE: Multistability.pm
#
#  DESCRIPTION: This is the module to scoring the multistable response accroding
#               to the ramping up with 5 steps from input signal
#
#        FILES: Dependent on Scoring.pm and Stimuli.pm
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Song Feng 
# ORGANIZATION: LifeWorks
#      VERSION: 1.0
#      CREATED: 12/04/2013 15:44:11
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use diagnostics;

package Multistability;
use Class::Std::Storable;
use base qw(Scoring);
{
    use Carp;
    use Utils;
    use Global qw($verbosity $TAG);

    use Stimulus;

    #==========================================================================
    # Class Attributes
    #==========================================================================

    #==========================================================================
    # Instance Methods
    #==========================================================================

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

                $stats_ref->{species_report_flag} = $self->anc_process_species_report("$matlab_work/$genome_name.species.rpt");
                if ($stats_ref->{species_report_flag} == 0) {
                    my @anc_species = $self->anc_get_species();
                    $stats_ref->{num_anc_species} = @anc_species;
                    printn "ANC NUM SPECIES: ".scalar(@anc_species) if $verbosity > 1;
                    printn "ANC SPECIES: @anc_species" if $verbosity > 2;

                    #---------------------------------------------------------
                    # OUTPUT KINASE AND PHOSPHATASE
                    #---------------------------------------------------------
                    my @adjacent_kinase_names = map {$_->get_name()} @tg_adjacent_kinases;
                    my @kinase_gene_names = map {$_->get_upper_ref()->get_upper_ref()->get_name()} @tg_adjacent_kinases;
                    my @adjacent_phosphatase_names = map {$_->get_name()} @tg_adjacent_phosphatases;
                    my @phosphatase_gene_names = map {$_->get_upper_ref()->get_upper_ref()->get_name()} @tg_adjacent_phosphatases;

                    my $min_phos = 0; my $max_phos = 0;
                    if (scalar @adjacent_kinase_names > 0) {
                        for (my $i = 0; $i < @adjacent_kinase_names; $i++) {
                            my $pd_name = $adjacent_kinase_names[$i];
                            my $gene_name = $kinase_gene_names[$i];
                            my $protein_concentration = 0;
                            if ($anc_model =~ /Init : \{\s+structure\s?=>\s?$gene_name,\s+IC\s?=>\s?(\S+),/g) {
                                $protein_concentration = $1 + 0;
                                $stats_ref->{$gene_name} = $protein_concentration;
                            }
                            my @phos_info = ();
                            while ($anc_model =~ /CanBindRule : \{\s+name\s?=>\s?\S$pd_name\s?(TPD\S+)\s?\(\s?(\S)\s?(\S)\s?(\S)\s?(\S)\s?\)\S,\n.*\n.*\n.*\s+kf\s?=>\s?\S+,\s+kb\s?=>\s?\S+,\s+kp\s?=>\s?(\S+),/g) {
                                my $rule_name = 'phos_'.$pd_name.'_'."$1".'_'."$2"."$3"."$4"."$5";
                                my $rule_rate = $6 + 0;
                                $stats_ref->{$rule_name} = $rule_rate;
                                push(@phos_info, $rule_rate);
                            }
                            my $min = 0; my $max = $min;
                            if (scalar @phos_info > 0) {
                                $min = $phos_info[0]; $max = $min;
                                for (my $i = 1; $i < @phos_info; $i++) {
                                    if ($phos_info[$i] < $min) {
                                        $min = $phos_info[$i];
                                    }
                                    if ($phos_info[$i] > $max) {
                                        $max = $phos_info[$i];
                                    }
                                }
                            } else {
                                die "didn't find the rate of phosphorylation rule $stats_ref->{rule_name}";
                            }
                            $min_phos += $min * $protein_concentration; $max_phos += $max * $protein_concentration;
                        }
                    }
                    $stats_ref->{tg_phosphorylation_min} = $min_phos;
                    $stats_ref->{tg_phosphorylation_max} = $max_phos;
                    printn "phosphorylation min: $min_phos; phosphosrylation max: $max_phos" if $verbosity >= 1;

                    my $min_dephos = 0; my $max_dephos = 0;
                    if (scalar @adjacent_phosphatase_names > 0) {
                        for (my $i = 0; $i < @adjacent_phosphatase_names; $i++) {
                            my $pd_name = $adjacent_phosphatase_names[$i];
                            my $gene_name = $phosphatase_gene_names[$i];
                            my $protein_concentration = 0;
                            if ($anc_model =~ /Init : \{\s+structure\s?=>\s?$gene_name,\s+IC\s?=>\s?(\S+),/g) {
                                $protein_concentration = $1 + 0;
                                $stats_ref->{$gene_name} = $protein_concentration;
                            }
                            my @dephos_info = ();
                            while ($anc_model =~ /CanBindRule : \{\s+name\s?=>\s?\S$pd_name\s?(TPD\S+)\s?\(\s?(\S)\s?(\S)\s?(\S)\s?(\S)\s?\)\S,\n.*\n.*\n.*\n.*\n.*\n\s+kp\s?=>\s?(\S+),/g) {
                                my $rule_name = 'dephos_'.$pd_name.'_'."$1".'_'."$2"."$3"."$4"."$5";
                                my $rule_rate = $6 + 0;
                                $stats_ref->{$rule_name} = $rule_rate;
                                push(@dephos_info, $rule_rate);
                            }
                            my $min = 0; my $max = $min;
                            if (scalar @dephos_info > 0) {
                                $min = $dephos_info[0]; $max = $min;
                                for (my $i = 1; $i < @dephos_info; $i++) {
                                    if ($dephos_info[$i] < $min) {
                                        $min = $dephos_info[$i];
                                    }
                                    if ($dephos_info[$i] > $max) {
                                        $max = $dephos_info[$i];
                                    }
                                }
                            } else {
                                die "didn't find the rate of dephosphorylation rule $stats_ref->{rule_name}";
                            }
                            $min_dephos += $min * $protein_concentration; $max_dephos += $max * $protein_concentration;
                        }
                    }
                    $stats_ref->{tg_dephosphorylation_min} = $min_dephos;
                    $stats_ref->{tg_dephosphorylation_max} = $max_dephos;
                    printn "dephosphorylation min: $min_dephos; dephosphosrylation max: $max_dephos" if $verbosity >= 1;


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
                    $network_connectivity += 200 * ($num_reactions_tg_0 > 1 ? 1 : $num_reactions_tg_0);
                    $network_connectivity += 200 * ($num_reactions_tg_1 > 1 ? 1 : $num_reactions_tg_1);


                    # check that number of species is less than maximum
                    if (!defined $config_ref->{max_species} || $config_ref->{max_species} < 0 || $stats_ref->{num_anc_species} < $config_ref->{max_species}) {
                        $network_connectivity += 100;
                    }
                }
            }

            #########################################################################################
            #---------------------------------------------------------
            # NETWORK SIMULATION
            #---------------------------------------------------------
            # sim_flag indicates that network was successfully simulated
            # and that calculated results are valid
            my $ANC_ok_flag = $stats_ref->{ANC_ok_flag} = ($network_connectivity >= 1000) ? 1 : 0;
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
                    $self->matlab_plot_all_complexes();
                }
                $self->matlab_cmd("disp('Done plotting')\n");
                $self->matlab_wait_on("Done plotting");

                #########################################################################
                #---------------------------------------------------------
                # SCORE STEADY STATE
                # Should be better if we caculate all event_times's delays
                # and add them togethoer
                #---------------------------------------------------------
                printn "computing steady-state slopes..." if $verbosity > 1;

                my $steady_state_threshold = $config_ref->{steady_state_threshold};
                my @substr_times = ($event_times[1]);
                for (my $i = 2; $i < @event_times; $i++) {
                   push(@substr_times, $event_times[$i] - $event_times[$i - 1]);
                }
                my $ss_time_max = 0;
                foreach my $ss_time (@substr_times) {
                    if ($ss_time_max < $ss_time) {
                        $ss_time_max = $ss_time;
                    }
                }
                $stats_ref->{steady_state_score} = n_hill(
                    $ss_time_max, $steady_state_threshold,1);

                #---------------------------------------------------------
                # REPORT RESULT VECTOR
                #---------------------------------------------------------
                printn "RESULT VECTOR: INPUT = LG OUTPUT = TG00001 DELAY = $config_ref->{LG_delay}" if $verbosity > 1;
                my (@pos_output_vector, @neg_output_vector);

                my @sampling_times = @event_times;
                for (my $i=0; $i < @sampling_times; $i++) {
                    my $t = $sampling_times[$i]; 
                    $pos_output_vector[$i] = $self->matlab_get_state(complex => "TG00001", t => $t);
                    $neg_output_vector[$i] = $self->matlab_get_state(complex => "TG00000", t => $t);

                    if ($config_ref->{round_values_flag}) {
                        $pos_output_vector[$i] = $self->matlab_round_value(
                            value => $pos_output_vector[$i],
                            AbsTol => $config_ref->{AbsTol},
                            RelTol => $config_ref->{RelTol},
                        );
                        $neg_output_vector[$i] = $self->matlab_round_value(
                            value => $neg_output_vector[$i],
                            AbsTol => $config_ref->{AbsTol},
                            RelTol => $config_ref->{RelTol},
                        );
                    }

                    printf("RESULT VECTOR:  t=%-6.2f input vector:  %8.4g pos_output_vector: %8.6g neg_output_vector: %8.6g\n",
                        $t, $input_vector[$i], $pos_output_vector[$i], $neg_output_vector[$i]) if $verbosity > 1;
                }

                if (!$timeout_flag) {
                    #---------------------------------------------------------
                    # SELECT BEST OUTPUT VECTOR
                    #---------------------------------------------------------
                    my $i_dy2_bottom = 1; # index at bottom of first middle step
                    my $i_dy2_top    = 2; # index at top of first middle step
                    my $i_dy4_bottom = 3;
                    my $i_dy4_bottom = 4;
                    my $i_dy2n_bottom = $#sampling_times - $i_dy2_bottom;   # index at bottom of first middle step
                    my $i_dy2n_top    = $#sampling_times - $i_dy2_top;   # index at top of first middle step
                    my $i_dy4n_bottom = $#sampling_times - $i_dy4_bottom;   # index at bottom of first middle step
                    my $i_dy4n_top    = $#sampling_times - $i_dy4_top;   # index at top of first middle step

                    my $t_bottom2 = $sampling_times[$i_dy2_bottom];
                    my $t_top2 = $sampling_times[$i_dy2_top];
                    my $t_bottom_n2 = $sampling_times[$i_dy2n_bottom];
                    my $t_top_n2 = $sampling_times[$i_dy2n_top];
                    my $t_bottom4 = $sampling_times[$i_dy4_bottom];
                    my $t_top4 = $sampling_times[$i_dy4_top];
                    my $t_bottom_n4 = $sampling_times[$i_dy4n_bottom];
                    my $t_top_n4 = $sampling_times[$i_dy4n_top];

                    my $pos2_t1 = $pos_output_vector[$i_dy2_bottom];
                    my $pos2_t2 = $pos_output_vector[$i_dy2_top];
                    my $neg2_t1 = $neg_output_vector[$i_dy2n_bottom];
                    my $neg2_t2 = $neg_output_vector[$i_dy2n_top];
                    my $pos4_t1 = $pos_output_vector[$i_dy4_bottom];
                    my $pos4_t2 = $pos_output_vector[$i_dy4_top];
                    my $neg4_t1 = $neg_output_vector[$i_dy4n_bottom];
                    my $neg4_t2 = $neg_output_vector[$i_dy4n_top];

                    if ($config_ref->{round_values_flag}) {
                        ($pos2_t1, $pos2_t2, $neg2_t1, $neg2_t2, $pos4_t1, $pos4_t2, $neg4_t1, $neg4_t2) = map {
                        $self->matlab_round_value(value => $_, AbsTol => $config_ref->{AbsTol}, RelTol => $config_ref->{RelTol})
                        } ($pos2_t1, $pos2_t2, $neg2_t1, $neg2_t2, $pos4_t1, $pos4_t2, $neg4_t1, $neg4_t2);
                    }
                    my $pos2_delta = ($pos2_t2 - $pos2_t1)/($pos2_t2 + $pos2_t1)*2; # relative measure the time points (how fast reach steady state?)
                    my $neg2_delta = ($neg2_t2 - $neg2_t1)/($neg2_t2 + $neg2_t1)*2;
                    my $pos4_delta = ($pos4_t2 - $pos4_t1)/($pos4_t2 + $pos4_t1)*2; # relative measure the time points (how fast reach steady state?)
                    my $neg4_delta = ($neg4_t2 - $neg4_t1)/($neg4_t2 + $neg4_t1)*2;

                    my $pos_output_select = ((($pos2_delta >= $neg2_delta) ? 1 : 0) + (($pos4_delta >= $neg4_delta) ? 1 : 0) > 0) ? 1 : 0;
                    my $output_complex = $pos_output_select ? "TG00001" : "TG00000";
                    my $delta = $stats_ref->{delta} = $pos_output_select ? ($pos2_delta + $pos4_delta)/2 : ($neg2_delta + $neg4_delta)/2;
                    my $delta_score = $stats_ref->{delta_score} = p_hill($stats_ref->{delta}, $config_ref->{delta_threshold}, 1);

                    if ($delta <= 0) {
                        printn "\nWARNING: pos_delta and neg_delta are both zero or negative\n";
                        $stats_ref->{sim_flag} = 0;
                    }
                    printn "pos_output_select = $pos_output_select" if $verbosity > 1;
                    printn "pos_delta = $pos_delta, neg_delta = $neg_delta" if $verbosity > 1;
                    my @output_vector = $pos_output_select ? @pos_output_vector : @neg_output_vector;

                    #---------------------------------------------------------
                    # PHASE PLOT
                    #---------------------------------------------------------
                    if (defined $config_ref->{plot_phase} && $config_ref->{plot_phase}) {
                        my $output_node = $pos_output_select ? "TG00001" : "TG00000";
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

                    ######### multistability measure should be multistability measure!

                    ######################################################################################################
                    # multistability measure
                    my @output_vector_slopes = map {$output_vector[$_] - $output_vector[$_-1]} (1..$#output_vector);
                    my $LG_steps = $config_ref->{LG_steps};
                    confess "ERROR: LG_steps must be odd" if (int $LG_steps/2) == ($LG_steps/2);
                    #my $dy1 = max_numeric(0, $output_vector_slopes[$i_dy2_bottom-1]);
                    my $dy2 = max_numeric(0, $output_vector_slopes[$i_dy2_top - 1]);
                    my $dy4 = max_numeric(0, $output_vector_slopes[$i_dy4_top - 1]);
                    #my $dy3 = max_numeric(0, $output_vector_slopes[$i_dy2_top]);
                    #my $dy1n = max_numeric(0, -$output_vector_slopes[$i_dy2n_bottom]);
                    my $dy2n = max_numeric(0, - $output_vector_slopes[$i_dy2n_top]);
                    my $dy4n = max_numeric(0, - $output_vector_slopes[$i_dy4n_top]);
                    #my $dy3n = max_numeric(0, -$output_vector_slopes[$i_dy2n_top-1]);

                    my ($dy1_min, $dy1_max)   = $self->matlab_get_state_range(
                        complex => $output_complex,
                        t1 => $sampling_times[$i_dy2_bottom-1],
                        t2 => $sampling_times[$i_dy2_bottom],
                    );
                    my ($dy1n_min, $dy1n_max) = $self->matlab_get_state_range(
                        complex => $output_complex,
                        t1 => $sampling_times[$i_dy2n_bottom],
                        t2 => $sampling_times[$i_dy2n_bottom+1],
                    );
                    my ($dy3_min, $dy3_max)   = $self->matlab_get_state_range(
                        complex => $output_complex,
                        t1 => $sampling_times[$i_dy2_top],
                        t2 => $sampling_times[$i_dy2_top+1],
                    );
                    my ($dy3n_min, $dy3n_max) = $self->matlab_get_state_range(
                        complex => $output_complex,
                        t1 => $sampling_times[$i_dy2n_top-1],
                        t2 => $sampling_times[$i_dy2n_top],
                    );
                    my ($dy5_min, $dy5_max)   = $self->matlab_get_state_range(
                        complex => $output_complex,
                        t1 => $sampling_times[$i_dy4_top],
                        t2 => $sampling_times[$i_dy4_top+1],
                    );
                    my ($dy5n_min, $dy5n_max) = $self->matlab_get_state_range(
                        complex => $output_complex,
                        t1 => $sampling_times[$i_dy4n_top-1],
                        t2 => $sampling_times[$i_dy4n_top],
                    );
                    printn "dy1_min/max = ($dy1_min, $dy1_max) dy1n_min/max = ($dy1n_min, $dy1n_max)" if $verbosity > 1;
                    printn "dy3_min/max = ($dy3_min, $dy3_max) dy3n_min/max = ($dy3n_min, $dy3n_max)" if $verbosity > 1;
                    printn "dy5_min/max = ($dy5_min, $dy5_max) dy5n_min/max = ($dy5n_min, $dy5n_max)" if $verbosity > 1;
                    my $dy1 =  $dy1_max  - $dy1_min;
                    my $dy1n = $dy1n_max - $dy1n_min;
                    my $dy3 =  $dy3_max  - $dy3_min;
                    my $dy3n = $dy3n_max - $dy3n_min;
                    my $dy5 =  $dy5_max  - $dy5_min;
                    my $dy5n = $dy5n_max - $dy5n_min;

                    printn "dy1 = $dy1 dy2 = $dy2 dy3 = $dy3 dy4 = $dy4 dy5 = $dy5" if $verbosity > 1;
                    printn "dy1n = $dy1n dy2n = $dy2n dy3n=$dy3n dy4n = $dy4n dy5n = $dy5n" if $verbosity > 1;

                    my $max_dy = $config_ref->{TG_init};

                    my $amplitude = 0;
                    my $step1_bottom_dy = 0;
                    my $step1_top_dy = 0;
                    my $step2_bottom_dy = 0;
                    my $step2_top_dy = 0;
                    if ($delta > 0 && $max_dy != 0 && $dy2 > 0 && $dy2n > 0 && $dy4 > 0 && $dy4n > 0) {
                        $amplitude = ($dy2 + $dy2n + $dy4 + $dy4n)/4/$max_dy;
                        my $mean_dy1 = ($dy1 + $dy1n)/2;
                        my $mean_dy2 = ($dy2 + $dy2n)/2;
                        my $mean_dy3 = ($dy3 + $dy3n)/2;
                        my $mean_dy4 = ($dy4 + $dy4n)/2;
                        my $mean_dy5 = ($dy5 + $dy5n)/2;
                        # adding 0.01*mean_dy2 means you score at most 100
                        $step1_bottom_dy = $mean_dy2/($mean_dy1 + 0.001);
                        $step1_top_dy = $mean_dy2/($mean_dy3 + 0.001);
                        $step2_bottom_dy = $mean_dy4/($mean_dy3 + 0.001);
                        $step2_top_dy = $mean_dy4/($mean_dy5 + 0.001);
                    }
                    $stats_ref->{amplitude_score} = $amplitude;
                    my $w_m1_bot = $config_ref->{w_m1_bot};
                    my $w_m1_top = $config_ref->{w_m1_top};
                    my $w_m2_bot = $config_ref->{w_m2_bot};
                    my $w_m2_top = $config_ref->{w_m2_top};
                    $stats_ref->{multistability_score} = (
                        (p_hill($step1_bottom_dy, $config_ref->{ultrasensitivity_threshold}, 1)**$w_m1_bot) *
                        (p_hill($step1_top_dy, $config_ref->{ultrasensitivity_threshold}, 1)**$w_m1_top) *
                        (p_hill($step2_bottom_dy, $config_ref->{ultrasensitivity_threshold}, 1)**$w_m2_bot) *
                        (p_hill($step2_top_dy, $config_ref->{ultrasensitivity_threshold}, 1)**$w_m2_top)
                    )**(1/($w_m1_bot + $w_m1_top + $w_m2_bot + $w_m2_top));
 
                    # check if any concentrations are negative
                    foreach my $sample (@output_vector) {
                        if ($sample < 0) {
                            printn "\nWARNING: detected negative concentrations in output vector\n";
                            $stats_ref->{sim_flag} = 0;
                        }
                    }

                    if (($stats_ref->{mean_squared_err_score} < 0) || ($stats_ref->{mean_squared_err_score} > $stats_ref->{max_mean_squared_err})) {
                        #		if (($stats_ref->{mean_squared_err_score} < 0) || ($stats_ref->{mean_squared_err_score} > 1)) {
                        # numerics got messed up, set score to zero
                        printn "\nWARNING: computed mean_squared_error_score out of bounds\n";
                        $stats_ref->{sim_flag} = 0;
                    }
                }
            }
        }  # if $parse_successful
        $stats_ref->{network_connectivity} = $network_connectivity;
        $stats_ref->{network_score} = p_hill($stats_ref->{network_connectivity}, 1000, 1);

        #---------------------------------------------------------
        # FINAL SCORE
        #---------------------------------------------------------
        my $final_score = 0;
        my $network_score = $stats_ref->{network_score} = $stats_ref->{network_score} || 0;
        my $complexity_score = $stats_ref->{complexity_score} = $stats_ref->{complexity_score} || 0;
        my $steady_state_score = $stats_ref->{steady_state_score} = $stats_ref->{steady_state_score} || 0;
        my $amplitude_score = $stats_ref->{amplitude_score} = $stats_ref->{amplitude_score} || 0;
        my $multistability_score = $stats_ref->{multistability_score} = $stats_ref->{multistability_score} || 0;


        if ($parse_successful) {
            my $w_n = $config_ref->{w_n};
            my $w_c = $config_ref->{w_c};
            my $w_s = $config_ref->{w_s};
            my $w_a = $config_ref->{w_a};
            my $w_m = $config_ref->{w_m};

            # is the input connected to the output?
            my $g0  = $stats_ref->{network_connected_flag} ? 1 : 0;
            my $g0n = !$g0 ? 1 : 0;
            # ANC network is OK (i.e. max species not reached) and was therefore simulated?
            # Also, no timeout occurred while simulating to steady-state?
            # No strange numerical problems?
            my $g1  = !$stats_ref->{timeout_flag} && $stats_ref->{sim_flag} ? 1 : 0;
            # is the output amplitude large enough that the output can be trusted artifact-free?
            my $g3  = $g1 && ($stats_ref->{delta_score} > 0.5) ? 1 : 0;
            my $g3n = !$g3 ? 1 : 0;

            # don't optimize network_score once the network is connected
            $final_score =  ($network_score * $g0n + $g0)**$w_n;
            # optimize complexity only if the network is connected
            $final_score *= (1e-3 + $complexity_score * $g0)**$w_c;
            # optimize amplitude if ANC output ok and no timeout during simulation
            $final_score *= (1e-6 + $amplitude_score * $g1)**$w_a;
            # optimize multistability if ANC output ok and no timeout during simulation
            $final_score *= (1e-6 + $multistability_score * $g1)**$w_m;

            $final_score = $final_score**(1/($w_n + $w_c + $w_a + $w_m));  # re-normalization
        }


        # prevent neg've scores
        if ($final_score < 0) {
            printn "\nWARNING: negative score detected, setting to 0\n";
            $final_score = 0;
        }

        $stats_ref->{score} = $final_score;

        $genome_model_ref->set_score($stats_ref->{score});

        printn $genome_model_ref->sprint_stats();
        printn "final_score = $final_score";

        #---------------------------------------------------------
        # REMOVE FILES
        #---------------------------------------------------------
        if (defined $local_dir) {
            `echo $local_dir/matlab/${genome_name}*     | xargs rm -f`;
            #my $file_glob = "$matlab_work/${genome_name}*";
            #my @files = glob($file_glob);
            #if (@files) {
            #    printn "Moving @files to $work_dir/matlab" if $verbosity > 1;
            #    system("mv @files $work_dir/matlab");
            #}
        }

        return 1;
    } ## --- end sub score_genome
}
 
