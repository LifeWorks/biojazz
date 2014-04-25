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
                # OUTPUT KINASE AND PHOSPHATASE
                #---------------------------------------------------------
                my @adjacent_kinase_names = map {$_->get_name()} @tg_adjacent_kinases;
                my @kinase_gene_names = map {$_->get_upper_ref()->get_upper_ref()->get_name()} @tg_adjacent_kinases;
                my @adjacent_phosphatase_names = map {$_->get_name()} @tg_adjacent_phosphatases;
                my @phosphatase_gene_names = map {$_->get_upper_ref()->get_upper_ref()->get_name()} @tg_adjacent_phosphatases;
                
                my $K1 = 0.0;
                if (scalar @adjacent_kinase_names > 0) {
                    for (my $i = 0; $i < @adjacent_kinase_names; $i++) {
                        my $pd_name = $adjacent_kinase_names[$i];
                        my $gene_name = $kinase_gene_names[$i];
                        my $protein_concentration = 0;
                        if ($anc_model =~ /Init : \{\s+structure\s?=>\s?$gene_name,\s+IC\s?=>\s?(\S+),/g) {
                            $protein_concentration = $1 + 0;
                            $stats_ref->{$gene_name} = $protein_concentration;
                        }
                        my @K1s = ();
                        while ($anc_model =~ /CanBindRule : \{\s+name\s?=>\s?\S$pd_name\s?(TPD\S+)\s?\(\s?(\S)\s?(\S)\s?(\S)\s?(\S)\s?\)\S,\n.*\n.*\n.*\s+kf\s?=>\s?(\S+),\s+kb\s?=>\s?(\S+),\s+kp\s?=>\s?(\S+),/g) {
                            my $rule_name = 'K1_'.$pd_name.'_'."$1".'_'."$2"."$3"."$4"."$5";
                            my $rule_rate = ($7 + $8) / $6;
                            $stats_ref->{$rule_name} = $rule_rate;
                            push(@K1s, $rule_rate);
                        }
                        if (scalar @K1s > 0) {
                            $K1 = $K1s[0] / $config_ref->{TG_init};
                            for (my $i = 1; $i < @K1s; $i++) {
                                $K1 *= ($K1s[$i] / $config_ref->{TG_init});
                            }
                        } else {
                            die "didn't find the rate of phosphorylation rule";
                        }
                        $K1 = $K1**(1/(scalar @K1s));
                    }
                }
                $stats_ref->{tg_K1} = $K1;

                my $K2 = 0.0;
                if (scalar @adjacent_phosphatase_names > 0) {
                    for (my $i = 0; $i < @adjacent_phosphatase_names; $i++) {
                        my $pd_name = $adjacent_phosphatase_names[$i];
                        my $gene_name = $phosphatase_gene_names[$i];
                        my $protein_concentration = 0;
                        if ($anc_model =~ /Init : \{\s+structure\s?=>\s?$gene_name,\s+IC\s?=>\s?(\S+),/g) {
                            $protein_concentration = $1 + 0;
                            $stats_ref->{$gene_name} = $protein_concentration;
                        }
                        my @K2s = ();
                        while ($anc_model =~ /CanBindRule : \{\s+name\s?=>\s?\S$pd_name\s?(TPD\S+)\s?\(\s?(\S)\s?(\S)\s?(\S)\s?(\S)\s?\)\S,\n.*\n.*\n.*\s+kf\s?=>\s?(\S+),\s+kb\s?=>\s?(\S+),\s+kp\s?=>\s?(\S+),/g) {
                            my $rule_name = 'K2_'.$pd_name.'_'."$1".'_'."$2"."$3"."$4"."$5";
                            my $rule_rate = ($7 + $8) / $6;
                            $stats_ref->{$rule_name} = $rule_rate;
                            push(@K2s, $rule_rate);
                        }
                        if (scalar @K2s > 0) {
                            $K2 = $K2s[0] / $config_ref->{TG_init};
                            for (my $i = 1; $i < @K2s; $i++) {
                                $K2 *= ($K2s[$i] / $config_ref->{TG_init});
                            }
                        } else {
                            die "didn't find the rate of phosphorylation rule";
                        }
                        $K2 = $K2**(1/(scalar @K2s));
                    }
                }
                $stats_ref->{tg_K2} = $K2;
            }
        }

        #---------------------------------------------------------
        # REMOVE FILES
        #---------------------------------------------------------
        if (defined $local_dir) {
            #`echo $local_dir/matlab/G*     | xargs rm -f`;
            #`echo $config_ref->{local_dir}/matlab/G*     | xargs rm -f`;
 
            my $file_glob = "$matlab_work/${genome_name}*";
            my @files = glob($file_glob);
            if (@files) {
                printn "Moving @files to $work_dir/matlab" if $verbosity > 1;
                system("mv @files $work_dir/matlab");
            }
        }
    }
}


sub run_testcases {
    $verbosity = 1;

    $TAG = "test";
    srand(33433);

    my $config_file = <<END;
#----------------------------------------
# CPU AND CLUSTER SETTINGS
#----------------------------------------
nice = 10
vmem = 200000000

#----------------------------------------
# WORKSPACE AND CUSTOM SCORING MODULE
#----------------------------------------
scoring_class = Ultrasensitive
work_dir = ultrasensitive

#----------------------------------------
# GENOME PARAMS
#----------------------------------------

# Scaling: all concentrations in uM, all 2nd-order rates in uM^-1 s^-1

# Genome class
radius = 2
kf_max = 1e3    # uM^-1 s^-1
kf_min = 1e-3
kb_max = 1e3
kb_min = 1e-3
kp_max = 1e3
kp_min = 1e-3

# Gene class
regulated_concentration_width = 10
gene_unused_width = 4
regulated_concentration_max = 1e3    # 1mM
regulated_concentration_min = 1e-3   # 1nM  ~ 1 molecule in prokaryote

# Domain class
RT_transition_rate_width = 10
TR_transition_rate_width = 10
RT_phi_width = 10
domain_unused_width = 4
RT_transition_rate_max = 1e2
RT_transition_rate_min = 1e-2
TR_transition_rate_max = 1e2
TR_transition_rate_min = 1e-2
RT_phi_max = 1.0
RT_phi_min = 0.0

# ProtoDomain class
binding_profile_width = 10
kf_profile_width = 20
kb_profile_width = 20
kp_profile_width = 10
Keq_profile_width = 10
protodomain_unused_width = 4
Keq_ratio_max = 1e2
Keq_ratio_min = 1e-2

#----------------------------------------
# ANC PARAMS
#----------------------------------------
max_external_iterations = -1
max_internal_iterations = -1
max_complex_size = 8
max_species = 256
max_csite_bound_to_msite_number = 1
default_max_count = 2          # this prevents polymerization (see ANC manual)
default_steric_factor = 1e3    # in micro-mol/L
#export_graphviz = nothing
export_graphviz = network,collapse_states,collapse_complexes

#----------------------------------------
# FACILE/MATLAB SETTINGS
#----------------------------------------
solver = ode23s
#solver = stoch

sampling_interval = 1
SS_timescale = 500.0

# MATLAB odeset params
InitialStep = 1e-8
AbsTol = 1e-9
RelTol = 1e-3
MaxStep = 500.0

#----------------------------------------
# SIMULATION/SCORING PARAMS
#----------------------------------------
plot_input = 1
plot_output = 1
plot_species = 0
plot_phase = 1
plot_min = -1

round_values_flag = 0

steady_state_threshold = 1000   # IC settling time
steady_state_score_threshold = 0.5

delta_threshold = 0.01          # relative measure of amplitude used to filter out integration noise
amplitude_threshold = 0.01      # absolute measure of amplitude
ultrasensitivity_threshold = 5  # ratio of 2nd step over 1st step

w_n = 1.0
w_c = 1.0
w_s = 1.0
w_a = 1.0
w_u = 1.0
w_u1 = 1.0
w_u3 = 1.0

LG_range = 10          # uM (about 6 molecules in 1e-18L vol ???)
LG_delay = ~
LG_strength = 4.0      # in Hz
LG_ramp_time = 3000
LG_steps = 3

LG_timeout = 20000

stimulus = ss_ramp_equation

hill_n = 40
hill_k = 5

TG_init = 1000  # uM
cell_volume = 1e-18             # 1e-18L --> sub-cellular volume

lg_binding_profile = 0110011010
tg_binding_profile = 1011101000

END

    burp_file("test/custom/Ultrasensitive.cfg", $config_file);

    my $scoring_ref = Ultrasensitive->new({
            node_ID => 98,
            config_file => "test/custom/Ultrasensitive.cfg",
            work_dir    => "test/custom",
            matlab_startup_options => "-nodesktop -nosplash",
        });

    printn $scoring_ref->_DUMP();

    my $config_ref = {};
    read_config($config_ref, "test/custom/Ultrasensitive.cfg");

    use GenomeModel;
    my $genome_model_ref = GenomeModel->new({
            name => "Ultrasensitive",
            Genome => {
                radius => $config_ref->{radius},
                kf_max => $config_ref->{kf_max},
                kf_min => $config_ref->{kf_min},
                kb_max => $config_ref->{kb_max},
                kb_min => $config_ref->{kb_min},
                kp_max => $config_ref->{kp_max},
                kp_min => $config_ref->{kp_min},
                Gene => {
                    regulated_concentration_width => $config_ref->{regulated_concentration_width},
                    unused_width => $config_ref->{gene_unused_width},
                    regulated_concentration_max => $config_ref->{regulated_concentration_max},
                    regulated_concentration_min => $config_ref->{regulated_concentration_min},
                    Domain => {
                        RT_transition_rate_width => $config_ref->{RT_transition_rate_width},
                        TR_transition_rate_width => $config_ref->{TR_transition_rate_width},
                        RT_phi_width => $config_ref->{RT_phi_width},
                        unused_width => $config_ref->{domain_unused_width},
                        RT_transition_rate_max => $config_ref->{RT_transition_rate_max},
                        RT_transition_rate_min => $config_ref->{RT_transition_rate_min},
                        TR_transition_rate_max => $config_ref->{TR_transition_rate_max},
                        TR_transition_rate_min => $config_ref->{TR_transition_rate_min},
                        RT_phi_max => $config_ref->{RT_phi_max},
                        RT_phi_min => $config_ref->{RT_phi_min},
                        ProtoDomain => {
                            binding_profile_width => $config_ref->{binding_profile_width},
                            kf_profile_width => $config_ref->{kf_profile_width},
                            kb_profile_width => $config_ref->{kb_profile_width},
                            kp_profile_width => $config_ref->{kp_profile_width},
                            Keq_profile_width => $config_ref->{Keq_profile_width},
                            unused_width => $config_ref->{protodomain_unused_width},
                            Keq_ratio_max => $config_ref->{Keq_ratio_max},
                            Keq_ratio_min => $config_ref->{Keq_ratio_min},
                        },
                    },
                },
            },
        });

    # CONFIGURE/CREATE GENOME
    my $lg_binding_profile = $config_ref->{lg_binding_profile};
    my $tg_binding_profile = $config_ref->{tg_binding_profile};
    my $sequence_ref = $genome_model_ref->get_genome_parser_ref()->create_sequence({
            PRE_JUNK => undef,   # undef == ""
            POST_JUNK => "0000",
            genes => [
                {
                    START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                    regulated_concentration => 1.0, # uM
                    UNUSED => "0000",
                    domains => [
                        {
                            allosteric_flag => 0,
                            RT_transition_rate => 1.00,
                            TR_transition_rate => 1.00,
                            RT_phi => 1.0,
                            protodomains => [
                                {
                                    type => "bsite",
                                    substrate_polarity => 0,
                                    binding_profile => BindingProfile->binding_complement($lg_binding_profile)->sprint(),
                                    kf_profile => "00101000100100010010",
                                    kb_profile => "11001000110000001000",
                                    kp_profile => "00011111000111110011",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11111100111111001110",
                                    kb_conformation_mask => "0",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                                {
                                    type => "bsite",
                                    substrate_polarity => 0,
                                    binding_profile => "0101111111",
                                    kf_profile => "11101101111101100000",
                                    kb_profile => "10010101100001100111",
                                    kp_profile => "11110000010000110010",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11111100111111001110",
                                    kb_conformation_mask => "0",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                            ],
                            UNUSED => "0",
                        },
                    ],
                },
                {
                    START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                    regulated_concentration => 1.0, # uM
                    UNUSED => "0000",
                    domains => [
                        {
                            allosteric_flag => 0,
                            RT_transition_rate => 1.00,
                            TR_transition_rate => 1.00,
                            RT_phi => 1.0,
                            protodomains => [
                                {
                                    type => "bsite",
                                    substrate_polarity => 0,
                                    binding_profile => BindingProfile->binding_complement("0101111111")->sprint(),
                                    kf_profile => "01111111010110111000",
                                    kb_profile => "10011001111111001000",
                                    kp_profile => "01110100110011000011",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11101101111111111000",
                                    kb_conformation_mask => "0",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                                {
                                    type => "csite",
                                    substrate_polarity => 0,
                                    binding_profile => BindingProfile->binding_complement($tg_binding_profile)->sprint(),
                                    kf_profile => "11010111000001100000",
                                    kb_profile => "11101010101110011000",
                                    kp_profile => "11001011010100010000",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11111100111111001110",
                                    kb_conformation_mask => "0",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                            ],
                            UNUSED => "0",
                        },
                    ],
                },
                {
                    START_CODE => undef, STOP_CODE => undef, # these fields will be filled in
                    regulated_concentration => 0.1, # uM
                    UNUSED => "0000",
                    domains => [
                        {
                            allosteric_flag => 0,
                            RT_transition_rate => 1.0,
                            TR_transition_rate => 1.0,
                            RT_phi => 0.0,
                            protodomains => [
                                {
                                    type => "bsite",
                                    substrate_polarity => 0,
                                    binding_profile => BindingProfile->binding_complement("0101111111")->sprint(),
                                    kf_profile => "01111111010110111000",
                                    kb_profile => "10011001111111001000",
                                    kp_profile => "01110100110011000011",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11101101111111111000",
                                    kb_conformation_mask => "0",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                                {
                                    type => "csite",
                                    substrate_polarity => 1,
                                    binding_profile => BindingProfile->binding_complement($tg_binding_profile)->sprint(),
                                    kf_profile => "11010111000001100000",
                                    kb_profile => "11101010101110011000",
                                    kp_profile => "11001011010100010000",
                                    Keq_ratio => 1.0,
                                    kf_polarity_mask => "0",
                                    kb_polarity_mask => "0",
                                    kf_conformation_mask => "11111100111111001110",
                                    kb_conformation_mask => "0",
                                    kp_conformation_mask => "0",
                                    UNUSED => "0",
                                },
                            ],
                            UNUSED => "0",
                        },
                    ],
                },
            ],
        });
    printn "sequence=".$sequence_ref->get_sequence();
    $genome_model_ref->set_sequence_ref($sequence_ref);

    # save the genome object
    use Storable qw(store retrieve);
    store($genome_model_ref, "test/custom/Ultrasensitive.obj");

    $scoring_ref->score_genome($genome_model_ref);
    printn $genome_model_ref->_DUMP();
    sleep 20;
}


# Package BEGIN must return true value
return 1;

