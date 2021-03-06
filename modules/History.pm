#-#####################################################################################
#- File:     History.pm
#- Synopsys:
#-#####################################################################################
#- Detailed Description:
#- ---------------------
#-#####################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package History;
use Class::Std::Storable;
use base qw();
{
    use Carp;

    use Text::CSV;
    use Utils;
    use Generation;
    use Matrix;

    use Globals qw ($verbosity $TAG $config_ref);

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %last_generation_of :ATTR(get => 'last_generation', set => 'last_generation', default => 0);

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
    # Function: collect_history_from_logfile
    # Synopsys: Collect history of genomes (i.e. their scores and other attributes).
    #--------------------------------------------------------------------------------------
    sub collect_history_from_logfile {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            logfile => undef,
            @_,
        );
        check_args(\%args,1);

        my $logfile = $args{logfile};

        open (LOGFILE, "< $logfile") or die "Can't open file $logfile\n";
        my $data_dir = "$config_ref->{work_dir}/$TAG/report/report_logfile";
        `mkdir -p $data_dir`;

        my $generation;
        my $individual;
        my @attribute_names = ();
        my @attribute_values = ();
        while (<LOGFILE>) {
            my $line = $_;
            if ($line =~ /report_current_generation:\s+generation\s+(\S+)/) {
                $generation = $1 + 0;
                @attribute_names = ();
            }
            if ($line =~ /individual\s+(\S+)\s+:/) {
                $individual = $1;
                @attribute_values = ();
                while ($line =~ /(\S+)\s*=\s*(\S*\d)/g) {
                    push(@attribute_values, $2);
                    push(@attribute_names, $1);
                }
                my $file_name = sprintf("$data_dir/generation%03d_I%02d.csv", $generation, $individual);
                open my $data_file, ">> $file_name" or die "$file_name: $!";
                my $csv = Text::CSV->new({binary => 1, eol => "\n"});
                $csv->print($data_file, \@attribute_names);
                $csv->print($data_file, \@attribute_values);
                close($data_file) || warn "close failed: $!";
            }
        }
        $last_generation_of{$obj_ID} = $generation;
        close LOGFILE;
    }

    #--------------------------------------------------------------------------------------
    # Function: collect_history_from_genomes
    # Synopsys: Collect history of genomes (i.e. their scores and other attributes).
    #--------------------------------------------------------------------------------------
    sub collect_history_from_genomes {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            genomes_dir => undef,
            max_generations => -1,
            @_,
        );
        check_args(\%args,2);
        my $max_generations = $args{max_generations};
        
        my @config_attribute_names = @{$config_ref->{genome_attribute_names}};
        confess "genome_attribute_names is not specified!" if (!@config_attribute_names);

        my $generation = $config_ref->{first_generation};
        my $fossil_epoch = (defined $config_ref->{fossile_epoch} && $config_ref->{selection_method} eq "population_based_selection") ? $config_ref->{fossil_epoch} : 1;

        while(1) {
            last if ($max_generations != -1) && ($generation >= $max_generations);
            my $gen_ref;
            $gen_ref = Generation->new({});
            $gen_ref->load_generation(
                dir => $args{genomes_dir},
                number => $generation,
            );
            last if ($gen_ref->get_num_elements() == 0);
            $last_generation_of{$obj_ID} = $generation;

            # this function scans for the keys in each stats_ref of
            # each individual in the generation
            my $data_dir = "$config_ref->{work_dir}/$TAG/stats";
            my $file_name = sprintf("$data_dir/report_genome_generation%03d.csv", $generation);
            open my $data_file, ">> $file_name" or die "$file_name: $!";
            my $csv = Text::CSV->new({binary => 1, eol => "\n"});

            my $second_attribute = 'population/mutationSteps';
            if ($config_ref->{selection_method} eq "kimura_selection") {
                $second_attribute = 'mutation_attempts';
            } elsif ($config_ref->{selection_method} eq "population_based_selection") {
                $second_attribute = 'population_per_mutant';
            }

            my @intrinsic_attribute_names;
            if ($config_attribute_names[0] eq "all") {
                @intrinsic_attribute_names = $gen_ref->get_attribute_names();
            } else {
                @intrinsic_attribute_names = @config_attribute_names;
            }
            my @genome_attribute_names = ('name', $second_attribute, 'accum_mutations', 
                'accum_point_mutations', 'Stepwise_mutations', 'stepwise_point_mutations', @intrinsic_attribute_names);
            $csv->print($data_file, \@genome_attribute_names);

            my @genomes = $gen_ref->get_elements();
            my @attributes = ();
            for (my $i=0; $i < @genomes; $i++) {
                my $genome_ref = $genomes[$i];
                push(@attributes, $genome_ref->get_name());
                push(@attributes, $genome_ref->get_number());
                push(@attributes, $genome_ref->get_accum_mutations());
                push(@attributes, $genome_ref->get_accum_point_mutations());
                push(@attributes, $genome_ref->get_stepwise_mutations());
                push(@attributes, $genome_ref->get_stepwise_point_mutations());
                # Here, we output each genome stats into a line 
                # of CSV file
                for (my $j = 0; $j < scalar @intrinsic_attribute_names; $j++) {
                    my $attribute_value = $genome_ref->get_stats_ref()->{$intrinsic_attribute_names[$j]};
                    push(@attributes, $attribute_value);
                }
                # process the attributes and output
                $csv->print($data_file, \@attributes);

                # destroy @attributes
                undef @attributes;

            }
            close($data_file) || warn "close failed: $!";


            if ($config_ref->{selection_method} eq "kimura_selection") {
                $generation++;
            } elsif ($config_ref->{selection_method} eq "population_based_selection") {
                $generation += $fossil_epoch;
            }

        }
    }

    #--------------------------------------------------------------------------------------
    # Function: collect_info_from_networks
    # Synopsys: Collect structure and reaction information from ANC mod files
    #--------------------------------------------------------------------------------------
    sub collect_info_from_networks {
        my $self = shift; my $obj_ID = ident $self;
        my %args = (
            analysis_dir => undef,
            @_,
        );
        check_args(\%args,1);
        my $analysis_dir = $args{analysis_dir};
        
        my $file_glob = "$analysis_dir/matlab/*.mod";
        my @anc_files = (glob $file_glob);

        foreach my $anc_file (@anc_files) {
            print "EXTRACTION: extract interaction information from ANC file $anc_file\n" if $verbosity >= 1;
            $self->extract_network_matrix(anc_file => $anc_file, );
        }
    }


#===  FUNCTION  ================================================================
#         NAME: extract_network_matrix
#      PURPOSE: extract the network infomation from an anc model
#   PARAMETERS: anc_model as a string
#      RETURNS: save the interation info into a CSV file
#  DESCRIPTION: 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================

    sub extract_network_matrix {
        my $self = shift; my $obj_ID = ident $self;

        my %args = (
            anc_file => undef,
            @_,
        );
        check_args(\%args, 1);

        my $anc_file = $args{anc_file};
        my $analysis_dir = defined $config_ref->{analysis_dir} ? $config_ref->{analysis_dir} : "analysis";
        `mkdir -p $config_ref->{work_dir}/$analysis_dir/$TAG/networkParam`;

        my $kd_matrix_ref = Matrix->new({});
        my $kp_matrix_ref = Matrix->new({});
        my @matrix_nodes = ();
        my @matrix_pds = ();
        my @node_domains = ();
        my @node_genes = ();
        my @concentrations = ();
        my $genome_name;

        if ($anc_file =~ /\/(G\S+?_I\S+?).mod/) {
            $genome_name = $1; 
        } else {
            die "can't extract the name of the genome/ANC file $anc_file\n";
        }
        open (ANC, "< $anc_file") or die "Can't open file $anc_file\n";

        my $anc_model = "";
        # preprocessing the model file
        LINES: while (<ANC>) {
            $_ =~ s/^\s+//;         # Strip out leading whitespace
            $_ =~ s/\s+$//;         # Strip out trailing whitespace (including \n)
            $_ =~ s/\s*\#.*//;      # Strip out trailing comment and whitespace
            $_ =~ s/\s*\/\/.*//;    # Strip out trailing comment and whitespace
            next if($_ =~ /^$/);    # Skip empty lines
            my $line = $_;

            $anc_model .= "$line\n"; # add \n since closing part will always follow
        }
        close ANC;

        # now, go over the first iteration to extract all the reactant information
        while ($anc_model =~ /CanBindRule : \{\s+?name\s?=>\s?"(\S+)\s?(\S+)\s?\(\s?(\S)\s?(\S)\s?(\S)\s?(\S)\s?\)",\n.*\n.*\n.*\s+kf\s?=>\s?(\S+),\s+kb\s?=>\s?(\S+),\n/g) {
            my $pd0 = $1; my $pd1 = $2;
            my $ms0 = $3; my $rt0 = $4;
            my $ms1 = $5; my $rt1 = $6;
            my $kf = $7; my $kb = $8;
            print "Found interaction: $pd0 ($ms0, $rt0) and $pd1 ($ms1, $rt1) with kf = $kf , kb = $kb \n" if $verbosity > 1;

            my $node0 = join("_", $pd0, $ms0, $rt0); my $index0;
            my $node1 = join("_", $pd1, $ms1, $rt1); my $index1;
            if (grep {$_ eq $node0} @matrix_nodes) {
                ($index0) = grep {$matrix_nodes[$_] eq $node0} 0..$#matrix_nodes;
            } else {
                push(@matrix_nodes, $node0);
                $index0 = $#matrix_nodes;
                push(@matrix_pds, $pd0);
            }
            if (grep {$_ eq $node1} @matrix_nodes) {
                ($index1) = grep {$matrix_nodes[$_] eq $node1} 0..$#matrix_nodes;
            } else {
                push(@matrix_nodes, $node1);
                $index1 = $#matrix_nodes;
                push(@matrix_pds, $pd1);
            }

            if ($index0 < $index1) {
                $kd_matrix_ref->set_element($index0, $index1, "kf=$kf");
                $kd_matrix_ref->set_element($index1, $index0, "kb=$kb");
            } elsif ($index0 > $index1) {
                $kd_matrix_ref->set_element($index1, $index0, "kf=$kf");
                $kd_matrix_ref->set_element($index0, $index1, "kb=$kb");
            } else {
                $kd_matrix_ref->set_element($index0, $index1, "kf=$kf;kb=$kb");
            }

        }
        my $node_num = scalar @matrix_nodes;
        if (!($kd_matrix_ref->get_element($node_num-1, $node_num-1))) {
            $kd_matrix_ref->set_element($node_num-1, $node_num-1, 0);
        }

        # now consider the kp rates
        while ($anc_model =~ /CanBindRule : \{\s+?name\s?=>\s?"(\S+)\s?(\S+)\s?\(\s?(\S)\s?(\S)\s?(\S)\s?(\S)\s?\)",\n.*\n.*\n.*\n.*\n.*\s+kp\s?=>\s?(\S+),\n/g) {
            my $pd0 = $1; my $pd1 = $2;
            my $ms0 = $3; my $rt0 = $4;
            my $ms1 = $5; my $rt1 = $6;
            my $kp = $7; 
            print "Found catalytic reaction: $pd0 ($ms0, $rt0) catalyse $pd1 ($ms1, $rt1) with kp = $kp\n" if $verbosity > 1;

            my $node0 = join("_", $pd0, $ms0, $rt0); my $index0;
            my $node1 = join("_", $pd1, $ms1, $rt1); my $index1;
            if (grep {$_ eq $node0} @matrix_nodes) {
                ($index0) = grep {$matrix_nodes[$_] eq $node0} 0..$#matrix_nodes;
            } else {
                die "ERROR: $node0 is not in node array!";
            }
            if (grep {$_ eq $node1} @matrix_nodes) {
                ($index1) = grep {$matrix_nodes[$_] eq $node1} 0..$#matrix_nodes;
            } else {
                die "ERROR: $node1 is not in node array!";
            }

            $kp_matrix_ref->set_element($index0, $index1, "kp=$kp");
        }
        if (!($kp_matrix_ref->get_element($node_num-1, $node_num-1))) {
            $kp_matrix_ref->set_element($node_num-1, $node_num-1, 0);
        }

        # now consider the concentrations
        for my $PD (@matrix_pds) {
            if ($anc_model =~ /AllostericStructure : \{\s+name\s?=>\s?"(\S+)",\n.*\n.*\n.*\s+?elements\s?=>\s?\[\S*?$PD\S*?\],\s+\}/) {
                push(@node_domains, $1);
            }
            else {
                die "Can't find the domain for $PD";
            }
        }
        for my $domain (@node_domains) {
            if ($anc_model =~ /\nStructure : \{\s+name\s?=>\s?"(\S+)",\s+elements\s?=>\s?\[\S*?$domain\S*?\],\n.*\s+\}/) {
                push(@node_genes, $1);
            }
            else {
                die "Can't find the gene for $domain";
            }
        }
        for my $gene (@node_genes) {
            if ($anc_model =~ /Init : \{\s+structure\s?=>\s?$gene,\s+IC\s?=>\s?(\S+),\s+\}/) {
                push(@concentrations, $1);
            }
            else {
                die "Can't find the concentration for $gene";
            }
        }

        my $csv = Text::CSV->new({binary => 1, eol => "\n"});
        my $export_dir = "$config_ref->{work_dir}/$analysis_dir/$TAG/networkParam";

        my ($rownum, $colnum) = $kd_matrix_ref->matdim();
        if (($rownum != $colnum) || ($rownum != scalar @matrix_nodes) || ($rownum != scalar @node_genes) || ($rownum != scalar @concentrations)) {
            die "Matrix dimension is not same as nodes number";
        }

        my $filename = "$export_dir/kd_$genome_name.csv";
        open my $kd_file, ">> $filename" or die "$filename: $!";
        $kd_matrix_ref->check();
        for (my $i = 0; $i < $rownum; $i++) {
            my $row_ref = $kd_matrix_ref->get_matrix_ref()->[$i];
            $csv->print($kd_file, $row_ref);
        }
        $csv->print($kd_file, \@matrix_nodes);
        $csv->print($kd_file, \@node_domains);
        $csv->print($kd_file, \@node_genes);
        $csv->print($kd_file, \@concentrations);

        close($kd_file) || warn "close failed: $!";

        $filename = "$export_dir/kp_$genome_name.csv";
        open my $kp_file, ">> $filename" or die "$filename: $!";
        $kp_matrix_ref->check();
        for (my $i = 0; $i < $rownum; $i++) {
            my $row_ref = $kp_matrix_ref->get_matrix_ref()->[$i];
            $csv->print($kp_file, $row_ref);
        }
        $csv->print($kp_file, \@matrix_nodes);
        $csv->print($kp_file, \@node_domains);
        $csv->print($kp_file, \@node_genes);
        $csv->print($kp_file, \@concentrations);

        close($kp_file) || warn "close failed: $!";


        return 1;
    } ## --- end sub extract_network_matrix


}

sub run_testcases {
    printn "NO TESTCASES!!!";
}


# Package BEGIN must return true value
return 1;

