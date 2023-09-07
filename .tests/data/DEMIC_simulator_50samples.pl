use strict;
use threads;
#usage: perl DEMIC_simulator_50samples.pl random_select_genus4simulate.list /dir/of/genomes/ /dir/of/output/

my $max_thread = 20;
my $sample_number = 5; # prev 50
my $sample_occurrence_ratio = 0.6;

my $precision = 100;
my $read_length = 100;
my $substitution_rate = .003;
my $indel_rate = .0002;
my $substitution_motif_rate = .2;
my $indel_motif_rate = .52;
my @substitution_motif = ('GGG', 'CGG', 'AGG');
my @substitution_base = ('A', 'C', 'T', 'T', 'G');
my @indel_base = ('A', 'C', 'T', 'G');
my @indel_motif = ('AAA', 'TTT', 'GGG', 'CCC');
my $sequencing_quality = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJK';
my @sequencing_qualities= split '', $sequencing_quality;

my $list = $ARGV[0];
my $dir = $ARGV[1];
my $dir_out = $ARGV[2];
my %sp_info;

open LIST, "<", $list or die "cannot open $list: $!";
while (<LIST>) {
    chomp;
    my @line = split /\t/;
    my ($NC, undef) = split /[.]/, $line[1];
    my ($genus, undef) = split /\s+/, $line[0];
    my ($start, undef) = split /[.]+/, $line[3];
    $sp_info{$NC} = [$genus, $line[2], $start];
}
close LIST;


my $fa_distr_ref = &thread_distribution( $max_thread, \%sp_info );

my @ths;
for my $th_ref (@{$fa_distr_ref}) {
    # $th_ref is the reference to the array of distributed fq files in this thread
    my $th = threads -> new({'context' => 'void'}, \&fq_sim_ref, $th_ref);
    my $th_id = $th->tid();
    if ($max_thread > 1) {
        print " Worker $th_id begins to scan: \n @{$th_ref}\n";
    }
    unshift @ths, $th;
}

# Join the threads

for (@ths) {
    my $th_id = $_->tid();
    $_ -> join();
    if ($max_thread > 1) {
        print " Worker $th_id finished reporting.\n";
     }
}
@ths = ();

sub fq_sim_ref {

    my @NC_th = @{$_[0]};
    for my $NC (@NC_th) {
        my $info_ref = $sp_info{$NC};
        #while (my ($NC, $info_ref) = each %sp_info) {
        my $fa;
        opendir GENUS, $dir or die "cannot opendir $dir: $!";
        for my $file (readdir GENUS) {
            if (-f $dir.$file and $file =~ /(^.+)\.fna$/) {
                print "$1\n";
                my @NC_ids = split '__', $1;
                if ($NC_ids[1] eq $NC) {
                    $fa = $file;
                    last;
                }
            }
        }
        closedir GENUS;

        my ($seq_raw, $seq) = ('', '');
        if (!defined $fa) {
            die "cannot find fasta file for $NC: $fa!";
        } else {
            open FA, "<", $dir.$fa or die "cannot open $dir$fa: $!";
            while (<FA>) {
                unless (/^>/) {
                    chomp;
                    $seq_raw .= $_;
                }
            }
        }
        my ($genome_length, $start) = (${$info_ref}[1], ${$info_ref}[2]);
        if (length($seq_raw) != $genome_length) {
            #print ${$info_ref}[1], "\t", length($seq_raw), "\n";
            die "different genome length for $NC: $genome_length";
        } else {
            $seq = &adjustment_genome_seq($seq_raw, $genome_length, $start);
            die "different genome length after adjustment for $NC" if length($seq) != $genome_length+$read_length-1;
        }

        my $GC_content = &GC_content_calculate($seq, $genome_length);   #added

        my $GC_bias_rand = 0; #sprintf ("%.2f", 0-2-rand(2));  #added
        #my $GC_bias = &GC_bias_calculate($seq, $GC_bias_rand);

        for my $m (1 .. $sample_number) {
            my $rand_sample = rand();
            next if $rand_sample > $sample_occurrence_ratio;

            my $PTR = sprintf ("%.3f", rand(19)/10+1.1);
            my $n_average = sprintf ("%.5f", 10**(0-0.3+rand(13)/10));

            my $accu_prob_ref = &accu_prob_half($genome_length, $PTR, $GC_content, $GC_bias_rand);  ## added $GC_content, $GC_bias_rand

            #my $PTR3 = sprintf("%.3f", $PTR);
            my $n2 = sprintf("%.2f", $n_average);
            print "$NC\t$m\t$n2\t$PTR\n";

            open FQ, ">", $dir_out.$NC."_${m}_${n2}_${PTR}_$GC_bias_rand.fq" or die "cannot write $dir_out${NC}_${m}_${n2}_${PTR}_$GC_bias_rand.fq: $!";
            for my $i ( 1 .. int($n_average*$genome_length/$read_length) ) {
                my $read_start;
                my ($read_seq, $read_raw, $read_quality);
                my ($sub_loci, $indel_loci) = ('NA', 'NA');

                my $rand_first_last = rand(2);
                if(int($rand_first_last) == 1) {
                    $read_start = &rand2read_start($genome_length, ${$accu_prob_ref}[0], ${$accu_prob_ref}[1], 1);  #1 for first half
                } else {
                    $read_start = &rand2read_start($genome_length, ${$accu_prob_ref}[2], ${$accu_prob_ref}[3], 0);  #0 for last half
                }
                

                my $rand_strand = rand(2);
                if ( int($rand_strand)==0 ) {
                    $read_raw = substr($seq, $read_start, $read_length+1);
                } else {
                    if ($read_start > 0){
                        $read_raw = &comp_rev(substr($seq, $read_start-1, $read_length+1));
                    } else {
                        $read_raw = &comp_rev(substr($seq, $read_start, $read_length+1));
                    }
                }

                my $rand_subindel = rand();

                if ($rand_subindel <= $substitution_rate*$read_length){
                    my $sub_ref = &substitution_read ($read_raw);
                    $read_seq = ${$sub_ref}[0];
                    $sub_loci = ${$sub_ref}[1];
                } elsif ($rand_subindel <= ($indel_rate+$substitution_rate)*$read_length) {
                    my $indel_ref = &indel_read ($read_raw); #pos: insertion; neg: deletion
                    $read_seq = ${$indel_ref}[0];
                    $indel_loci = ${$indel_ref}[1];
                } else {
                    $read_seq = substr($read_raw, 0, $read_length);
                }
                
                print FQ '@',"${NC}_$m.$i $read_start:", int($rand_strand), ":${sub_loci}:$indel_loci\n";
                print FQ $read_seq,"\n+\n";

                $read_quality = &quality_read(length($read_seq), $sub_loci, $indel_loci);
                print FQ "$read_quality\n";

                #print $read_start, "\n";
            }
            close FQ;
            
        }

    }

}

sub substitution_read {
    my $read = substr($_[0], 0, $read_length);
    my $rand_motif = rand();
    my ($index_motif, $err_loci) = (-1, -1);
    if ($rand_motif <= $substitution_motif_rate) {
        for my $motif (@substitution_motif) {
            $index_motif = rindex($read, $motif);
            if ($index_motif >= 0) {
                my $rand_position = rand(3);             
                $err_loci = $index_motif + int($rand_position);
                last;
            }
        }
    } 
    if ($index_motif < 0) {
        my $rand_base_index = rand(scalar(@substitution_base));
        $err_loci = rindex($read, $substitution_base[int($rand_base_index)]);
        $err_loci = length($read) - 1 if $err_loci == -1;
    }
    my $ori_base = substr( $read, $err_loci, 1 );
    substr( $read, $err_loci, 1 ) = &simulate_seq_error($ori_base);
    [$read, $err_loci+1];
}

sub indel_read {
    my $read;
    my $rand_motif = rand();
    my ($index_motif, $err_loci) = (-1, -1);
    my $rand_indel = rand(2);
    if(int($rand_indel) == 0) {
        #insertion
        $read = substr($_[0], 0, $read_length-1);
    } else {
        #deletion
        $read = $_[0];
    }
    if ($rand_motif <= $indel_motif_rate) {
        for my $motif (@indel_motif) {
            $index_motif = rindex($read, $motif);
            if ($index_motif >= 0) {
                my $rand_position = rand(3);           
                $err_loci = $index_motif + int($rand_position);
                last;
            }
        }
    }
    if ($index_motif < 0) {
        my $rand_base_index = rand(scalar(@indel_base));
        $err_loci = rindex($read, $indel_base[int($rand_base_index)]);
        $err_loci = length($read) - 1 if $err_loci == -1;
    }
    
    if(int($rand_indel) == 0) {
        #insertion
        my $insert_base_rand = rand(scalar(@indel_base));
        $read = substr($read, 0, $err_loci).$indel_base[int($insert_base_rand)].substr($read, $err_loci);
        [$read, $err_loci+1];
    } else {
        #deletion
        $read = substr($read, 0, $err_loci).substr($read, $err_loci+1);
        [$read, ($err_loci+1)*(-1)];
    }
}

sub quality_read_pre {
    my ($seq_length, $sub_loci, $indel_loci) = @_;
    my @rand = map { int(rand(15))+28 } ( 1 .. $seq_length );
    if ($sub_loci ne 'NA') {
        $rand[$sub_loci-1] -= 15+int(rand(14));
    } elsif ($indel_loci>=1) {
        $rand[$indel_loci-1] -= 15+int(rand(14));
    }
    my @seq_qualities = map {$sequencing_qualities[$_]} @rand;
    my $seq_quality = join '', @seq_qualities;
    $seq_quality;
}

sub quality_read {
    my ($seq_length, $sub_loci, $indel_loci) = @_;
    my @rand = map { int(rand(15))-7+($_**2)*(0-.00706)+.606*$_+30 } ( 1 .. $seq_length );
    if ($sub_loci ne 'NA') {
        $rand[$sub_loci-1] -= 20+int(rand(14));
    } elsif ($indel_loci>=1) {
        $rand[$indel_loci-1] -= 20+int(rand(14));
    }
    @rand = map { $_ > 42 ? 42 : ($_< 0 ? 0 : $_) } @rand;
    my @seq_qualities = map {$sequencing_qualities[$_]} @rand;
    my $seq_quality = join '', @seq_qualities;
    $seq_quality;
}

sub simulate_seq_error{
    my $ori_base = $_[0];
    my @base = ('A', 'T', 'C', 'G');
    my $err_base_index;
    for my $i( 0 .. $#base ){
        if ($base[$i] =~ /$ori_base/i){
            while(1){
                $err_base_index = int(rand(4));
                last unless $err_base_index  == $i;
            }
            last;
        }
    }
    $base[$err_base_index];
}

sub adjustment_genome_seq {
    my ( $genome, $length, $oriC ) = @_;
    my $genome_ad;
    if ($oriC > int($length/2)+1) {
        my $genome_ad_first = substr ($genome, 0, $oriC - int($length/2)-1);
        $genome_ad = substr($genome, $oriC - int($length/2)-1).$genome_ad_first;
    } else {
        my $genome_ad_last = substr ($genome, $genome-(int($length/2)+1-$oriC), int($length/2)+1-$oriC);
        $genome_ad = $genome_ad_last.substr($genome, 0, $genome-(int($length/2)+1-$oriC));
    }
    $genome_ad .= substr($genome, 0, $read_length-1);

    $genome_ad;
}


sub GC_content_calculate {
    my $seq = $_[0];
    my $L = $_[1];
    my (@GC_content_prec1, @GC_content_prec2);
    my $GC_content_total = 0;
    for my $i (0 .. int($L/2/$precision-1)) {
        #the exponential probability distribution function after *100

        my $seq_local1 = substr ($seq, $i*$precision, $precision);
        $GC_content_prec1[$i] = &nt_count($seq_local1);

        my $seq_local2 = substr ($seq, $L-($i+1)*$precision, $precision);
        $GC_content_prec2[$i] = &nt_count($seq_local2);

        $GC_content_total += ($GC_content_prec1[$i]+$GC_content_prec2[$i])*$precision;
    }
    my $remaining_nt = int($L/2)-$precision*$#GC_content_prec1;
    my $seq_local_last1 = substr ($seq, $precision*($#GC_content_prec1+1), $remaining_nt);
    push @GC_content_prec1, &nt_count($seq_local_last1);
    my $seq_local_last2 = substr ($seq, $L-($precision*($#GC_content_prec1+1)+$remaining_nt), $remaining_nt);
    push @GC_content_prec2, &nt_count($seq_local_last2);
    $GC_content_total += ($GC_content_prec1[-1]+$GC_content_prec2[-1])*$remaining_nt;
    [\@GC_content_prec1, \@GC_content_prec2, $GC_content_total/length($seq)];
}

sub nt_count {
    my $seq_local = $_[0];
    my @seq_array = split ('', "\U$seq_local");
    my %nt4count;
    for my $nt (@seq_array) {
        $nt4count{$nt} ++;
    }
    ($nt4count{"G"} + $nt4count{"C"})/length($seq_local);
}


### Subroutine to calculate accumulative probability of chromosome (partial or complete) starts along half of bacterial reference genome 
### Return referrence of array for accumulative probability with index as read start/precision in each bin
### Return referrence of array for accumulative probability groups of bins with each corresponding to 1 percent of probability for speed optimization
### The last bin is usually not a complete bin (< length of precision)

sub accu_prob_half {
    #L: total length of genome
    #n_ave: average sequencing depth
    #m: designated PTR
    my ($L, $m, $GC_content, $GC_bias_rand) = @_;
    my ($GC_content1, $GC_content2, $ave_GC) = @{$GC_content};
    my (@dist_prec, @dist_prec1, @dist_prec2);
    my ($total1, $total2);
    my (@accu_prob1, @accu_prob2);
    my (@accu_prob_group1, @accu_prob_group2);

    for my $i (0 .. int($L/2/$precision-1)) {
        #the exponential probability distribution function after *100
        $dist_prec[$i] = $m**(2*($i+.5)*$precision/$L); #*($m-1)*log($m)*
        $dist_prec1[$i] = $dist_prec[$i]*( 2**($GC_bias_rand*(${$GC_content1}[$i]-$ave_GC)) );
        $dist_prec2[$i] = $dist_prec[$i]*( 2**($GC_bias_rand*(${$GC_content2}[$i]-$ave_GC)) );
        $total1 += $dist_prec1[$i];
        $total2 += $dist_prec2[$i];
        #print "$i\t$dist_prec[$i]\n";
    }

    $accu_prob1[0] = $dist_prec1[0]/$total1;
    $accu_prob2[0] = $dist_prec2[0]/$total2;
    push @{$accu_prob_group1[0]}, 0;
    push @{$accu_prob_group2[0]}, 0;
    #print "0\t$accu_prob[0]\n";
    for my $i (1 .. $#dist_prec1) {
        $accu_prob1[$i] = $accu_prob1[$i-1] + $dist_prec1[$i]/$total1;
        my $pre_n100_1 = int($accu_prob1[$i]*1000);
        push @{$accu_prob_group1[$pre_n100_1]}, $i;

        $accu_prob2[$i] = $accu_prob2[$i-1] + $dist_prec2[$i]/$total2;
        my $pre_n100_2 = int($accu_prob2[$i]*1000);
        push @{$accu_prob_group2[$pre_n100_2]}, $i;
        #print "$i\t$accu_prob[$i]\t$pre_n100\n";
    }

    [\@accu_prob1, \@accu_prob_group1, \@accu_prob2, \@accu_prob_group2];
}

### Subroutine to generate read start according to given accumulative probability and accumulative probability groups
### Input genome length, two output of subroutine accu_prob_half and whether first half
### A loop used in the subroutine to ensure generation of read start when random number was within last incomplete bin
### Return read start

sub rand2read_start {
    my $L = $_[0];
    my $tag_first_half = $_[3];
    my @accu_prob = @{$_[1]};
    my @accu_prob_group = @{$_[2]};
    
    while(1){
        my $read_start;
        my $i2;
        my $rand1 = rand();
        my @right_group = @{$accu_prob_group[int($rand1*1000)]};

        if ($rand1 < $accu_prob[ $right_group[0] ]) {
            $i2 = $right_group[0];
        } else {
            for my $n (1 .. $#right_group) {
                my $i = $right_group[$n];
                if ($rand1 >= $accu_prob[$i-1] and $rand1 < $accu_prob[$i]) {
                    $i2 = $i;
                }
            }
        }
        if (!defined $i2){
            $i2 = $right_group[-1]+1;
        }

        my $rand2 = rand();
        if ($i2 == $#accu_prob) {
            my $i_bp = int($rand2*$precision);
            if ( $i_bp <= ($L/2)%($precision) -1 ) {
                $read_start = $i2*$precision + $i_bp;
            } elsif ( $tag_first_half == 1 and $L%2 == 1 and $i_bp == (($L+1)/2)%($precision)-1 ) {
                $read_start = $i2*$precision + $i_bp;
            }

        } else {
            $read_start = $i2*$precision + int($rand2*$precision);
        }
        
        if ($tag_first_half == 0 and defined $read_start) {
            #print "$rand1\t";
            return ($L - ($read_start+1));
        } elsif ($tag_first_half == 1 and defined $read_start) {
            #print "$rand1\t";
            return $read_start;
        } 
    }
}

sub comp_rev {
    my $seq = reverse($_[0]);
    $seq =~ tr/ATCG/TAGC/;
    $seq;
}

### Subroutine to distribute files to threads according to their sizes
### Input number of threads, and a hash with file names as keys and their sizes as values
### Return the reference of an array with each element as the reference of distributed files for a thread

sub thread_distribution {
    my $t = $_[0];
    my %file_size = %{$_[1]};   # e.g. 'SRR10000' => 1999999 my @sp_length_sort = sort { ${$sp_info{$b}}[1] <=> ${$sp_info{$a}}[1] } keys %sp_info;
    my @file_distributed;   # the array to return
    my @ordered_file = sort { ${$file_size{$b}}[1] <=> ${$file_size{$a}}[1] } (keys %file_size);
    my @thread_sum;

    for my $file (@ordered_file) {
        my ($min_thread, undef) = sort { $thread_sum[$a] <=> $thread_sum[$b] } (0 .. $t-1);
        push @{$file_distributed[ $min_thread ]}, $file;
        $thread_sum[ $min_thread ] += ${$file_size{$file}}[1];
    }
    \@file_distributed;
}