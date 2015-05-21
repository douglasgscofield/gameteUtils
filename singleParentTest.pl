#!/usr/bin/env perl

# Copyright (c) 2012,2015 Douglas G. Scofield, Uppsala University
# douglasgscofield@gmail.com
#
# No warranty is implied or assumed by this code.  Please send bugs, suggestions etc.
#
use strict;
use warnings;
use POSIX qw/isdigit/;
use Getopt::Long;
use List::Util;
# http://stackoverflow.com/questions/21204733/a-better-chi-square-test-for-perl
use Statistics::Distributions qw/ chisqrprob /;

my $o_sample1;
my $o_sample2;
my $o_mincov = 10;
my $o_minalt = 2;
my $o_allsites = 0;
my $o_show_zerodiv = 1;

my $n_homozygous = 0;
my $o_max_p_val = 0.05;


my $which_bams = "";
my @which_bams = ();
my $o_quiet;
my $o_help;
my $o_N;
my $o_fai;
my $current_reference = ""; # the name of the current reference sequence
my $N_references = 0; # the number of reference sequences seen
my $N_coordinates = 0; # the total number of coordinates seen

my $usage = "
NAME

  $0 - Convert pileup to profile2 format

SYNOPSIS

  samtools mpileup -B -C50 -q1 -f test.fa test.bam | pileup2pro2.pl > test_pro2.txt

OPTIONS

    -                          read input from stdin, write to stdout
    --in FILE                  read input from FILE, else from first non-argument
                               on the command line or from stdin
    --out FILE                 write output to FILE, else to stdout
    --which-bams INT[,INT...]  produce profile output for the INT-th BAM(s) in 
                               order as provided on the samtools mpileup command
                               line, starting with 1; otherwise produce profile
                               output for all BAMs 
    --N                        include 5th count column for Ns
    --quiet                    don't print progress to stderr
    --help, -?                 help message

  Profile2 format lists bases present at each position in a reference sequence,
  with columns sequence, position, A, C, G, T:

     contig_1	1	0	2	0	0
     contig_1	2	2	0	0	0
     contig_1	3	0	2	0	0
     contig_1	4	2	0	0	0
     contig_1	5	0	0	0	2
     contig_1	6	0	0	2	0

  This script simply converts the format so any filtering on base or mapping 
  quality, etc. that you may wish to do should be done when generating the pileup.
  Pileup format created from multiple BAM files has 3 columns per BAM file; this
  script will merge all columns while creating profile output up line unless the
  --which-bam option is given.

";

sub print_usage_and_exit($) {
    my $msg = shift;
    print "$msg\n" if $msg;
    print $usage;
    exit 0;
}

GetOptions(
    "1=s"              => \$o_sample1,
    "2=s"              => \$o_sample2,
    "mincov=i"         => \$o_mincov,
    "minalt=i"         => \$o_minalt,
    "allsites"         => \$o_allsites,
    "show-zerodiv"     => sub { $o_show_zerodiv = 1 },
    "no-show-zerodiv"  => sub { $o_show_zerodiv = 0 },
    "fai=s"            => \$o_fai,
    "max-p-val=f"      => \$o_max_p_val,
    "help|?"           => \$o_help,
) or print_usage_and_exit("");

print_usage_and_exit("") if $o_help;

# fill reference sequence order hash from fai file
my %ref_order;
my $ref_index = 0;
open (my $ref, "<", $o_fai) or die "cannot open Fasta index file '$o_fai': $!";
while (<$ref>) {
    my @l = split /\t/;
    $ref_order{$l[0]} = ++$ref_index if not exists $ref_order{$l[0]};
}
print STDERR "Found ".scalar(keys(%ref_order))." reference sequences in $o_fai\n";

open (my $s1, "<", $o_sample1) or die "cannot open sample 1 profile '$o_sample1': $!";
open (my $s2, "<", $o_sample2) or die "cannot open sample 2 profile '$o_sample2': $!";

sub read_line($) {
    my $fh = shift;
    my $l = <$fh>;
    return () if ! $l;
    chomp $l;
    return split /\t/, $l;
}

my @base = qw/ A C G T /;

sub do_test($$) {
    my ($l1, $l2) = @_;
    die "inconsistent ref/pos" if $l1->[0] ne $l2->[0] or $l1->[1] != $l2->[1];
    my @d = ( [ $l1->[2], $l2->[2], $l1->[2] + $l2->[2], 0 ],   # A
              [ $l1->[3], $l2->[3], $l1->[3] + $l2->[3], 1 ],   # C
              [ $l1->[4], $l2->[4], $l1->[4] + $l2->[4], 2 ],   # G
              [ $l1->[5], $l2->[5], $l1->[5] + $l2->[5], 3 ] ); # T
    @d = sort { $a->[2] <=> $b->[2] } @d;
    my $cov = $d[2]->[2] + $d[3]->[2];

    my $o = "$base[$d[3]->[3]]/$base[$d[2]->[3]] : " .  # allele 1 / allele 2
            "$d[3]->[0]/$d[2]->[0] , $d[3]->[1]/$d[2]->[1]";  # counts sample 2 , counts sample 2

    return ($l1->[0], $l1->[1], $o, ".", ".") if $cov < $o_mincov;  # insufficient coverage
    return ($l1->[0], $l1->[1], $o, "0", "-1") if $d[2]->[1] >= $o_minalt;  # there is an allele 3
    return ($l1->[0], $l1->[1], $o, "0", "1") if $d[2]->[2] < $o_minalt;  # homozygous site
    # http://stackoverflow.com/questions/21204733/a-better-chi-square-test-for-perl
    # 
    # observed @o:
    # sample 1 low | sample 1 high
    # sample 2 low | sample 2 high

    my @o = ( [ $d[2]->[0], $d[3]->[0] ],
              [ $d[2]->[1], $d[3]->[1] ] );
    my @rowsum = ( $o[0]->[0] + $o[0]->[1],  # sample 1 low + sample 1 high
                   $o[1]->[0] + $o[1]->[1] );  # sample 2 low + sample 2 high
    my @colsum = ( $o[0]->[0] + $o[1]->[0],  # sample 1 low + sample 2 low
                   $o[0]->[1] + $o[1]->[1] );  # sample 1 high + sample 2 high
    # expected @e:
    my @e = ( [ $rowsum[0] * $colsum[0] / $cov, $rowsum[0] * $colsum[1] / $cov ],
              [ $rowsum[1] * $colsum[0] / $cov, $rowsum[1] * $colsum[1] / $cov ] );

    return ($l1->[0], $l1->[1], $o, "-1", "-1") if $e[0]->[0] == 0 or
                                                   $e[0]->[1] == 0 or
                                                   $e[1]->[0] == 0 or
                                                   $e[1]->[1] == 0;

    my $chi2 = ($o[0]->[0] - $e[0]->[0])**2 / $e[0]->[0] +
               ($o[0]->[1] - $e[0]->[1])**2 / $e[0]->[1] +
               ($o[1]->[0] - $e[1]->[0])**2 / $e[1]->[0] +
               ($o[1]->[1] - $e[1]->[1])**2 / $e[1]->[1];
    my $df = 1;
    my $prob = chisqrprob($df, $chi2);
    $chi2 = sprintf("%.6f", $chi2);
    $prob = sprintf("%.4f", $prob);
    return ($l1->[0], $l1->[1], $o, $chi2, $prob);
 }

my @l1 = read_line($s1);
my @l2 = read_line($s2);

while (@l1 and @l2) {
    if ($l1[0] eq $l2[0]) { # same reference
        if ($l1[1] == $l2[1]) { # same position
            # perform the TEST
            my @result = do_test(\@l1, \@l2); 
            @l1 = read_line($s1);
            @l2 = read_line($s2);
            if (not $o_allsites) {
                next if $result[3] =~ m/^[\.0]$/;
                next if $result[4] > $o_max_p_val;
                next if $result[3] == -1 && $result[4] == -1 && ! $o_show_zerodiv;
            }
            print STDOUT join("\t", @result), "\n";
        } elsif ($l1[1] < $l2[1]) { # advance $l1
            @l1 = read_line($s1);
        } elsif ($l2[1] < $l1[1]) { # advance $l2
            @l2 = read_line($s2);
        } else {
            die "unknown condition involving positions";
        }
    } elsif ($ref_order{$l1[0]} < $ref_order{$l2[0]}) { # advance $l1
        @l1 = read_line($s1);
    } elsif ($ref_order{$l2[0]} < $ref_order{$l1[0]}) { # advance $l2
        @l2 = read_line($s2);
    } else {
        die "unknown condition involving reference order";
    }
}

