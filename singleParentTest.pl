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
use FindBin;
use lib $FindBin::RealBin;  # add script directory to @INC to find BinomialTest
use BinomialTest qw/ binomial_test /;

my $o_sample1;
my $o_sample2;
my $o_fai;
my $o_total_freq_test = 1;
my $o_total_freq_prob = 0.10;
my $o_mincov = 6; # 10
my $o_minalt = 2;
my $o_max_allele_3 = 0.1;

my $o_pool_freq_prob = 0.10;
my $o_pool_test_type = "binom";

my $o_allsites = 0;
my $o_show_zerodiv = 1;

my $n_homozygous = 0;
my $o_max_p_val = 0.05;

my $o_help;
my $current_reference = ""; # the name of the current reference sequence
my $N_references = 0; # the number of reference sequences seen
my $N_coordinates = 0; # the total number of coordinates seen

my $usage = "
NAME

  $0 - test for differences in allele frequencies between 2 gamete pools from a single parent


SYNOPSIS

    singleParentTest.pl --fai reference.fa.fai --1 pool1.pro --2 pool2.pro 


OPTIONS

    --1 FILE                Sample 1 profile
    --2 FILE                Sample 2 profile
    --fai FILE              Fasta index file for reference, required to specify the order of
                            sequences in the profile files

  SITE ELIGIBILITY

    --total-freq-test       Test for deviation of total allele frequency (pool 1 plus pool 2) away
                            from 1:1, using a two-sided binomial test with expected frequency 0.5.
                            Only sites which fail this test (and thus have total ratio of
                            1:1) will be given the differential frequency test via chi-squared.
                            This test is performed by default.
    --total-freq-prob FLOAT Probability threshold to consider a site violating 1:1; a site with
                            total frequency test above this threshold is accepted for testing
                            [default $o_total_freq_prob]
    --no-total-freq-test    Suppress the above test
    --mincov INT            Minimum read coverage to consider a site [default $o_mincov]
    --max-allele-3 FLOAT    Maximum accepted frequency count for a 3rd allele [default $o_max_allele_3]

    An alternative means of determining site eligibility depends on specifying some filter limits
    based on coverage.  It seems as if the above test will be better in all circumstances, but
    perhaps that will not hold up.  The two options below are ignored unless --no-total-freq-test
    is specified.

    --minalt INT            Minimum reads for alternate allele to consider a site

  ALLELE FREQUENCY TESTS

    --pool-test-type STRING Type of test for pool frequencies: chisq or binom [default $o_pool_test_type]

    --pool-freq-prob FLOAT  Probability threshold below which a pool is considered to violate 1:1
                            [default $o_pool_freq_prob]

    --max-p-val FLOAT       Maximum P value to report, ignored for 'binom'

  OTHER OPTIONS

    --allsites              Show test result for all sites, not just those in each sample
    --show-zerodiv          Whether to show sites that had zero division errors during testing
    --no-show-zerodiv 

    --help, -?              Help message

";

sub print_usage_and_exit($) {
    my $msg = shift;
    print "$msg\n" if $msg;
    print $usage;
    exit 0;
}

GetOptions(
    "1=s"                 => \$o_sample1,
    "2=s"                 => \$o_sample2,
    "fai=s"               => \$o_fai,
    "total-freq-test"     => sub { $o_total_freq_test = 1 },
    "no-total-freq-test"  => sub { $o_total_freq_test = 0 },
    "total-freq-prob=f"   => \$o_total_freq_prob,
    "mincov=i"            => \$o_mincov,
    "max-allele-3=f"      => \$o_max_allele_3,
    "minalt=i"            => \$o_minalt,
    "pool-test-type=s"    => \$o_pool_test_type,
    "pool-test-prob=f"    => \$o_pool_freq_prob,
    "allsites"            => \$o_allsites,
    "show-zerodiv"        => sub { $o_show_zerodiv = 1 },
    "no-show-zerodiv"     => sub { $o_show_zerodiv = 0 },
    "max-p-val=f"         => \$o_max_p_val,
    "help|?"              => \$o_help,
) or print_usage_and_exit("");

print_usage_and_exit("") if $o_help or not $o_sample1 or not $o_sample2;

# fill reference sequence order hash from fai file
my %REF_ORDER;
my $ref_index = 0;
open (my $ref, "<", $o_fai) or die "cannot open Fasta index file '$o_fai': $!";
while (<$ref>) {
    my @l = split /\t/;
    $REF_ORDER{$l[0]} = ++$ref_index if not exists $REF_ORDER{$l[0]};
}
print STDERR "Found ".scalar(keys(%REF_ORDER))." reference sequences in $o_fai\n";

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

sub sorted_total_counts($$) {
    my ($l1, $l2) = @_;
    die "inconsistent ref/pos" if $l1->[0] ne $l2->[0] or $l1->[1] != $l2->[1];
    my @d = ( [ $l1->[2], $l2->[2], $l1->[2] + $l2->[2], 0 ],   # A
              [ $l1->[3], $l2->[3], $l1->[3] + $l2->[3], 1 ],   # C
              [ $l1->[4], $l2->[4], $l1->[4] + $l2->[4], 2 ],   # G
              [ $l1->[5], $l2->[5], $l1->[5] + $l2->[5], 3 ] ); # T
    return sort { $a->[2] <=> $b->[2] } @d;
}

# for binom, output is:
# 0: reference
# 1: position
# 2: test description string
# 3: total frequency probability
# both 4 and 5 are "." if 3 indicates not 1:1 allele ratio
# 4: two-pool test P value: pool 1 P * pool 2 P, "." if either pool's P > $o_pool_freq_prob
# 5: one-pool test P value: pool 1 P * pool 2 P, "." if both pools' P <= $o_pool_freq_prob
sub do_binom_test($$) {
    my ($l1, $l2) = @_;
    my @d = sorted_total_counts($l1, $l2);
    my $cov = $d[2]->[2] + $d[3]->[2];

    my $o = "$base[$d[3]->[3]]/$base[$d[2]->[3]] " .  # allele 1 / allele 2
            "binom $d[3]->[0]/$d[2]->[0] $d[3]->[1]/$d[2]->[1]";  # counts sample 1 , counts sample 2

    return ($l1->[0], $l1->[1], $o, ".", ".", ".") if $cov < $o_mincov;  # insufficient coverage

    my $totprob = sprintf("%.3f", binomial_test($d[2]->[2], $cov));

    if ($totprob < $o_total_freq_prob) {  # effectively homozygous site
        return ($l1->[0], $l1->[1], $o, $totprob, ".", ".");
    }

    if ($d[1]->[2] >= $cov * $o_max_allele_3) {  # if coverage of allele 3 too high
        return ($l1->[0], $l1->[1], $o, $totprob, "-1", "-1");
    }

    my $pool1cov = $d[2]->[0] + $d[3]->[0];
    my $pool2cov = $d[2]->[1] + $d[3]->[1];
    my $pool1prob = binomial_test($d[2]->[0], $pool1cov);
    my $pool2prob = binomial_test($d[2]->[1], $pool2cov);
    my ($twopoolprob, $onepoolprob) = sprintf("%.5f", $pool1prob * $pool2prob) x 2;

    # now check to see if we can perform a two-pool test
    if ($d[2]->[0] == $d[3]->[0] or # equal counts
        (($d[2]->[0] <=> $d[3]->[0]) + ($d[2]->[1] <=> $d[3]->[1])) != 0) { # not opposite "signs"
        if ($pool1prob > $o_pool_freq_prob or $pool2prob > $o_pool_freq_prob) {
            $twopoolprob = ".";
        }
    }
    if ($pool1prob > $o_pool_freq_prob and $pool2prob > $o_pool_freq_prob) {
        $onepoolprob = ".";
    }
    return ($l1->[0], $l1->[1], $o, $totprob, $twopoolprob, $onepoolprob);
}

# for chisq, result is:
# 0: reference
# 1: position
# 2: test description string
# 3: test statistic
# 4: probability
sub do_chisq_test($$) {
    my ($l1, $l2) = @_;
    my @d = sorted_total_counts($l1, $l2);
    my $cov = $d[2]->[2] + $d[3]->[2];

    my $o = "$base[$d[3]->[3]]/$base[$d[2]->[3]] " .  # allele 1 / allele 2
            "chisq $d[3]->[0]/$d[2]->[0] , $d[3]->[1]/$d[2]->[1]";  # counts sample 2 , counts sample 2

    return ($l1->[0], $l1->[1], $o, ".", ".") if $cov < $o_mincov;  # insufficient coverage
    return ($l1->[0], $l1->[1], $o, "0", "-1") if $d[1]->[2] >= $cov * $o_max_allele_3;  # there is an allele 3
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
            my @result;
            if ($o_pool_test_type eq 'binom') {
                @result = do_binom_test(\@l1, \@l2); 
                @l1 = read_line($s1);
                @l2 = read_line($s2);
                if (not $o_allsites) {
                    next if $result[3] < $o_total_freq_prob;
                    next if $result[4] == -1;
                }
            } elsif ($o_pool_test_type eq 'chisq') {
                @result = do_chisq_test(\@l1, \@l2); 
                @l1 = read_line($s1);
                @l2 = read_line($s2);
                if (not $o_allsites) {
                    next if $result[3] =~ m/^[\.0]$/;
                    next if $result[4] > $o_max_p_val;
                    next if $result[3] == -1 && $result[4] == -1 && ! $o_show_zerodiv;
                }
            } else {
                die "unknown test type $o_pool_test_type";
            }
            print STDOUT join("\t", @result), "\n";
        } elsif ($l1[1] < $l2[1]) { # advance $l1
            @l1 = read_line($s1);
        } elsif ($l2[1] < $l1[1]) { # advance $l2
            @l2 = read_line($s2);
        } else {
            die "unknown condition involving positions";
        }
    } elsif ($REF_ORDER{$l1[0]} < $REF_ORDER{$l2[0]}) { # advance $l1
        @l1 = read_line($s1);
    } elsif ($REF_ORDER{$l2[0]} < $REF_ORDER{$l1[0]}) { # advance $l2
        @l2 = read_line($s2);
    } else {
        die "unknown condition involving reference order";
    }
}

