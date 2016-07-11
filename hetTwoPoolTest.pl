#!/usr/bin/env perl

# DONE: consistent allele ordering across pools
# DONE: after the above, add allele frequencies
# TODO: chr only
# TODO: two-tailed unselected probability
# TODO: probability output: was %0.5f but adjust during output??
# TODO: handle non-0/1 heterozygotes
#
# BUG? chr1 82219 seeming three alleles

# Copyright (c) 2016 Douglas G. Scofield, Uppsala University
# douglas.scofield@ebc.uu.se, douglasgscofield@gmail.com
#
# No warranty is implied or assumed by this code.  Please send bugs, suggestions etc.
#
use strict;
use warnings;
use Carp;
use POSIX qw/isdigit log10/;
use Getopt::Long;
use List::Util;
use FindBin;
use lib $FindBin::RealBin;  # add script directory to @INC to find BinomialTest
use BinomialTest qw/ binomial_test /;

my $o_fai_file;
my $o_hets_file;       # $h will be file handle for $o_hets_file
my $o_unselected_file; # $u will be file handle for $o_unselected_file
my $o_unselected_test = 1;
my $o_unselected_prob = 0.01;
my $o_selected_1_file; # $s1 will be file handle for $o_selected_1_file
my $o_selected_2_file; # $s2 will be file handle for $o_selected_2_file
my $o_selected_prob = 0.10;

my $o_null_fraction = 0.000;

my $o_mincov = 8; # 10
my $o_genotype = 1;  # check for incompatible genotypes?
my $o_max_allele_3 = 0.1;


my $o_log10_p = 1;
my $o_sort_bases = 1;
my $o_allsites = 1;
my $o_chr_only = 1;

my $n_homozygous = 0;
my $o_max_p_val = 0.05;

my $o_help;
my $current_reference = ""; # the name of the current reference sequence
my $N_references = 0; # the number of reference sequences seen
my $N_coordinates = 0; # the total number of coordinates seen
my $N_indels = 0; # total number of indels skipped
my $o_indels = 0;
my $indels;  # $indels is file handle for indels file
my %INDELS;

my $usage_short = "
    $0 --fai reference.fa.fai --hets h.vcf --unselected u.pro --selected-1 s1.pro --selected-2 s2.pro

        --fai FILE
        --hets FILE
        --unselected FILE
        --unselected-test FLOAT  [default $o_unselected_prob]
        --selected-1 FILE
        --selected-2 FILE
        --selected-test FLOAT    [default $o_selected_prob]
        --null-fraction FLOAT    [default $o_null_fraction]
        --mincov INT             [default $o_mincov]
        --max-allele-3 FLOAT     [default $o_max_allele_3]
        --log10-p                [default $o_log10_p]
        --chr-only               [default $o_chr_only]
        --indels                 [default $o_indels]
        --allsites               [default $o_allsites]
        --sort-bases             [default $o_sort_bases]
        --no-sort-bases
        --help, -?               More detailed help message

";

my $usage = "
NAME

  hetTwoPoolTest.pl - test list of heterozygous sites for deviations from 1:1 
                      in an 'unselected' pool and two 'selected' pools


SYNOPSIS

    $0 --fai ref.fa.fai --hets h.vcf --unselected u.pro --selected-1 s1.pro --selected-2 s2.pro


OPTIONS

    --fai FILE    Fasta index file for reference, required to specify the order of
                  sequences in the VCF and profile files

    --hets FILE   VCF-format file specifying heterozygous sites. Use care in
                  producing this file; for simplicity, this script only uses
                  presence-absence positional information from the file to
                  determine whether a site is heterozygous, it explicitly does
                  *not* use any sort of genotype quality, PASS/FAIL, etc.  Any
                  sort of filtering you would wish to provide must be done
                  before this file is read by this script.  One advantage of
                  the simplification is that the file needn't be in a fully
                  compliant VCF format, so long as the first five tab-separated
                  columns of lines that do not begin with '#' contain the same
                  information as is in the '#CHROM POS ID REF ALT' columns of a
                  VCF.  The ID isn't used and may be '.', but REF and ALT are
                  used to check the alleles in the unselected pool.

    --unselected FILE   Profile base content for 'unselected' pool
 
    --unselected-test FLOAT

                  Test for deviation of 'unselected' pool allele frequency away
                  from 1:1, using a two-sided binomial test with expected
                  frequency 0.5.  This is applied to the pool specified with
                  the --unselected option.  This test is always performed and
                  its result is always part of the output, so it is unnecessary
                  to specify this option unless you want to change FLOAT, which
                  is the probability threshold for the test. Sites with
                  binomial probabilities <= FLOAT are marked as having failed
                  this test and thus show allele frequency skew.  The 6th
                  column of output indicates the test disposition, and is
                  prefixed with 'het_' for each site that is a called
                  heterozygote.  [default $o_unselected_prob]

    --selected-1 FILE   Profile base content for selected pool 1
 
    --selected-2 FILE   Profile base content for selected pool 2

    --selected-test FLOAT

                  Test for deviation of 'selected' pool allele frequencies away
                  from 1:1, using a two-sided binomial test with expected
                  frequency 0.5.  This is applied to each of the selected pools
                  specified with the --selected-# options.
                  Sites with binomial probabilities <= FLOAT are marked
                  as having failed this test and thus show allele frequency
                  skew.  The 6th column of output indicates the test
                  disposition, and is prefixed with 'het_' for each site that
                  is a called heterozygote.  [default $o_selected_prob]


    --null-fraction FLOAT

                  For fraction FLOAT sites that are *not* called heterozygotes,
                  perform the unselected and selected tests as above.  Output
                  This fraction is also subject to the --mincov and
                  --max-allele-3 settings.  The 5th column of output indicates
                  the test disposition, and is prefixed with 'null_' for each
                  site tested as a result of this option.  [default
                  $o_null_fraction]

    --mincov INT  Minimum read coverage to consider a site, applied to each of
                  the unselected and selected pools [default $o_mincov].

    --max-allele-3 FLOAT

                  Maximum accepted frequency count for a 3rd allele, applied to
                  each of the unselected and selected pools [default $o_max_allele_3]

  OTHER OPTIONS

    --log10-p     Output log10-d P values for two-pool test [default $o_log10_p]

    --chr-only    Only sites on reference sequences beginning with 'chr'

    --indels      Write apparent indels to an indels output file [default $o_indels]

    --allsites    Show test result for all heterozygous sites, not just those
                  for which the unselected pool passes the unselected test

    --sort-bases
    --no-sort-bases
                  In output, sort, or do not sort, bases in lexical order (A,C,G,T)
                  The default is ".($o_sort_bases?"":"not")." to sort

    --help, -?    Help message

OUTPUT

Fourteen tab-separated columns:

   1 reference
   2 position
   3 read coverage in unselected pool
   4 test description string, unselected pool
   5 unselected probability
   6 unselected test result string
   7 read coverage in selected 1 pool
   8 read coverage in selected 2 pool
   9 test description strings from selected 1 and 2 pools
  10 unselected test results from selected 1 and 2 pools
  11 probabilities from selected 1 and 2 pools
  12 -1/0/+1 comparisons of allele counts from selected pools
  13 twopool test result string
  14 twopool test multiplied probability

'Unselected test result string' contains one of the following, prefixed with
'het_' for a heterozygous site and 'null_' for a non-heterozygous site selected
via --null-fraction.

    zerocov    : no pool coverage at this site; test probability is reported as '.'
    mincov     : pool coverage is below the --mincov value
    genotype   : genotype incompatibility (REF and ALT do not match top pool alleles)
    allele3    : coverage frequency of a 3rd allele is greater than --max-allele-3
    binom_fail : binomial test probability is less than --unselected-prob
    binom_pass : binomial test probability is >= --unselected-prob

'Twopool test result string' contains one of the following, prefixed with
'het_' for a heterozygous site and 'null_' for a non-heterozygous site selected
via --null-fraction.  For those named similarly to those above, the meaning is similar.

    zerocov      : one/both selected pools has zero coverage; probability reported as '.'
    mincov       : one/both selected pools has coverage below --mincov
    genotype     : one/both selected pools has genotype incompatibility
    allele3      : one/both selected pools has 3rd allele frequency greater than --max-allele-3
    counts       : the selected pools do not have opposed allele counts (see field 12)
    {fail,pass}prob{1,2} : whether the binomial probabilities of the selected pools
                   are > --selected-prob or <= --selected-prob, if at least one fails
    twopool_pass : both selected pools have opposed allele counts and have binomial
                   probabilities <= --selected-prob

";

sub print_usage_and_exit {
    my $msg = join(" ", @_);
    print "$msg\n" if $msg;
    print $usage;
    exit 0;
}

if (scalar(@ARGV) == 0) {
    print $usage_short;
    exit 1;
}

GetOptions(
    "fai=s"               => \$o_fai_file,
    "hets=s"              => \$o_hets_file,
    "unselected=s"        => \$o_unselected_file,
    "unselected-test=f"   => \$o_unselected_prob,
    "selected-1=s"        => \$o_selected_1_file,
    "selected-2=s"        => \$o_selected_2_file,
    "selected-test=f"     => \$o_selected_prob,
    "null-fraction=f"     => \$o_null_fraction,
    "mincov=i"            => \$o_mincov,
    "max-allele-3=f"      => \$o_max_allele_3,
    "log10-p"             => \$o_log10_p,
    "chr-only"            => \$o_chr_only,
    "no-chr-only"         => sub { $o_chr_only = 0 },
    "allsites"            => \$o_allsites,
    "sort-bases"          => \$o_sort_bases,
    "no-sort-bases"       => sub { $o_sort_bases = 0 },
    "help|?"              => \$o_help,
) or print_usage_and_exit();

print_usage_and_exit() if $o_help 
                          or not $o_fai_file
                          or not $o_hets_file
                          or not $o_unselected_file
                          or not $o_selected_1_file
                          or not $o_selected_2_file;

my @base = qw/ A C G T /;
my %base; $base{A} = 0; $base{C} = 1; $base{G} = 2; $base{T} = 3;

# fill reference sequence order hash from fai file
my %REF_ORDER;
my $ref_index = 0;
open (my $ref, "<", $o_fai_file) or die "cannot open Fasta index file '$o_fai_file': $!";
while (<$ref>) {
    my @l = split /\t/;
    $REF_ORDER{$l[0]} = ++$ref_index if not exists $REF_ORDER{$l[0]};
}
print STDERR "Found ".scalar(keys(%REF_ORDER))." reference sequences in $o_fai_file\n";

sub open_possibly_gzipped($) {
    my ($file, $fh) = shift;
    if ($file =~ /\.gz$/) {
        open ($fh, "-|", "gzip -dc $file") or die "cannot open '$file': $!";
    } else {
        open ($fh, "<", $file) or die "cannot open '$file': $!";
    }
    return $fh;
}

my $h = open_possibly_gzipped($o_hets_file);

if ($o_indels) {
    $o_indels = $o_hets_file . "_indels.txt";
    open ($indels, ">", $o_indels) or die "cannot open indels output file '$o_indels': $!"; 
}

my $u = open_possibly_gzipped($o_unselected_file);
my $s1 = open_possibly_gzipped($o_selected_1_file);
my $s2 = open_possibly_gzipped($o_selected_2_file);

sub read_hets_line() {
    READ_LINE:
    my $l = <$h>;
    while ($l and $l =~ /^#/) {
        $l = <$h>;
    }
    if ($o_chr_only) {
        while ($l and $l =~ /^chr/) {
            $l = <$h>;
        }
    }
    return () if ! $l;
    chomp $l;
    #print STDERR "read_hets_line: '$l'\n";
    my @l = split /\t/, $l;
    if (length($l[3]) > 1 or length($l[4]) > 1) {  # ref or alt not a single base
        ++$N_indels;
        ++$INDELS{"$l[0]:$l[1]"};
        print { $indels } "INDEL\t", join("\t", @l[(0,1,3,4)]), "\n" if $o_indels;
        goto READ_LINE;
    }
    # return CHROM POS REF ALT
    ($l[3], $l[4]) = ($l[4], $l[3]) if $o_sort_bases and $base{$l[3]} > $base{$l[4]};
    return @l[(0,1,3,4)];
}

sub read_profile_line($$) {
    my ($fh, $tag) = @_;
    my $l = <$fh>;
    if ($o_chr_only) {
        while ($l and $l =~ /^chr/) {
            $l = <$h>;
        }
    }
    return () if ! $l;
    chomp $l;
    #print STDERR "read_profile_line:$tag: '$l'\n";
    return split /\t/, $l;
}

sub read_unselected_line() { return read_profile_line($u, "unselected"); }

sub read_selected_1_line() { return read_profile_line($s1, "selected_1"); }

sub read_selected_2_line() { return read_profile_line($s2, "selected_2"); }

sub sorted_pool_counts($) {
    my $l1 = shift;
    my @d = ( [ $l1->[2], 0 ],   # A
              [ $l1->[3], 1 ],   # C
              [ $l1->[4], 2 ],   # G
              [ $l1->[5], 3 ] ); # T
    @d = sort { $a->[0] <=> $b->[0] } @d;
    ($d[2], $d[3]) = ($d[3], $d[2]) if $o_sort_bases and $d[3]->[1] > $d[2]->[1];
    return @d;
}

# for do_singlepool_test(), return list is:
# 0: reference
# 1: position
# 2: read coverage in unselected pool
# 3: test description string
# 4: unselected probability
# 5: test result string
sub do_unselected_test {
    my ($l, $h_ref, $h_alt) = @_;
    my @d = sorted_pool_counts($l);
    my $cov = $d[2]->[0] + $d[3]->[0];

    my $allele1 = $base[$d[3]->[1]];
    my $allele2 = $base[$d[2]->[1]];
    my $otot = "binom,$allele1/$allele2:$d[3]->[0]/$d[2]->[0]";

    return ($l->[0], $l->[1], $cov, $otot, ".", "zerocov")
        if ($cov == 0); # zero coverage

    my $prob = binomial_test($d[2]->[0], $cov);

    return ($l->[0], $l->[1], $cov, $otot, $prob, "mincov")
        if ($cov < $o_mincov); # insufficient coverage

    if ($o_genotype and defined($h_ref) and defined($h_alt)) {
        # check for incompatible genotype

        carp "unreasonable ref '$h_ref' and alt '$h_alt' alleles"
            if !defined($base{$h_ref}) or !defined($base{$h_alt});

        if (($h_ref ne $allele1 and $h_alt ne $allele2) and 
            ($h_ref ne $allele2 and $h_alt ne $allele1)) {
            $otot .= ",$h_ref/$h_alt";
            return ($l->[0], $l->[1], $cov, $otot, $prob, "genotype");
        }
    }

    return ($l->[0], $l->[1], $cov, $otot, $prob, "allele3")
        if ($d[1]->[0] >= $cov * $o_max_allele_3); # allele 3 coverage too high

    return ($l->[0], $l->[1], $cov, $otot, $prob, "binom_fail")
        if ($prob < $o_unselected_prob); # effectively homozygous site

    return ($l->[0], $l->[1], $cov, $otot, $prob, "binom_pass");
}

sub sorted_twopool_counts($$) {
    my ($l1, $l2) = @_;
    die "inconsistent ref/pos" if $l1->[0] ne $l2->[0] or $l1->[1] != $l2->[1];
    my @d = ( [ $l1->[2], $l2->[2], $l1->[2] + $l2->[2], 0 ],   # A
              [ $l1->[3], $l2->[3], $l1->[3] + $l2->[3], 1 ],   # C
              [ $l1->[4], $l2->[4], $l1->[4] + $l2->[4], 2 ],   # G
              [ $l1->[5], $l2->[5], $l1->[5] + $l2->[5], 3 ] ); # T
    @d = sort { $a->[2] <=> $b->[2] } @d;
    ($d[2], $d[3]) = ($d[3], $d[2]) if $o_sort_bases and $d[2]->[3] < $d[3]->[3];
    return @d;
}
sub rounddot($) {
    my $f = shift;
    return $f eq "." ? $f : sprintf("%.7f", $f);
}
sub multdot($$) {
    my ($p1, $p2) = @_;
    my $ans = ($p1 eq "." or $p2 eq ".") ? "." : ($p1 * $p2);
    return $ans;
}

# for do_twopool_test, return list is:
# 0: reference
# 1: position
# 2: read coverage of allele1 and allele2 in selected 1 pool
# 3: read coverage of allele1 and allele2 in selected 2 pool
# 4: test description string; combination of strings for each pool
# 5: selected 1 and 2 pool unselected test results
# 6: selected 1 and 2 pool binomial probabilities
# 7: selected 1 and 2 pool result of <=> comparisons of allele counts
# 8: twopool_test result string
# 9: multiplied probabilities

sub do_twopool_test($$$$) {
    my ($l1, $l2, $h_ref, $h_alt) = @_;

    my @d      = sorted_twopool_counts($l1, $l2);
    my @utest1 = do_unselected_test($l1, $h_ref, $h_alt);
    my @utest2 = do_unselected_test($l2, $h_ref, $h_alt);

    my $cov1 = $utest1[2];
    my $cov2 = $utest2[2];

    # test; unselected test pool 1 description; unselected test pool 2 description
    my $desc = "hetTwoPool;$utest1[3];$utest2[3]";

    # unselected test pool 1 test result; unselected test pool 2 test result
    my $result12 = "$utest1[5];$utest2[5]";

    # unselected test pool 1 probability; unselected test pool 2 probability;
    my $prob12 = rounddot($utest1[4]) . ";" . rounddot($utest2[4]);

    # two selected pools in allele count opposition
    my $counts1 = $d[3]->[0] <=> $d[2]->[0];  # -1 if 3 < 2, 0 if 3 == 3, 1 if 3 > 2
    my $counts2 = $d[3]->[1] <=> $d[2]->[1];
    my $counts = sprintf("%+d;%+d", $counts1, $counts2);

    # what is going to be the result of this test?
    my $result;
    if ($utest1[5] eq "zerocov" or $utest2[5] eq "zerocov") {
        $result = "zerocov";
    } elsif ($utest1[5] eq "mincov" or $utest2[5] eq "mincov") {
        $result = "mincov";
    } elsif ($utest1[5] eq "genotype" or $utest2[5] eq "genotype") {
        $result = "genotype";
    } elsif ($utest1[5] eq "allele3" or $utest2[5] eq "allele3") {
        $result = "allele3";
    } elsif (! $counts1 or ! $counts2 or $counts1 + $counts2 != 0) {
        $result = "counts";
    } elsif ($utest1[4] > $o_selected_prob or $utest2[4] > $o_selected_prob) {
        $result = $utest1[4] > $o_selected_prob ? "failprob1" : "passprob1";
        $result .= ";";
        $result .= $utest2[4] > $o_selected_prob ? "failprob2" : "passprob2";
    } else {
        $result = "twopool_pass";
    }

    # two selected pools likelihood
    my $prob = multdot($utest1[4], $utest2[4]);

    return ($l1->[0], $l1->[1], $cov1, $cov2, $desc, $result12, $prob12, $counts, $result, $prob);
}

sub do_combined_test($$$$$) {
    my ($u, $l1, $l2, $h_ref, $h_alt) = @_;
    my @result_u = do_unselected_test($u, $h_ref, $h_alt); 
    # result_u has 6 fields, last is test state
    my @result_12 = do_twopool_test($l1, $l2, $h_ref, $h_alt); 
    # result_12 has 10 fields, last is test state
    die "do_combined_test inconsistent ref/pos" if $result_u[0] ne $result_12[0] or $result_u[1] != $result_12[1];
    # combined result has 
    #  0: reference
    #  1: position
    #  2: read coverage in unselected pool
    #  3: test description string, unselected pool
    #  4: unselected probability
    #  5: test result string
    #  6: read coverage in selected 1 pool
    #  7: read coverage in selected 2 pool
    #  8: test description strings from selected 1 and 2 pools
    #  9: unselected test results from selected 1 and 2 pools
    # 10: probabilities from selected 1 and 2 pools
    # 11: <=> comparisons of allele counts from selected 1 and 2 pools
    # 12: twopool test result
    # 13: twopool test multiplied probability
    return (@result_u, @result_12[2..9]);
}

my $config = "Starting hetTwoPool test
fai file                  : $o_fai_file
het sites                 : $o_hets_file
unselected                : $o_unselected_file
selected 1                : $o_selected_1_file
selected 2                : $o_selected_2_file

unselected-test crit prob : $o_unselected_prob
selected-test crit prob   : $o_selected_prob

null-fraction             : $o_null_fraction

mincov                    : $o_mincov
max-allele-3              : $o_max_allele_3

log10-p                   : $o_log10_p
allsites                  : $o_allsites
chr-only                  : $o_chr_only
indels                    : $o_indels
sort-bases                : $o_sort_bases

";
my $outconfig = $config;
$outconfig =~ s/^/#/mg;
print STDERR $config;
print STDOUT $outconfig;

my @h = read_hets_line();  # CHROM POS REF ALT
my @p = read_unselected_line();  # CHROM POS A C G T N
my @s1 = read_selected_1_line();
my @s2 = read_selected_2_line();

sub print_result($) {  # do rounding, etc. for output
    my $r = shift;
    $r->[4] = sprintf("%.5f", $r->[4]) if $r->[4] ne ".";
    $r->[13] = log10($r->[13]) if $r->[13] ne "." and $o_log10_p;
    $r->[13] = sprintf("%.7f", $r->[13]) if $r->[13] ne ".";
    print STDOUT join("\t", @$r), "\n";
}

while (@h and @p and @s1 and @s2) {
    if ($h[0] eq $p[0] and $h[0] eq $s1[0] and $h[0] eq $s2[0]) { # same reference
        if ($h[1] == $p[1] and $h[1] == $s1[1] and $h[1] == $s2[1]) { # same position
            #my @result_un = do_unselected_test(\@p, $h[2], $h[3]); 
            ## result_un has 6 fields, last is test state
            #my @result_two = do_twopool_test(\@s1, \@s2, $h[2], $h[3]); 
            ## result_two has 10 fields, last is test state
            my @result = do_combined_test(\@p, \@s1, \@s2, $h[2], $h[3]);
            @h = read_hets_line();
            @p = read_unselected_line();
            @s1 = read_selected_1_line();
            @s2 = read_selected_2_line();
            next if $result[5] ne "binom_pass" and not $o_allsites;
            # prefix the test states with "het_" to show this was a called het site
            $result[5] = "het_$result[5]";
            $result[12] = "het_$result[12]";
            #print STDOUT join("\t", @result), "\n";
            print_result(\@result);
        } elsif ($h[1] < $p[1]) {
            @h = read_hets_line(); # advance hets file
        } elsif ($p[1] != $s1[1] or $p[1] != $s2[1]) {
            # sync up pools
            do {
                if ($p[1] < $s1[1]) {
                    @p = read_unselected_line();
                } elsif ($s1[1] < $s2[1]) {
                    @s1 = read_selected_1_line();
                } elsif ($s2[1] < $p[1] or $s2[1] < $s1[1]) {
                    @s2 = read_selected_2_line();
                }
            } until ($p[1] == $s1[1] and $p[1] == $s2[1]);
        } elsif ($h[1] > $p[1]) {
            # hets is ahead of the pools, and we know the pools are synced, so
            # we know that this common pool site is not in hets check to see if
            # we can treat it as a null site
            if ($o_null_fraction > 0.0 and rand() < $o_null_fraction) {
                my @result = do_combined_test(\@p, \@s1, \@s2, undef, undef);
                $result[5] = "null_$result[5]";
                $result[12] = "null_$result[12]";
                print STDOUT join("\t", @result), "\n";
                print_result(\@result);
            }
            @p = read_unselected_line(); # advance pools
            @s1 = read_selected_1_line();
            @s2 = read_selected_2_line();
        } else {
            die "unknown condition involving positions";
        }
    } elsif ($REF_ORDER{$h[0]} < $REF_ORDER{$p[0]}) {
        @h = read_hets_line(); # advance $h because wrong ref sequence
    } elsif ($REF_ORDER{$h[0]} > $REF_ORDER{$p[0]}) {
        if ($o_null_fraction > 0.0 and !$INDELS{"$p[0]:$p[1]"} and rand() < $o_null_fraction) {
            my @result = do_combined_test(\@p, \@s1, \@s2, undef, undef);
            $result[5] = "null_$result[5]";
            $result[12] = "null_$result[12]";
            #print STDOUT join("\t", @result), "\n";
            print_result(\@result);
        }
        @p = read_unselected_line(); # advance pools because wrong ref sequence
        @s1 = read_selected_1_line();
        @s2 = read_selected_2_line();
    } else {
        die "unknown condition involving reference order";
    }
}


print STDERR "Skipped $N_indels indels\n";

print STDERR "Finished hetTwoPool test\n";

