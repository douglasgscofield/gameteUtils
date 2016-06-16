#!/usr/bin/env perl

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
# http://stackoverflow.com/questions/21204733/a-better-chi-square-test-for-perl
use Statistics::Distributions qw/ chisqrprob /;
use FindBin;
use lib $FindBin::RealBin;  # add script directory to @INC to find BinomialTest
use BinomialTest qw/ binomial_test /;

my $o_fai_file;
my $o_hets_file;
my $o_pool_file;
my $o_unselected_test = 1;
my $o_unselected_prob = 0.01;
my $o_null_fraction = 0.0;
my $o_mincov = 16; # 10
my $o_genotype = 1;
my $o_max_allele_3 = 0.1;

my $o_pool_freq_prob = 0.10;
my $o_pool_test_type = "binom";

my $o_log10_p = 0;
my $o_allsites = 0;

my $n_homozygous = 0;
my $o_max_p_val = 0.05;

my $o_help;
my $current_reference = ""; # the name of the current reference sequence
my $N_references = 0; # the number of reference sequences seen
my $N_coordinates = 0; # the total number of coordinates seen

my $usage = "
NAME

  unselectedTest.pl - test list of heterozygous sites for deviations from 1:1 in an 'unselected' pool


SYNOPSIS

    $0 --fai reference.fa.fai --hets h.vcf --pool pool.pro


OPTIONS

    --fai FILE              Fasta index file for reference, required to specify the order of
                            sequences in the VCF and profile files

    --hets FILE             VCF-format file specifying heterozygous sites. Use care in producing
                            this file; for simplicity, this script only uses presence-absence 
                            positional information from the file to determine whether a site is
                            heterozygous, it explicitly does *not* use any sort of genotype
                            quality, PASS/FAIL, etc.  Any sort of filtering you would wish to
                            provide must be done before this file is read by this script.  One
                            advantage of the simplification is that the file needn't be in a
                            fully compliant VCF format, so long as the first five tab-separated
                            columns of lines that do not begin with '#' contain the same
                            information as is in the '#CHROM POS ID REF ALT' columns of a VCF.
                            The ID isn't used and may be '.', but REF and ALT are used to check
                            the alleles in the unselected pool.

    --pool FILE             Pool profile.  If this is produced from a pileup, it will contain
                            all 

    --unselected-test FLOAT Test for deviation of pool allele frequency away from 1:1, using a
                            two-sided binomial test with expected frequency 0.5.
                            Only sites which fail this test (and thus have total ratio of
                            1:1) will be given the differential frequency test via chi-squared.
                            This test is performed by default, so it is unnecessary to specify
                            this option, unless you want to change FLOAT, which is the probability
                            threshold for the test; sites with binomial probabilities <= FLOAT
                            are considered to fail this test and show allele frequency skew.  The
                            5th column of output indicates the test disposition, and is prefixed
                            with 'het_' for each site that is a called heterozygote.
                            [default $o_unselected_prob]

    --null-fraction FLOAT   For fraction FLOAT sites that are *not* called heterozygotes, perform
                            the unselected test as above.  Output This fraction is also subject to the
                            --mincov and --max-allele-3 settings.  The 5th column of output indicates
                            the test disposition, and is prefixed with 'null_' for each site tested as
                            a result of this option.  [default $o_null_fraction]

    --no-unselected-test    Do NOT perform the unselected test

    --mincov INT            Minimum read coverage to consider a site [default $o_mincov]

    --max-allele-3 FLOAT    Maximum accepted frequency count for a 3rd allele [default $o_max_allele_3]

  OTHER OPTIONS

    --log10-p               Output log10-d P values [default $o_log10_p]  NOT FUNCTIONAL

    --allsites              Show test result for all heterozygous sites, not just those that pass the test

    --help, -?              Help message

OUTPUT

Five tab-separated columns:

    reference, position, test description, test probability, test result

'Test result' contains one of the following, prefixed with 'het_' for a
heterozygous site and 'null_' for a non-heterozygous site selected via
--null-fraction.

    zerocov    : no pool coverage at this site; test probability is reported as '.'
    mincov     : pool coverage is below the --mincov value
    genotype   : genotype incompatibility (REF and ALT do not match top pool alleles)
    allele3    : the coverage frequency of a 3rd allele is higher than --max-allele-3
    binom_fail : the binomial test is less than --unselected-prob
    binom_pass : the binomial test is >= --unselected-prob

";

sub print_usage_and_exit {
    my $msg = join(" ", @_);
    print "$msg\n" if $msg;
    print $usage;
    exit 0;
}

GetOptions(
    "fai=s"               => \$o_fai_file,
    "hets=s"              => \$o_hets_file,
    "pool=s"              => \$o_pool_file,
    "unselected-test=f"   => \$o_unselected_prob,
    "null-fraction=f"     => \$o_null_fraction,
    "no-unselected-test"  => sub { $o_unselected_test = 0 },
    "mincov=i"            => \$o_mincov,
    "max-allele-3=f"      => \$o_max_allele_3,
    "log10-p"             => \$o_log10_p,
    "allsites"            => \$o_allsites,
    "help|?"              => \$o_help,
) or print_usage_and_exit();

print_usage_and_exit() if $o_help or not $o_fai_file or not $o_hets_file or not $o_pool_file;

# fill reference sequence order hash from fai file
my %REF_ORDER;
my $ref_index = 0;
open (my $ref, "<", $o_fai_file) or die "cannot open Fasta index file '$o_fai_file': $!";
while (<$ref>) {
    my @l = split /\t/;
    $REF_ORDER{$l[0]} = ++$ref_index if not exists $REF_ORDER{$l[0]};
}
print STDERR "Found ".scalar(keys(%REF_ORDER))." reference sequences in $o_fai_file\n";

open (my $h, "<", $o_hets_file) or die "cannot open hets file '$o_hets_file': $!";
open (my $p, "<", $o_pool_file) or die "cannot open pool profile '$o_pool_file': $!";

sub read_hets_line() {
    my $l = <$h>;
    while ($l and $l =~ /^#/) {
        $l = <$h>;
    }
    return () if ! $l;
    chomp $l;
    print STDERR "read_hets_line: '$l'\n";
    my @l = split /\t/, $l;
    # return CHROM POS REF ALT
    return @l[(0,1,3,4)];
}

sub read_pool_line() {
    my $l = <$p>;
    return () if ! $l;
    chomp $l;
    print STDERR "read_pool_line: '$l'\n";
    return split /\t/, $l;
}

my @base = qw/ A C G T /;
my %base = map { $_ => 1 } @base;

sub sorted_pool_counts($) {
    my $l1 = shift;
    my @d = ( [ $l1->[2], 0 ],   # A
              [ $l1->[3], 1 ],   # C
              [ $l1->[4], 2 ],   # G
              [ $l1->[5], 3 ] ); # T
    return sort { $a->[2] <=> $b->[2] } @d;
}

# for do_unselected_test(), output is:
# 0: reference
# 1: position
# 2: test description string
# 3: unselected probability
# 4: test result string
sub do_unselected_test {
    my ($l, $href, $halt) = @_;
    my @d = sorted_pool_counts($l);
    my $cov = $d[2]->[0] + $d[3]->[0];

    my $major = $base[$d[3]->[1]];  # major allele
    my $minor = $base[$d[2]->[1]];  # minor allele
    my $otot = "binom,$major/$minor:$d[3]->[0]/$d[2]->[0]";

    return ($l->[0], $l->[1], $otot, ".", "zerocov") if ($cov == 0); # zero coverage

    my $prob = sprintf("%.5f", binomial_test($d[2]->[0], $cov));

    return ($l->[0], $l->[1], $otot, $prob, "mincov") if ($cov < $o_mincov); # insufficient coverage

    if ($o_genotype and defined($href) and defined($halt)) { # check for incompatible genotype
        carp "unreasonable ref '$href' and alt '$halt' alleles" if !defined($base{$href}) or !defined($base{$halt});
        if (($href eq $major and $halt eq $minor) or ($href eq $minor and $halt eq $major)) {
            return ($l->[0], $l->[1], $otot, $prob, "genotype");
        }
    }
    return ($l->[0], $l->[1], $otot, $prob, "allele3") if ($d[1]->[0] >= $cov * $o_max_allele_3); # if coverage of allele 3 too high
    return ($l->[0], $l->[1], $otot, $prob, "binom_fail") if ($prob < $o_unselected_prob); # effectively homozygous site
    return ($l->[0], $l->[1], $otot, $prob, "binom_pass");
}

print STDERR "Starting unselected test with '$o_hets_file' and '$o_pool_file'\n";

my @h = read_hets_line();  # CHROM POS REF ALT
my @p = read_pool_line();  # CHROM POS A C G T N

while (@h and @p) {
    if ($h[0] eq $p[0]) { # same reference
        if ($h[1] == $p[1]) { # same position
            my @result = do_unselected_test(\@p, $h[2], $h[3]); 
            # result has 5 fields, last is test state
            @h = read_hets_line();
            @p = read_pool_line();
            next if $result[4] ne "binom_pass" and not $o_allsites;
            # prefix the test state with "het_" to show this was a called het site
            $result[4] = "het_$result[4]";
            print STDOUT join("\t", @result), "\n";
        } elsif ($h[1] < $p[1]) {
            @h = read_hets_line(); # advance hets file
        } elsif ($h[1] > $p[1]) {
            # hets is ahead of the pool, so we know that this pool site is not in hets
            # check to see if we can treat it as a null site
            if ($o_null_fraction > 0.0 and rand() < $o_null_fraction) {
                my @result = do_unselected_test(\@p, undef, undef); 
                $result[4] = "null_$result[4]";
                print STDOUT join("\t", @result), "\n";
            }
            @p = read_pool_line(); # now advance pool file
        } else {
            die "unknown condition involving positions";
        }
    } elsif ($REF_ORDER{$h[0]} < $REF_ORDER{$p[0]}) {
        @h = read_hets_line(); # advance $h because wrong ref sequence
    } elsif ($REF_ORDER{$h[0]} > $REF_ORDER{$p[0]}) {
        if ($o_null_fraction > 0.0 and rand() < $o_null_fraction) {
            my @result = do_unselected_test(\@p); 
            $result[4] = "null_$result[4]";
            print STDOUT join("\t", @result), "\n";
        }
        @p = read_pool_line(); # advance $p because wrong ref sequence
    } else {
        die "unknown condition involving reference order";
    }
}

print STDERR "Finished unselected test\n";

