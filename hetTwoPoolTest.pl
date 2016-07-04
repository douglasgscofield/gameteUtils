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
my $o_null_fraction = 0.001;
my $o_mincov = 8; # 10
my $o_genotype = 1;
my $o_max_allele_3 = 0.1;

my $o_pool_freq_prob = 0.10;
my $o_pool_test_type = "binom";

my $o_log10_p = 0;
my $o_allsites = 1;

my $n_homozygous = 0;
my $o_max_p_val = 0.05;

my $o_help;
my $current_reference = ""; # the name of the current reference sequence
my $N_references = 0; # the number of reference sequences seen
my $N_coordinates = 0; # the total number of coordinates seen
my $N_indels = 0; # total number of indels skipped
my $o_indels;
my %INDELS;

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
$o_indels = $o_hets_file . "_indels.txt";
open (my $indels, ">", $o_indels) or die "cannot open indels output file '$o_indels': $!"; 

sub read_hets_line() {
    READ_LINE:
    my $l = <$h>;
    while ($l and $l =~ /^#/) {
        $l = <$h>;
    }
    return () if ! $l;
    chomp $l;
    #print STDERR "read_hets_line: '$l'\n";
    my @l = split /\t/, $l;
    if (length($l[3]) > 1 or length($l[4]) > 1) {  # ref or alt not a single base
        ++$N_indels;
        ++$INDELS{"$l[0]:$l[1]"};
        print { $indels } "INDEL\t", join("\t", @l[(0,1,3,4)]), "\n";
        goto READ_LINE;
    }
    # return CHROM POS REF ALT
    return @l[(0,1,3,4)];
}

sub read_pool_line() {
    my $l = <$p>;
    return () if ! $l;
    chomp $l;
    #print STDERR "read_pool_line: '$l'\n";
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
    return sort { $a->[0] <=> $b->[0] } @d;
}

# for do_unselected_test(), output is:
# 0: reference
# 1: position
# 2: read coverage
# 3: test description string
# 4: unselected probability
# 5: test result string
sub do_unselected_test {
    my ($l, $href, $halt) = @_;
    my @d = sorted_pool_counts($l);
    my $cov = $d[2]->[0] + $d[3]->[0];

    my $major = $base[$d[3]->[1]];  # major allele
    my $minor = $base[$d[2]->[1]];  # minor allele
    my $otot = "binom,$major/$minor:$d[3]->[0]/$d[2]->[0]";

    return ($l->[0], $l->[1], $cov, $otot, ".", "zerocov") if ($cov == 0); # zero coverage

    my $prob = sprintf("%.5f", binomial_test($d[2]->[0], $cov));

    return ($l->[0], $l->[1], $cov, $otot, $prob, "mincov") if ($cov < $o_mincov); # insufficient coverage

    if ($o_genotype and defined($href) and defined($halt)) { # check for incompatible genotype
        carp "unreasonable ref '$href' and alt '$halt' alleles" if !defined($base{$href}) or !defined($base{$halt});
        if (($href ne $major and $halt ne $minor) and ($href ne $minor and $halt ne $major)) {
            $otot .= ",$href/$halt";
            return ($l->[0], $l->[1], $cov, $otot, $prob, "genotype");
        }
    }
    return ($l->[0], $l->[1], $cov, $otot, $prob, "allele3") if ($d[1]->[0] >= $cov * $o_max_allele_3); # if coverage of allele 3 too high
    return ($l->[0], $l->[1], $cov, $otot, $prob, "binom_fail") if ($prob < $o_unselected_prob); # effectively homozygous site
    return ($l->[0], $l->[1], $cov, $otot, $prob, "binom_pass");
}

print STDERR "Starting unselected test with '$o_hets_file' and '$o_pool_file'\n";

my @h = read_hets_line();  # CHROM POS REF ALT
my @p = read_pool_line();  # CHROM POS A C G T N

while (@h and @p) {
    if ($h[0] eq $p[0]) { # same reference
        if ($h[1] == $p[1]) { # same position
            my @result = do_unselected_test(\@p, $h[2], $h[3]); 
            # result has 6 fields, last is test state
            @h = read_hets_line();
            @p = read_pool_line();
            next if $result[5] ne "binom_pass" and not $o_allsites;
            # prefix the test state with "het_" to show this was a called het site
            $result[5] = "het_$result[5]";
            print STDOUT join("\t", @result), "\n";
        } elsif ($h[1] < $p[1]) {
            @h = read_hets_line(); # advance hets file
        } elsif ($h[1] > $p[1]) {
            # hets is ahead of the pool, so we know that this pool site is not in hets
            # check to see if we can treat it as a null site
            if ($o_null_fraction > 0.0 and rand() < $o_null_fraction) {
                my @result = do_unselected_test(\@p, undef, undef); 
                $result[5] = "null_$result[5]";
                print STDOUT join("\t", @result), "\n";
            }
            @p = read_pool_line(); # now advance pool file
        } else {
            die "unknown condition involving positions";
        }
    } elsif ($REF_ORDER{$h[0]} < $REF_ORDER{$p[0]}) {
        @h = read_hets_line(); # advance $h because wrong ref sequence
    } elsif ($REF_ORDER{$h[0]} > $REF_ORDER{$p[0]}) {
        if ($o_null_fraction > 0.0 and !$INDELS{"$p[0]:$p[1]"} and rand() < $o_null_fraction) {
            my @result = do_unselected_test(\@p); 
            $result[5] = "null_$result[5]";
            print STDOUT join("\t", @result), "\n";
        }
        @p = read_pool_line(); # advance $p because wrong ref sequence
    } else {
        die "unknown condition involving reference order";
    }
}


print STDERR "Skipped $N_indels indels\n";

print STDERR "Finished unselected test\n";

#!/usr/bin/env perl

# Copyright (c) 2012,2015 Douglas G. Scofield, Uppsala University
# douglasgscofield@gmail.com
#
# No warranty is implied or assumed by this code.  Please send bugs, suggestions etc.
#
use strict;
use warnings;
use POSIX qw/isdigit log10/;
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
my $o_total_freq_prob = 0.01;
my $o_mincov = 16; # 10
my $o_minalt = 2;
my $o_max_allele_3 = 0.1;

my $o_pool_freq_prob = 0.10;
my $o_pool_test_type = "binom";

my $o_log10_p = 1;
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

    --mincov INT            Minimum total read coverage to consider a site
    --minalt INT            Minimum reads for alternate allele to consider a site

  ALLELE FREQUENCY TESTS

    --pool-test-type STRING Type of test for pool frequencies: chisq or binom [default $o_pool_test_type]

    --pool-freq-prob FLOAT  Probability threshold below which a pool is considered to violate 1:1
                            [default $o_pool_freq_prob]

    --max-p-val FLOAT       Maximum P value to report, ignored for 'binom'

  OTHER OPTIONS

    --log10-p               For the two-pool test, output the log10-d P value [default $o_log10_p]
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

    my $otot = "binom,$base[$d[3]->[3]]/$base[$d[2]->[3]]:$d[3]->[2]/$d[2]->[2]";

    return ($l1->[0], $l1->[1], $otot, ".", ".", ".", ".") if $cov < $o_mincov;  # insufficient coverage

    my $totprob = sprintf("%.3f", binomial_test($d[2]->[2], $cov));

    if ($totprob < $o_total_freq_prob) {  # effectively homozygous site
        return ($l1->[0], $l1->[1], $otot, $totprob, ".", ".", ".");
    }

    if ($d[1]->[2] >= $cov * $o_max_allele_3) {  # if coverage of allele 3 too high
        return ($l1->[0], $l1->[1], $otot, $totprob, "-1", "-1", "-1");
    }

    my $opool1 = "$d[3]->[0]/$d[2]->[0]";  # counts sample 1
    my $opool2 = "$d[3]->[1]/$d[2]->[1]";  # counts sample 2

    my $pool1cov = $d[3]->[0] + $d[2]->[0];
    my $pool2cov = $d[3]->[1] + $d[2]->[1];
    my $pool1prob = binomial_test($d[2]->[0], $pool1cov);
    my $pool2prob = binomial_test($d[2]->[1], $pool2cov);
    $opool1 = "$opool1:" . sprintf("%.3f", $pool1prob);
    $opool2 = "$opool2:" . sprintf("%.3f", $pool2prob);
    my $twopoolprob = $pool1prob * $pool2prob;
    my $onepoolprob = sprintf("%.6f", $twopoolprob);
    $twopoolprob = log10($twopoolprob) if $o_log10_p;
    $twopoolprob = sprintf("%.6f", $twopoolprob);

    # now check to see if we can perform a two-pool test
    my $counts1 = $d[3]->[0] <=> $d[2]->[0];  # -1 of 3 < 2, 0 if 3 == 3, 1 if 3 > 2
    my $counts2 = $d[3]->[1] <=> $d[2]->[1];
    # for a two-pool test, the counts must differ and be of opposite sign and the 
    # test probability must be significant for both
    if (! $counts1 or ! $counts2 or $counts1 + $counts2 != 0 or
        $pool1prob > $o_pool_freq_prob or $pool2prob > $o_pool_freq_prob) {
        $twopoolprob = ".";
    }
    # for a one-pool test, one must pass and one must fail
    if (not $pool1prob > $o_pool_freq_prob xor $pool2prob > $o_pool_freq_prob) {
        $onepoolprob = ".";
    }
    return ($l1->[0], $l1->[1], $otot, $totprob, "$opool1,$opool2", $twopoolprob, $onepoolprob);
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

    my $o = "$base[$d[3]->[3]]/$base[$d[2]->[3]]," .  # allele 1 / allele 2
            "chisq,$d[3]->[0]/$d[2]->[0],$d[3]->[1]/$d[2]->[1]";  # counts sample 2 , counts sample 2

    return ($l1->[0], $l1->[1], $o, ".", ".") if $cov < $o_mincov;  # insufficient coverage
    return ($l1->[0], $l1->[1], $o, "-1", "-1") if $d[1]->[2] >= $cov * $o_max_allele_3;  # there is an allele 3
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
    $prob = sprintf("%.6f", $prob);
    return ($l1->[0], $l1->[1], $o, $chi2, $prob);
}

print STDERR "Starting single parent test with '$o_sample1' and '$o_sample2'\n";

my @l1 = read_line($s1);
my @l2 = read_line($s2);

while (@l1 and @l2) {
    if ($l1[0] eq $l2[0]) { # same reference
        if ($l1[1] == $l2[1]) { # same position
            my @result;
            if ($o_pool_test_type eq 'binom') {
                # result has 7 fields
                @result = do_binom_test(\@l1, \@l2); 
                @l1 = read_line($s1);
                @l2 = read_line($s2);
                if (not $o_allsites) {
                    next if $result[3] eq "." or $result[3] < $o_total_freq_prob;
                    next if $result[4] eq "." or $result[4] eq "-1";
                    next if $result[5] eq "." and $result[6] eq ".";
                }
            } elsif ($o_pool_test_type eq 'chisq') {
                # result has 5 fields
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

print STDERR "Finished single parent test\n";
