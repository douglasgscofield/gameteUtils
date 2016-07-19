#!/usr/bin/env perl

# TODO: two-tailed unselected probability ?
# TODO: probability output: was %0.5f but adjust during output??
# TODO: handle non-0/1 heterozygotes
# DONE: create module to reduce lots of duplicated code between here and hetTwoPoolTest.pl
# DONE: add usage of hetsites-reads profile
# DONE: add --indels-annotate

# Copyright (c) 2016 Douglas G. Scofield, Uppsala University
# douglas.scofield@ebc.uu.se, douglasgscofield@gmail.com
#
# No warranty is implied or assumed by this code.  Please send bugs, suggestions etc.
#
use strict;
use warnings;
use Carp;
use Getopt::Long;
use List::Util;
use FindBin;
use lib $FindBin::RealBin;  # add script directory to @INC to find BinomialTest and/or GameteUtils
use GameteUtils qw/ fill_ref_index open_possibly_gzipped read_hetsites_line read_profile_line sorted_pool_counts rounddot multdot log10dot do_unselected_test /;
# multdot() not actually used in this file ...

# 'our' variables below are also used within the GameteUtils package
my $o_fai_file;
my $o_hetsites_file;
my $o_hetprofile_file;
my $o_unselected_file;
our $o_unselected_prob = 0.01;  # --unselected-test FLOAT
our $o_mincov = 8; # 10   # --mincov INT
our $o_genotype = 1;  # check for incompatible genotype, always done
our $o_max_allele_3 = 0.1;  # --max-allele-3 FLOAT
my $o_null_fraction = 0.000;

my $o_log10_p = 1;
our $o_sort_bases = 1;
our $o_hets_genotype = 1;  # GameteUtils::do_unselected_test() includes genotype columns
my $o_allsites = 1;

my $o_help;
my $current_reference = ""; # the name of the current reference sequence
my $N_references = 0; # the number of reference sequences seen
my $N_coordinates = 0; # the total number of coordinates seen
our $o_digits = 7; # number of digits to round output
our $N_indels = 0; # total number of indels skipped
our $o_indels = 1; # produce file of indels
our $o_indels_annotate = 1;  # annotate indels with quality, INFO and genotype columns from VCF
our $indels;  # $indels is file handle for indels file
our %INDELS;

my $usage_short = "
    $0 --fai reference.fa.fai --hetsites h.vcf [ --hetprofile het.pro ] --unselected u.pro 

        --fai FILE
        --hetsites FILE
        --hetprofile FILE
        --unselected FILE
        --unselected-test FLOAT  [default $o_unselected_prob]
        --mincov INT             [default $o_mincov]
        --max-allele-3 FLOAT     [default $o_max_allele_3]
        --null-fraction FLOAT    [default $o_null_fraction]
        --log10-p                [default $o_log10_p]
        --indels                 [default $o_indels]
        --no-indels
        --indels-annotate        [default $o_indels_annotate]
        --no-indels-annotate
        --allsites               [default $o_allsites]
        --sort-bases             [default $o_sort_bases]
        --no-sort-bases
        --help, -?               More detailed help message

";

my $usage = "
NAME

  unselectedTest.pl - test list of heterozygous sites for deviations from 1:1 in an 'unselected' pool


SYNOPSIS

    $0 --fai reference.fa.fai --hetsites h.vcf [ --hetprofile het.pro ] --unselected unselected.pro


OPTIONS

    --fai FILE    Fasta index file for reference, required to specify the order of
                  sequences in the VCF and profile files

    --hetsites FILE
                  VCF-format file specifying heterozygous sites. Use care in
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

    --hetprofile FILE
                  Profile base content for reads underlying hetsites VCF.  This
                  file is OPTIONAL.  If it is included, then an unselected test
                  is performed on the content of this file too, and those results
                  columns precede those for the unselected test performed on the
                  profile specified with --unselected.
     

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

    --mincov INT  Minimum read coverage to consider a site, applied to each of
                  the unselected and selected pools [default $o_mincov].

    --max-allele-3 FLOAT
                  Maximum accepted frequency count for a 3rd allele, applied to
                  each of the unselected and selected pools [default $o_max_allele_3]

    --null-fraction FLOAT
                  For fraction FLOAT sites that are *not* called heterozygotes,
                  perform the unselected and selected tests as above.  Output
                  This fraction is also subject to the --mincov and
                  --max-allele-3 settings.  The 5th column of output indicates
                  the test disposition, and is prefixed with 'null_' for each
                  site tested as a result of this option.  [default
                  $o_null_fraction]

  OTHER OPTIONS

    --log10-p     Output log10-d P values for unselected tests [default $o_log10_p]

    --indels
    --no-indels   Write (or don't) apparent indels to an indels output file.  If the
                  output file cannot be created, this option is disabled [default $o_indels].

    --indels-annotate
    --no-indels-annotate
                  Annotate indels written above with quality, INFO and genotype
                  columns from the input VCF [default $o_indels_annotate]

    --allsites    Show test result for all heterozygous sites, not just those
                  for which the unselected pool passes the unselected test

    --sort-bases
    --no-sort-bases
                  In output, sort, or do not sort, bases in lexical order (A,C,G,T)
                  The default is ".($o_sort_bases?"":"not")." to sort

    --help, -?    Help message

OUTPUT

Eight tab-separated columns:

  1 reference
  2 position
  3 genotype from the hets VCF file
  4 genotype quality from the hets VCF file
  5 unselected pool read coverage
  6 unselected pool test description string
  7 unselected pool probability
  8 unselected pool test result

'Test result' contains one of the following, prefixed with 'het_' for a
heterozygous site and 'null_' for a non-heterozygous site selected via
--null-fraction.

  zerocov    : no unselected coverage at this site; test probability is reported as '.'
  mincov     : unselected coverage is below the --mincov value
  genotype   : genotype incompatibility (REF and ALT do not match top unselected alleles)
  allele3    : the coverage frequency of a 3rd allele is higher than --max-allele-3
  binom_fail : the binomial test is less than --unselected-prob
  binom_pass : the binomial test is >= --unselected-prob

If a --hetprofile file has been provided, then the output contains 12
tab-separated columns, with the unselected test columns (described above) moved
to columns 7-10, and columns 3-6 now containing the results of an unselected
test applied to the contents of the --hetprofile file.  No special prefix is
added to the test result in these columns.  The test result for the unselected
test still contains the 'het_' or 'null_' prefixes as described above.

  1 reference
  2 position
  3 genotype from the hets VCF file
  4 genotype quality from the hets VCF file
  5 hetprofile read coverage
  6 hetprofile test description string
  7 hetprofile probability
  8 hetprofile test result
  9 unselected pool read coverage
 10 unselected pool test description string
 11 unselected pool probability
 12 unselected pool test result


";

sub print_usage_and_exit {
    my $x = shift;
    my $msg = join(" ", @_);
    print "$msg\n" if $msg;
    print $usage;
    $x ||= 0;
    exit $x;
}

if (scalar(@ARGV) == 0) {
    print $usage_short;
    exit 1;
}

GetOptions(
    "fai=s"               => \$o_fai_file,
    "hetsites=s"          => \$o_hetsites_file,
    "hetprofile=s"        => \$o_hetprofile_file,
    "unselected=s"        => \$o_unselected_file,
    "unselected-test=f"   => \$o_unselected_prob,
    "null-fraction=f"     => \$o_null_fraction,
    "mincov=i"            => \$o_mincov,
    "max-allele-3=f"      => \$o_max_allele_3,
    "log10-p"             => \$o_log10_p,
    "indels"              => \$o_indels,
    "no-indels"           => sub { $o_indels = 0 },
    "indels-annotate"     => \$o_indels_annotate,
    "no-indels-annotate"  => sub { $o_indels_annotate = 0 },
    "allsites"            => \$o_allsites,
    "help|?"              => \$o_help,
) or print_usage_and_exit(1);

print_usage_and_exit(1) if $o_help 
                          or not $o_fai_file
                          or not $o_hetsites_file
                          or not $o_unselected_file;

$o_indels ||= $o_indels_annotate;

# fill reference sequence order hash from fai file
my %REF_ORDER = fill_ref_index($o_fai_file);
print STDERR "Found ".scalar(keys(%REF_ORDER))." reference sequences in $o_fai_file\n";

my $h = open_possibly_gzipped($o_hetsites_file);

if ($o_indels) {
    $o_indels = $o_hetsites_file . "_unsel_indels.txt";
    if (! open ($indels, ">", $o_indels)) {
        carp "cannot open indels output file '$o_indels', disabling option\n"; 
        $o_indels = 0;
        $o_indels_annotate = 0;
    }
}

my $u = open_possibly_gzipped($o_unselected_file);
my $hp;
if ($o_hetprofile_file) {
    $hp = open_possibly_gzipped($o_hetprofile_file);
}

sub read_hetprofile_line() { return read_profile_line($hp, "hetprofile"); }

sub read_unselected_line() { return read_profile_line($u, "unselected"); }

sub do_combined_test($$$$$) {  # called if --hetprofile file given
    my ($hp, $u, $h_ref, $h_alt, $h_qual) = @_;
    my @result_hp = do_unselected_test($hp, $h_ref, $h_alt, $h_qual); 
    my @result_u = do_unselected_test($u, $h_ref, $h_alt, $h_qual); 
    croak "do_combined_test inconsistent ref/pos" if $result_hp[0] ne $result_u[0] or $result_hp[1] != $result_u[1];
    # combined result has 
    #  1 reference
    #  2 position
    #  3 genotype from the hets VCF file
    #  4 genotype quality from the hets VCF file
    #  5 hetprofile read coverage
    #  6 hetprofile test description string
    #  7 hetprofile probability
    #  8 hetprofile test result
    #  9 unselected pool read coverage
    # 10 unselected pool test description string
    # 11 unselected pool probability
    # 12 unselected pool test result
    return (@result_hp, @result_u[4..7]);
}


my $config = "Starting unselected test
fai file                  : $o_fai_file
het sites                 : $o_hetsites_file
het profile               : $o_hetprofile_file
unselected                : $o_unselected_file

unselected-test crit prob : $o_unselected_prob

null-fraction             : $o_null_fraction

mincov                    : $o_mincov
max-allele-3              : $o_max_allele_3

log10-p                   : $o_log10_p
allsites                  : $o_allsites
indels                    : $o_indels
indels-annotate           : $o_indels_annotate
sort-bases                : $o_sort_bases
hets-genotype             : $o_hets_genotype

";
my $outconfig = $config;
$outconfig =~ s/^/#/mg;
print STDERR $config;
print STDOUT $outconfig;

my @h = read_hetsites_line($h);  # CHROM POS REF ALT QUAL
my @p = read_unselected_line(); # CHROM POS A C G T N
my @hp; @hp = read_hetprofile_line() if $o_hetprofile_file;

sub print_result($) {  # do rounding, etc. for output
    my $r = shift;
    # col 7 is always a probability even though its meaning shifts with --hetprofile
    $r->[6] = log10dot($r->[6]) if $o_log10_p;
    $r->[6] = rounddot($r->[6]);
    if ($o_hetprofile_file) {
        $r->[10] = log10dot($r->[10]) if $o_log10_p;
        $r->[10] = rounddot($r->[10]);
    }
    print STDOUT join("\t", @$r), "\n";
}

if ($o_hetprofile_file) {  # --hetprofile file given, 12-column output
    while (@h and @p and @hp) {
        if ($h[0] eq $p[0] and $h[0] eq $hp[0]) { # same reference
            if ($h[1] == $p[1] and $h[1] == $hp[1]) { # same position
                my @result = do_combined_test(\@hp, \@p, $h[2], $h[3], $h[4]);
                @h = read_hetsites_line($h);
                @p = read_unselected_line();
                @hp = read_hetprofile_line();
                next if $result[11] ne "binom_pass" and not $o_allsites;
                # prefix the unselected test states with "het_" to show this was a called het site
                $result[11] = "het_$result[11]";
                #print STDOUT join("\t", @result), "\n";
                print_result(\@result);
            } elsif ($h[1] < $p[1]) {
                @h = read_hetsites_line($h); # advance hetsites file
            } elsif ($p[1] != $hp[1]) {
                # sync up hetprofile and unselected pool
                do {
                    if    ($p[1] < $hp[1]) { @p  = read_unselected_line(); }
                    elsif ($p[1] > $hp[1]) { @hp = read_hetprofile_line(); }
                } until ($p[1] == $hp[1]);
            } elsif ($h[1] > $p[1]) {
                # hetsites is ahead of hetprofile and the unselected pools, and
                # we know they are synced, so we know that this common pool
                # site is not in hetsites. check to see if we can treat it as a
                # null site
                if ($o_null_fraction > 0.0 and !$INDELS{"$p[0]:$p[1]"} and rand() < $o_null_fraction) {
                    my @result = do_combined_test(\@hp, \@p, undef, undef, undef);
                    $result[11] = "null_$result[11]";
                    print_result(\@result);
                }
                @p = read_unselected_line(); # advance pools
                @hp = read_hetprofile_line();
            } else {
                croak "unknown condition involving positions";
            }
        } elsif ($REF_ORDER{$h[0]} < $REF_ORDER{$p[0]}) {
            @h = read_hetsites_line($h); # advance $h because wrong ref sequence
        } elsif ($REF_ORDER{$p[0]} != $REF_ORDER{$hp[0]}) {
            # sync up hetprofile and unselected pool
            do {
                if    ($REF_ORDER{$p[0]} < $REF_ORDER{$hp[0]}) { @p  = read_unselected_line(); }
                elsif ($REF_ORDER{$p[0]} > $REF_ORDER{$hp[0]}) { @hp = read_hetprofile_line(); }
            } until ($REF_ORDER{$p[0]} == $REF_ORDER{$hp[0]});
        } elsif ($REF_ORDER{$h[0]} > $REF_ORDER{$p[0]}) {
            # hetsites is ahead of hetprofile and the unselected pools, and
            # we know they are synced, so we know that this common pool
            # site is not in hetsites. check to see if we can treat it as a
            # null site
            if ($o_null_fraction > 0.0 and !$INDELS{"$p[0]:$p[1]"} and rand() < $o_null_fraction) {
                my @result = do_combined_test(\@hp, \@p, undef, undef, undef);
                $result[11] = "null_$result[11]";
                print_result(\@result);
            }
            @p = read_unselected_line(); # advance pools
            @hp = read_hetprofile_line();
        } else {
            croak "unknown condition involving reference order";
        }
    }
} else {  # no --hetprofile file, 8-column output
    while (@h and @p) {
        if ($h[0] eq $p[0]) { # same reference
            if ($h[1] == $p[1]) { # same position
                my @result = do_unselected_test(\@p, $h[2], $h[3], $h[4]); 
                # result has 8 fields, 8th is test state
                @h = read_hetsites_line($h);
                @p = read_unselected_line();
                next if $result[7] ne "binom_pass" and not $o_allsites;
                # prefix the test state with "het_" to show this was a called het site
                $result[7] = "het_$result[7]";
                print_result(\@result);
            } elsif ($h[1] < $p[1]) {
                @h = read_hetsites_line($h); # advance hetsites file
            } elsif ($h[1] > $p[1]) {
                # hetsites is ahead of the pool, so we know that this pool site is not in hetsites
                # check to see if we can treat it as a null site
                if ($o_null_fraction > 0.0 and !$INDELS{"$p[0]:$p[1]"} and rand() < $o_null_fraction) {
                    my @result = do_unselected_test(\@p, undef, undef); 
                    $result[7] = "null_$result[7]";
                    print_result(\@result);
                }
                @p = read_unselected_line(); # now advance unselected file
            } else {
                croak "unknown condition involving positions";
            }
        } elsif ($REF_ORDER{$h[0]} < $REF_ORDER{$p[0]}) {
            @h = read_hetsites_line($h); # advance $h because wrong ref sequence
        } elsif ($REF_ORDER{$h[0]} > $REF_ORDER{$p[0]}) {
            if ($o_null_fraction > 0.0 and !$INDELS{"$p[0]:$p[1]"} and rand() < $o_null_fraction) {
                my @result = do_unselected_test(\@p, undef, undef); 
                $result[7] = "null_$result[7]";
                print_result(\@result);
            }
            @p = read_unselected_line(); # advance $p because wrong ref sequence
        } else {
            croak "unknown condition involving reference order";
        }
    }
}


print STDERR "Skipped $N_indels indels\n";

print STDERR "Finished unselected test\n";

