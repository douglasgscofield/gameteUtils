#!/usr/bin/env perl

# Given heterozygous sites and profiles of Fin, Nonselected, Central, and
# Outer, produce output which produces Likelihoods for each heterozygous site
# to perform likelihood ratio tests.  Values marked '*' are coming and their
# places are marked with '.' in the output.

# 01. Ref
# 02. Pos
# 03. Genotype
# 04. Genotype likelihood
# 05. LLp, likelihood of data given model
# 06. LLm, likelihood of monomorphism
# 07. LRTm, LRT of monomorphism
# 08. 
#
#
# Heterozygous 0.5/0.5 null vs. Fin
# Fin vs. Nonselected
# Nonselected vs. Central
# Central vs. Outer
#
# Given the above files 01. ref 02. pos 03. genotype 04. genotype quality
#

# TODO: call out sites that are het calls but unavailable in one of the pools;
# currently these are simply skipped
# TODO: figure out why negative value passed to log0() for chr1:961063, skipped it for now
# TODO: produce solid estimate of e_hat, the best from all 239183 chr1 sites is 0.00343111
# TODO: filter based on coverage: mean Tcov = 246.938. 7264 sites with > 2*Tcov on chr1
# TODO: Why all the zero Maf estimates, eg, chr1:756{84,86,91,99} ?  Not a new feature

# Copyright (c) 2016 Douglas G. Scofield, Uppsala University
# douglas.scofield@ebc.uu.se, douglasgscofield@gmail.com
#
# Methods are inspired by Lynch et al. Genome Biol Evol 2014 doi:10.1093/gbe/evu085
#
my $starttime = localtime();
use strict;
use warnings;
use feature 'say';
use Carp;
use Getopt::Long;
use List::Util;
use FindBin;
use lib $FindBin::RealBin;  # add script directory to @INC to find BinomialTest and/or GameteUtils
use GameteUtils qw/ @base_a %base_h log0 choose fill_ref_index open_possibly_gzipped read_hetsites_line read_profile_line sorted_pool_counts rounddot multdot log0dot log10dot do_unselected_test /;

# 'our' variables below are also used within the GameteUtils package
my $o_fai_file;
my $o_hetsites_file;   # $h will be file handle for $o_hetsites_file
my $o_finclip_file;   # $f will be file handle for $o_finclip_file
my $o_unselected_file; # $u will be file handle for $o_unselected_file
my $o_central_file; # $c will be file handle for $o_central_file
my $o_outer_file; # $o will be file handle for $o_outer_file

my $o_ploidy = 1;
my $o_poolsize = 1000;
use constant MAX_POOLSIZE => 1000;
my $o_maxcov = 0;  # 2*mean total coverage (Tcov) of all sites across chr1
our $o_mincov = 32;  # total coverage (Tcov) of all sites across chr1
our $o_genotype = 1;  # check for incompatible genotypes?
our $o_max_allele_3 = 0.1;
my $o_null_fraction = 0.000;

our $o_sort_bases = 0;
our $o_hets_genotype = 0;  # GameteUtils::do_unselected_test() results do not include genotype columns
my $o_allsites = 1;
my $o_max_sites = 0;

my $o_help;
my $N_sites = 0;  # number of sites examined
my $N_maxcov = 0;   # number of sites excluded for being > $o_maxcov
my $N_mincov = 0;   # number of sites excluded for being < $o_mincov
my $current_reference = ""; # the name of the current reference sequence
my $N_references = 0; # the number of reference sequences seen
my $N_coordinates = 0; # the total number of coordinates seen
our $o_digits = 7; # number of digits to round output
our $N_indels = 0; # total number of indels skipped
our $o_indels = 1; # produce file of indels skipped
our $o_indels_annotate = 1; # annotate indels with quality, INFO and genotype columns from VCF
our $indels;  # $indels is file handle for indels file
our %INDELS;  # track indel-containing positions

my $usage_short = "
    $0 --fai reference.fa.fai --hetsites h.vcf --unselected u.pro --central central.pro --selected-2 outer.pro

        --fai FILE
        --hetsites FILE
        --unselected FILE
        --central FILE
        --outer FILE
        --ploidy INT             [default $o_ploidy]
        --pool-size INT          [default $o_poolsize, maximum ".MAX_POOLSIZE."]
        --maxcov FLOAT           [no default, use the following]
             male 31: 2*246.938 = 493.876
             male 32: 2*266.425 = 532.850
             male 34: 2*277.560 = 555.120
        --mincov INT             [default $o_mincov]
        --max-allele-3 FLOAT     [default $o_max_allele_3]
        --indels                 [default $o_indels]
        --no-indels
        --max-sites INT          [default $o_max_sites]
        --allsites               [default $o_allsites]
        --help, -?               More detailed help message

";

my $usage = "
NAME

  hetTwoPoolTest.pl - test list of heterozygous sites for deviations from 1:1 
                      in an 'unselected' pool and two 'selected' pools


SYNOPSIS

    $0 --fai ref.fa.fai --hetsites h.vcf --unselected u.pro --central central.pro --selected-2 outer.pro


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

    --finclip FILE     Profile base content for finclip reads used to call het sites

    --unselected FILE  Profile base content for 'unselected' pool
 
    --central FILE     Profile base content for central pool
 
    --outer FILE       Profile base content for outer pool

    --maxcov INT       Maximum total read coverage to consider a site, applied to the
                       sum of coverage across all four read pools [default $o_maxcov].

    --mincov INT       Minimum total read coverage to consider a site, applied to the
                       sum of coverage across all four read pools [default $o_mincov].

  OPTIONS FOR LIKELIHOOD ESTIMATION; default appropriate for sperm pools
    
    --ploidy INT       Ploidy of individuals within input pools [default $o_ploidy]

    --pool-size INT    Number of individuals within input pools [default $o_poolsize]

  OTHER OPTIONS

    --max-sites INT    Maximum number of sites (for debugging) [default $o_max_sites]

    --indels
    --no-indels   Write (or don't) apparent indels to an indels output file.  If the
                  output file cannot be created, this option is disabled [default $o_indels].

    --indels-annotate
    --no-indels-annotate
                  Annotate indels written above with quality, INFO and genotype
                  columns from the input VCF [default $o_indels_annotate]

    --help, -?    Help message

OUTPUT

A bunch of tab-separated columns:

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
    "finclip=s"           => \$o_finclip_file,
    "unselected=s"        => \$o_unselected_file,
    "central=s"           => \$o_central_file,
    "outer=s"             => \$o_outer_file,
    "ploidy=i"            => \$o_ploidy,
    "pool-size=i"         => \$o_poolsize,
    "maxcov=f"            => \$o_maxcov,
    "mincov=f"            => \$o_mincov,
    "indels"              => \$o_indels,
    "no-indels"           => sub { $o_indels = 0 },
    "indels-annotate"     => \$o_indels_annotate,
    "no-indels-annotate"  => sub { $o_indels_annotate = 0 },
    "allsites"            => \$o_allsites,
    "max-sites=i"         => \$o_max_sites,
    "help|?"              => \$o_help,
) or print_usage_and_exit(1);

print_usage_and_exit(1) if $o_help 
                          or not $o_fai_file
                          or not $o_hetsites_file
                          or not $o_unselected_file
                          or not $o_central_file
                          or not $o_outer_file
                          or $o_ploidy <= 0
                          or $o_poolsize <= 0
                          or $o_maxcov <= 0
                          or $o_max_sites < 0;
if (not defined($o_maxcov) or ! $o_maxcov) {
    print STDERR "No sensible --maxcov value provided.  Consider these, based on chr1:

    male 31: 2*246.938 = 493.876
    male 32: 2*266.425 = 532.850
    male 34: 2*277.560 = 555.120
";
    exit(1);
}

if ($o_poolsize > MAX_POOLSIZE) {
    print STDERR "*** Using large pool sizes has not yet been solved, please choose a value <= ",MAX_POOLSIZE,"\n";
    exit(1);
}

$o_indels ||= $o_indels_annotate;

# fill reference sequence order hash from fai file
my %REF_ORDER = fill_ref_index($o_fai_file);
print STDERR "Found ".scalar(keys(%REF_ORDER))." reference sequences in $o_fai_file\n";

my $h = open_possibly_gzipped($o_hetsites_file);

if ($o_indels) {
    $o_indels = $o_hetsites_file . "_hetPool_indels.txt";
    if (! open ($indels, ">", $o_indels)) {
        print STDERR "cannot open indels output file '$o_indels', disabling option\n"; 
        $o_indels = 0;
        $o_indels_annotate = 0;
    }
}

my $f = open_possibly_gzipped($o_finclip_file);
my $u = open_possibly_gzipped($o_unselected_file);
my $c = open_possibly_gzipped($o_central_file);
my $o = open_possibly_gzipped($o_outer_file);

sub read_finclip_line()    { return read_profile_line($f,  "finclip"); }

sub read_unselected_line() { return read_profile_line($u,  "unselected"); }

sub read_central_line() { return read_profile_line($c, "central"); }

sub read_outer_line() { return read_profile_line($o, "outer"); }

sub do_combined_test($$$$$) {
    my ($u, $l1, $l2, $h_ref, $h_alt) = @_;
    my @result_u = do_unselected_test($u, $h_ref, $h_alt); 
    # result_u has 6 fields, last is test state
    my @result_12 = do_twopool_test($l1, $l2, $h_ref, $h_alt); 
    # result_12 has 10 fields, last is test state
    croak "do_combined_test inconsistent ref/pos" if $result_u[0] ne $result_12[0] or $result_u[1] != $result_12[1];
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
start time                : $starttime
fai file                  : $o_fai_file
het sites                 : $o_hetsites_file
finclip profile           : $o_finclip_file
unselected profile        : $o_unselected_file
central profile           : $o_central_file
outer profile             : $o_outer_file
ploidy                    : $o_ploidy
pool-size                 : $o_poolsize
maxcov                    : $o_maxcov
mincov                    : $o_mincov
allsites                  : $o_allsites
indels                    : $o_indels
max_sites                 : $o_max_sites

";
my $outconfig = $config;
$outconfig =~ s/^/#/mg;
print STDERR $config;
print STDOUT $outconfig;

my @h = read_hetsites_line($h); # CHROM POS REF ALT
my @f = read_finclip_line();    # CHROM POS A C G T
my @u = read_unselected_line();
my @c = read_central_line();
my @o = read_outer_line();

sub print_result($) {  # do rounding, etc. for output
    my $r = shift;
    $r->[4] = rounddot($r->[4]);
    $r->[13] = rounddot($r->[13]);
    print STDOUT join("\t", @$r), "\n";
}

say join("\t", "ref", "pos", "Tcov", "R_A", "qual", "M_m", "Maf_F", "Maf_U", "Maf_C", "Maf_O", "e_hat", "lnL_Maf_F", "lnL_Maf_U", "lnL_Maf_C", "lnL_Maf_O", "LRT_FU", "LRT_UC", "LRT_CO");

while (@h and @f and @u and @c and @o) {
    if ($h[0] eq $f[0] and $h[0] eq $u[0] and $h[0] eq $c[0] and $h[0] eq $o[0]) { # same reference
        if ($h[1] == $f[1] and $h[1] == $u[1] and $h[1] == $c[1] and $h[1] == $o[1]) { # same position

            # Major and minor and error allele bases and profile indices, same for all at each site
            # Major and minor allele assigned according to counts within unselected pool
            my @Mm = sorted_pool_counts(\@u);
            my $Major = $base_a[ $Mm[3]->[1] ];  # most abundant base
            my $minor = $base_a[ $Mm[2]->[1] ];  # penultimate base
            my $i_Major = $base_h{$Major} + 2;   # indices into profile arrays
            my $i_minor = $base_h{$minor} + 2;
            my %error = %base_h; delete @error{($Major, $minor)};
            my @i_error = map { $base_h{$_} + 2 } keys %error;
            my $i_error = $#f + 1; # index to error read count added at end of profile arrays
            $f[$i_error] = $f[$i_error[0]] + $f[$i_error[1]];
            $u[$i_error] = $u[$i_error[0]] + $u[$i_error[1]];
            $c[$i_error] = $c[$i_error[0]] + $c[$i_error[1]];
            $o[$i_error] = $o[$i_error[0]] + $o[$i_error[1]];
            my $n_all_Major = $f[$i_Major] + $u[$i_Major] + $c[$i_Major] + $o[$i_Major];
            my $n_all_minor = $f[$i_minor] + $u[$i_minor] + $c[$i_minor] + $o[$i_minor];
            my $n_all_error = $f[$i_error] + $u[$i_error] + $c[$i_error] + $o[$i_error];
            my $n_all_Total = $n_all_Major + $n_all_minor + $n_all_error;

            if ($n_all_Total < $o_mincov) {  # ignore this site, coverage too low
                ++$N_mincov;
                goto READ_NEXT_SITE;
            }

            if ($o_maxcov and $n_all_Total > $o_maxcov) {  # ignore this site, coverage too high
                ++$N_maxcov;
                goto READ_NEXT_SITE;
            }

            sub L_Maf($$$$$) {
                my ($tag, $p_M_hat, $n_Major, $n_minor, $e_hat) = @_;
                return '.' if $p_M_hat eq '.';
                my $this_n = $o_ploidy * $o_poolsize;
                my $L = 0;
                for my $i (0..$this_n) {
                    $L += choose($this_n, $i) * 
                          ($p_M_hat ** $i) *
                          ((1.0 - $p_M_hat) ** ($this_n - $i)) *
                          (phi_Major(($i / $this_n), $e_hat) ** $n_Major) *
                          (phi_minor(($i / $this_n), $e_hat) ** $n_minor);
                    # say "$tag: after iter $i L=$L";
                }
                return $L;
            }
            sub e_hat($$) {
                # error rate estimate; Lynch et al. Eq 3a
                ### NOTE: $e_hat will be back-substituted with the $e_hat estimate averaged across all sites
                # first estimate from 17788 sites on chr1 is 0.00314945
                # estimate from all 239183 sites on chr1 is 0.00343111
                # estimate from all sites on chr1 without too much coverage is ????????
                # estimate from all sites on chr1 without too much coverage and Maf < 0.9 is ????????
                my ($n_error, $n_Total) = @_;
                return 0.0 if $n_Total == 0;
                return 3.0 * ($n_error / $n_Total) / 2.0; # error rate estimate; Lynch et al. Eq 3a
            }
            sub phi_Major($$) {
                # Lynch et al. Eq 1a
                my ($p, $e) = @_;
                return ($p * (1.0 - (4.0 * $e / 3.0))) + ($e / 3.0);
            }
            sub phi_minor($$) {
                # Lynch et al. Eq 1b
                my ($p, $e) = @_;
                return ($p * ((4.0 * $e / 3.0) - 1.0)) + (1 - $e);
            }
            sub p_hat($$$) {
                # estimated Major allele frequency in the pool; Lynch et al. Eq 3b
                my ($n_Major, $n_minor, $this_e_hat) = @_;
                # estimated of Major read fraction in error-free reads; Lynch et al. p.1211
                return '.' if $n_Major + $n_minor == 0;
                my $p_Major_hat = $n_Major / ($n_Major + $n_minor);
                # estimated Major allele frequency in the pool; Lynch et al. Eq 3b
                my $denom = (1.0 - (4.0 * $this_e_hat / 3.0));
                return '.' if $denom == 0.0;
                return (($p_Major_hat * (1.0 - (2.0 * $this_e_hat / 3.0))) - ($this_e_hat / 3.0)) / $denom;
            }
            sub LL_p($$$$$) {
                # Log-likelihood of data given read model; Lynch et al. Eq 4a
                my ($n_Major, $n_minor, $n_error, $p_hat, $this_e_hat) = @_;
                my $phi_Major_hat = phi_Major($p_hat, $this_e_hat);
                my $phi_minor_hat = phi_minor($p_hat, $this_e_hat);
                my $ans = 0;
                return ($n_Major * log0($phi_Major_hat)) +
                       ($n_minor * log0($phi_minor_hat)) + 
                       ($n_error * log0($this_e_hat / 3.0));
            }
            sub LL_m($$$) {
                # Log-likelihood of monomorphism of (site fixed for) Major allele; Lynch et al. Eq 4b
                my ($n_Major, $n_minor, $n_error) = @_;
                my $n_Total = $n_Major + $n_minor + $n_error;
                my $e_r_hat = ($n_Total - $n_Major) / $n_Total; # estimated proportion non-Major allele reads
                return ($n_Major * log0(1.0 - $e_r_hat)) +
                       (($n_Total - $n_Major) * log0($e_r_hat / 3.0));
            }
            sub LRT_Maf($$$$$$$) {
                my ($a_L_Maf, $a_Major, $a_minor, $b_L_Maf, $b_Major, $b_minor, $this_e_hat) = @_;
                return '.' if $a_L_Maf eq '.' or $b_L_Maf eq '.';
                # summed profile quartet over both pools
                my $n_T_Major = $a_Major + $b_Major;
                my $n_T_minor = $a_minor + $b_minor;
                # estimate $p_T_hat for the summed quartet (for assuming homogeneity)
                my $p_T_hat = p_hat($n_T_Major, $n_T_minor, $this_e_hat);
                # Likelihoods of both given homogeneity
                my $L_01 = L_Maf("L_01", $p_T_hat, $a_Major, $a_minor, $this_e_hat);
                my $L_02 = L_Maf("L_02", $p_T_hat, $b_Major, $b_minor, $this_e_hat);
                # Likelihoods of both given heterogeneity are in $a_L_Maf and $b_L_Maf
                return 2.0 * (log0($a_L_Maf * $b_L_Maf) - log0($L_01 * $L_02));
            }

            # e_hat in output is estimated from all pooled base counts
            my $e_hat = e_hat($n_all_error, $n_all_Total);

            my $f_p_hat = p_hat($f[$i_Major], $f[$i_minor], $e_hat);
            my $u_p_hat = p_hat($u[$i_Major], $u[$i_minor], $e_hat);
            my $c_p_hat = p_hat($c[$i_Major], $c[$i_minor], $e_hat);
            my $o_p_hat = p_hat($o[$i_Major], $o[$i_minor], $e_hat);

            my $f_L_Maf = L_Maf("f", $f_p_hat, $f[$i_Major], $f[$i_minor], $e_hat);
            my $u_L_Maf = L_Maf("u", $u_p_hat, $u[$i_Major], $u[$i_minor], $e_hat);
            my $c_L_Maf = L_Maf("c", $c_p_hat, $c[$i_Major], $c[$i_minor], $e_hat);
            my $o_L_Maf = L_Maf("o", $o_p_hat, $o[$i_Major], $o[$i_minor], $e_hat);

            my $LRT_FU = LRT_Maf($f_L_Maf, $f[$i_Major], $f[$i_minor], $u_L_Maf, $u[$i_Major], $u[$i_minor], $e_hat);
            my $LRT_UC = LRT_Maf($u_L_Maf, $u[$i_Major], $u[$i_minor], $c_L_Maf, $c[$i_Major], $c[$i_minor], $e_hat);
            my $LRT_CO = LRT_Maf($c_L_Maf, $c[$i_Major], $c[$i_minor], $o_L_Maf, $o[$i_Major], $o[$i_minor], $e_hat);

            #say join("\t", "ref", "pos", "Tcov", "R_A", "qual", "M_m", "Maf_F", "Maf_U", "Maf_C", "Maf_O", "e_hat", "lnL_Maf_F", "lnL_Maf_U", "lnL_Maf_C", "lnL_Maf_O", "LRT_FU", "LRT_UC", "LRT_CO");

            say join("\t", 
                $h[0],  # ref
                $h[1],   # pos
                $n_all_Total,  # total coverage
                "$h[2]_$h[3]", # REF_ALT alleles
                $h[4],  # het genotype quality
                "${Major}_${minor}", # Major_minor alleles
                rounddot($f_p_hat, 4),
                rounddot($u_p_hat, 4),
                rounddot($c_p_hat, 4),
                rounddot($o_p_hat, 4),
                rounddot($e_hat,   5),  # e_hat estimate from all reads
                rounddot(log0dot($f_L_Maf), 8),
                rounddot(log0dot($u_L_Maf), 8),
                rounddot(log0dot($c_L_Maf), 8),
                rounddot(log0dot($o_L_Maf), 8),
                rounddot($LRT_FU,  5),  # LRT U vs F
                rounddot($LRT_UC,  5),  # LRT C vs U
                rounddot($LRT_CO,  5),  # LRT O vs C
            );

            print STDERR "completed $N_sites @ $h[0] : $h[1]\n" if (++$N_sites % 1000) == 0;

            if ($o_max_sites and $N_sites >= $o_max_sites) {
                print STDERR "reached --max-sites $o_max_sites\n";
                goto TERMINATE_SCRIPT;
            }

            #
            READ_NEXT_SITE:    # to make quitting early easier
            #

            @h = read_hetsites_line($h);
            @f = read_finclip_line();
            @u = read_unselected_line();
            @c = read_central_line();
            @o = read_outer_line();

        } elsif ($h[1] < $f[1]) {

            @h = read_hetsites_line($h); # advance hetsites file because wrong pos

        } elsif ($f[1] != $u[1] or $f[1] != $c[1] or $f[1] != $o[1]) {

            # sync up pools
            do {
                if    ($f[1] < $u[1]) { @f = read_finclip_line(); }
                elsif ($u[1] < $c[1]) { @u = read_unselected_line(); }
                elsif ($c[1] < $o[1]) { @c = read_central_line(); }
                elsif ($o[1] < $f[1] or
                       $o[1] < $u[1] or
                       $o[1] < $c[1]) { @o = read_outer_line(); }
            } until ($f[1] == $u[1] and $f[1] == $c[1] and $f[1] == $o[1]);

        } elsif ($h[1] > $f[1]) {

            # hetsites is ahead of the pools, and we know the pools are synced

            @f = read_finclip_line();    # advance pools because wrong pos
            @u = read_unselected_line();
            @c = read_central_line();
            @o = read_outer_line();

        } else {
            croak "unknown condition involving positions";
        }

    } elsif ($REF_ORDER{$h[0]} < $REF_ORDER{$f[0]}) {

        @h = read_hetsites_line($h); # advance $h because wrong ref

    } elsif ($REF_ORDER{$f[0]} != $REF_ORDER{$u[0]} or
             $REF_ORDER{$f[0]} != $REF_ORDER{$c[0]} or
             $REF_ORDER{$f[0]} != $REF_ORDER{$o[0]}) {

        # sync up pools
        do {
            if    ($REF_ORDER{$f[0]} < $REF_ORDER{$u[0]}) { @f = read_finclip_line(); }
            elsif ($REF_ORDER{$u[0]} < $REF_ORDER{$c[0]}) { @u = read_unselected_line(); }
            elsif ($REF_ORDER{$c[0]} < $REF_ORDER{$o[0]}) { @c = read_central_line(); }
            elsif ($REF_ORDER{$o[0]} < $REF_ORDER{$f[0]} or
                   $REF_ORDER{$o[0]} < $REF_ORDER{$u[0]} or
                   $REF_ORDER{$o[0]} < $REF_ORDER{$c[0]}) { @o = read_outer_line(); }
        } until ($REF_ORDER{$f[0]} == $REF_ORDER{$u[0]} and
                 $REF_ORDER{$f[0]} == $REF_ORDER{$c[0]} and
                 $REF_ORDER{$f[0]} == $REF_ORDER{$o[0]});

    } elsif ($REF_ORDER{$h[0]} > $REF_ORDER{$f[0]}) {

        # hetsites is ahead of the pools, and we know the pools are synced

        @f = read_finclip_line(); # advance pools because wrong ref
        @u = read_unselected_line();
        @c = read_central_line();
        @o = read_outer_line();

    } else {
        croak "unknown condition involving reference order";
    }
}

TERMINATE_SCRIPT:

print STDERR "Skipped $N_indels indels\n";
print STDERR "Processed $N_sites total sites\n";
print STDERR "Skipped $N_maxcov sites for coverage > $o_maxcov\n";
print STDERR "Skipped $N_mincov sites for coverage < $o_mincov\n";

print STDERR "End time : ".localtime()."\n";
print STDERR "Finished hetPoolLikelihoods test\n";

