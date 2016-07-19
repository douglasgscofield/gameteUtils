#!/usr/bin/env perl

package GameteUtils;

use strict;
use warnings;
use Carp;
use POSIX;  # no export, use POSIX::log10()
use FindBin;
use lib $FindBin::RealBin;  # add script directory to @INC to find BinomialTest
use BinomialTest;  # no export, use BinomialTest::binomial_test()


# package packaging

use Exporter qw/ import /;
our @EXPORT_OK = qw/ fill_ref_index open_possibly_gzipped read_hetsites_line read_profile_line sorted_pool_counts rounddot multdot log10dot do_unselected_test /;


# package contents


our @base = qw/ A C G T /;
our %base = ( A => 0, C => 1, G => 2, T => 3 );


# sub fill_ref_index($)
#
# Argument 1: Fasta index (.fai) filename to open
# Return    : hash filled with reference sequence order (list context), or
#             reference to this hash (scalar context)
#
# Main script 'our' variables: NONE
sub fill_ref_index($) {
    my $fai_file = shift;
    my %REF_ORDER;
    my $ref_index = 0;
    open (my $ref, "<", $fai_file) or croak "cannot open Fasta index file '$fai_file': $!";
    while (<$ref>) {
        my @l = split /\t/;
        $REF_ORDER{$l[0]} = ++$ref_index if not exists $REF_ORDER{$l[0]};
    }
    return wantarray ? %REF_ORDER : \%REF_ORDER;
}


# sub open_possibly_gzipped($)
#
# Argument 1: filename to open, must end with .gz if gzipped
# Return    : file handle to opened file
#
# Main script 'our' variables: NONE
sub open_possibly_gzipped($) {
    my ($file, $fh) = shift;
    if ($file =~ /\.gz$/) {
        open ($fh, "-|", "gzip -dc $file") or croak "cannot open '$file': $!";
    } else {
        open ($fh, "<", $file) or croak "cannot open '$file': $!";
    }
    return $fh;
}


# sub read_hetsites_line($)
#
# Argument 1: file handle to open hetsites VCF file
# Return    : list containing het site info: CHROM POS REF ALT, with REF and ALT
#             switched if $main::o_sort_bases and REF is lexically smaller than ALT
#
# Writes indels to $main::indels file if appropriate options set in main script.
#
# Main script 'our' variables:
#     $main::o_indels
#     $main::o_indels_annotate
#     $main::indels
#     $main::N_indels
#     $main::INDELS
#     $main::o_sort_bases
sub read_hetsites_line($) {
    my $fh = shift;
    READ_LINE:
    my $l = <$fh>;
    while ($l and $l =~ /^#/) {
        $l = <$fh>;
    }
    return () if ! $l;
    chomp $l;
    #print STDERR "GameteUtils::read_hetsites_line: '$l'\n";
    my @l = split /\t/, $l;
    if (length($l[3]) > 1 or length($l[4]) > 1) {  # ref or alt not a single base
        ++$main::N_indels;
        ++$main::INDELS{"$l[0]:$l[1]"};  # mark this position as an indel
        if ($main::o_indels) {
            my @i = ("INDEL", @l[(0,1,3,4)]);
            push @i, @l[5,7,8,9] if $main::o_indels_annotate;
            print { $main::indels } join("\t", @i), "\n";
        }
        goto READ_LINE;
    }
    # return CHROM POS REF ALT
    ($l[3], $l[4]) = ($l[4], $l[3]) if $main::o_sort_bases and $base{$l[3]} > $base{$l[4]};
    return @l[(0,1,3,4)];
}


# sub read_profile_line($$)
#
# Argument 1: file handle to open profile file
# Argument 2: tag string to apply to output messages
# Return    : list containing profile line contents, CHROM POS A C G T [N]
#
# Main script 'our' variables: NONE
sub read_profile_line($$) {
    my ($fh, $tag) = @_;
    my $l = <$fh>;
    return () if ! $l;
    chomp $l;
    #print STDERR "GameteUtils::read_profile_line:$tag: '$l'\n";
    return split /\t/, $l;
}


# sub sorted_pool_counts($)
#
# Argument 1: reference to list holding contents of profile line
# Return    : list of anonymous arrays of read counts and base id,
#             sorted by base count:
#             ( [ lowest-count , base-id ],
#               [ low-count    , base-id ],
#               [ high-count   , base-id ],
#               [ highest-count, base-id ] )
#             The order of ONLY the high-count and highest-count bases
#             is swapped if $main::o_sort_bases and the bases are lexically
#             out of order
#
# Main script 'our' variables:
#     $main::o_sorted_bases
sub sorted_pool_counts($) {
    my $l1 = shift;
    my @d = ( [ $l1->[2], 0 ],   # A
              [ $l1->[3], 1 ],   # C
              [ $l1->[4], 2 ],   # G
              [ $l1->[5], 3 ] ); # T
    @d = sort { $a->[0] <=> $b->[0] } @d;
    ($d[2], $d[3]) = ($d[3], $d[2]) if $main::o_sort_bases and $d[3]->[1] > $d[2]->[1];
    return @d;
}


# sub sorted_twopool_counts($)
#
# Argument 1: reference to list holding contents of profile line
# Argument 2: reference to list holding contents of profile line
# Return    : list of anonymous arrays of read counts and base id,
#             sorted by base count:
#             ( [ count-1, count-2, lowest-count-sum , base-id ],
#               [ count-1, count-2, low-count-sum    , base-id ],
#               [ count-1, count-2, high-count-sum   , base-id ],
#               [ count-1, count-2, highest-count-sum, base-id ] )
#             The order of ONLY the high-count-sum and highest-count-sum bases
#             is swapped if $main::o_sort_bases and the bases are lexically
#             out of order
#
# Main script 'our' variables:
#     $main::o_sorted_bases
sub sorted_twopool_counts($$) {
    my ($l1, $l2) = @_;
    croak "inconsistent ref/pos" if $l1->[0] ne $l2->[0] or $l1->[1] != $l2->[1];
    my @d = ( [ $l1->[2], $l2->[2], $l1->[2] + $l2->[2], 0 ],   # A
              [ $l1->[3], $l2->[3], $l1->[3] + $l2->[3], 1 ],   # C
              [ $l1->[4], $l2->[4], $l1->[4] + $l2->[4], 2 ],   # G
              [ $l1->[5], $l2->[5], $l1->[5] + $l2->[5], 3 ] ); # T
    @d = sort { $a->[2] <=> $b->[2] } @d;
    ($d[2], $d[3]) = ($d[3], $d[2]) if $main::o_sort_bases and $d[2]->[3] < $d[3]->[3];
    return @d;
}


# sub rounddot($)
#
# Argument 1: floating point value, or "."
# Return    : argument rounded to $main::o_digits, or "."
#             if argument was "."
#
# Main script 'our' variables:
#     $main::o_digits
sub rounddot($) {
    my $f = shift;
    return $f eq "." ? $f : sprintf("%.".$main::o_digits."f", $f);
}


# sub multdot($$)
#
# Argument 1: floating point value, or "."
# Argument 2: floating point value, or "."
# Return    : arguments 1 and 2 multiplied, or "." if either was "."
#
# Main script 'our' variables: NONE
sub multdot($$) {
    my ($p1, $p2) = @_;
    my $ans = ($p1 eq "." or $p2 eq ".") ? "." : ($p1 * $p2);
    return $ans;
}


# sub log10dot($)
#
# Argument 1: floating point value, or "."
# Return    : argument with log10() applied, or "." if argument was "."
#
# Main script 'our' variables: NONE
sub log10dot($) {
    my $f = shift;
    return $f eq "." ? $f : POSIX::log10($f);
}


# do_unselected_test
#
# Argument 1: reference to list holding contents of profile line
# Argument 2: reference base from het VCF [optional, may be undef]
# Argument 3: alternate base from het VCF [optional, may be undef]
# Return    : 6-element list containing results of the test
#             0: reference
#             1: position
#             2: read coverage
#             3: test description string
#             4: unselected probability
#             5: test result string
#
# Main script 'our' variables:
#     $main::o_mincov
#     $main::o_genotype
#     $main::o_max_allele_3
#     $main::o_unselected_prob
sub do_unselected_test {
    my ($l, $h_ref, $h_alt) = @_;
    my @d = sorted_pool_counts($l);
    my $cov = $d[2]->[0] + $d[3]->[0];

    my $allele1 = $base[$d[3]->[1]];
    my $allele2 = $base[$d[2]->[1]];
    my $otot = "unsel,$allele1/$allele2:$d[3]->[0]/$d[2]->[0]";

    return ($l->[0], $l->[1], $cov, $otot, ".", "zerocov")
        if ($cov == 0); # zero coverage

    my $prob = BinomialTest::binomial_test($d[2]->[0], $cov);

    return ($l->[0], $l->[1], $cov, $otot, $prob, "mincov")
        if ($cov < $main::o_mincov); # insufficient coverage

    if ($main::o_genotype and defined($h_ref) and defined($h_alt)) {
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
        if ($d[1]->[0] >= $cov * $main::o_max_allele_3); # allele 3 coverage too high

    return ($l->[0], $l->[1], $cov, $otot, $prob, "binom_fail")
        if ($prob < $main::o_unselected_prob); # effectively homozygous site

    return ($l->[0], $l->[1], $cov, $otot, $prob, "binom_pass");
}


# sub do_twopool_test($$$$)
#
# Argument 1: reference to list holding contents of pool 1 profile line
# Argument 2: reference to list holding contents of pool 2 profile line
# Argument 3: reference base from het VCF
# Argument 4: alternate base from het VCF
# Return    : 10-element list containing results of the test
#             0: reference
#             1: position
#             2: read coverage of allele1 and allele2 in selected 1 pool
#             3: read coverage of allele1 and allele2 in selected 2 pool
#             4: test description string; combination of strings for each pool
#             5: selected 1 and 2 pool unselected test results
#             6: selected 1 and 2 pool binomial probabilities
#             7: selected 1 and 2 pool result of <=> comparisons of allele counts
#             8: twopool_test result string
#             9: multiplied probabilities
#
# Main script 'our' variables: several others used by called subs
#     $main::o_selected_prob
sub do_twopool_test($$$$) {
    my ($l1, $l2, $h_ref, $h_alt) = @_;

    my @d      = sorted_twopool_counts($l1, $l2);
    my @utest1 = do_unselected_test($l1, $h_ref, $h_alt);
    my @utest2 = do_unselected_test($l2, $h_ref, $h_alt);

    my $cov1 = $utest1[2];
    my $cov2 = $utest2[2];

    # test; unselected test pool 1 description; unselected test pool 2 description
    my $desc = "twoPool;$utest1[3];$utest2[3]";

    # unselected test pool 1 test result; unselected test pool 2 test result
    my $result12 = "$utest1[5];$utest2[5]";

    # unselected test pool 1 probability; unselected test pool 2 probability;
    my $prob12 = rounddot($utest1[4]) . ";" . rounddot($utest2[4]);

    # two selected pools in allele count opposition
    my $counts1 = $d[3]->[0] <=> $d[2]->[0];  # -1 if 3 < 2, 0 if 3 == 3, 1 if 3 > 2
    my $counts2 = $d[3]->[1] <=> $d[2]->[1];
    my $counts = sprintf("%+d;%+d", $counts1, $counts2);

    # what is going to be the result of this test? check returns of unselected tests
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
    } elsif ($utest1[4] > $main::o_selected_prob or $utest2[4] > $main::o_selected_prob) {
        $result = $utest1[4] > $main::o_selected_prob ? "failprob1" : "passprob1";
        $result .= ";";
        $result .= $utest2[4] > $main::o_selected_prob ? "failprob2" : "passprob2";
    } else {
        $result = "twopool_pass";
    }

    # two selected pools likelihood
    my $prob = multdot($utest1[4], $utest2[4]);

    return ($l1->[0], $l1->[1], $cov1, $cov2, $desc, $result12, $prob12, $counts, $result, $prob);
}


1;
