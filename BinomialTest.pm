#!/usr/bin/env perl

package BinomialTest;

use strict;
use warnings;
use Carp;
use Math::CDF qw/ pbinom /;
use Exporter qw/ import /;
our @EXPORT_OK = qw/ binomial_test binomial_test2 /;

our %BT;
our %BT2;

sub binomial_test($$) {
    my ($x, $n) = @_;
    # x  number of successes
    # n  number of trials
    if (exists $BT{$n}{$x}) {
        #print STDERR "found [ $x, $n ]\n";
        return $BT{$n}{$x};
    }
    if ($x > $n or $n <= 0 or $x < 0) {
        carp "invalid number of successes or trials (x = $x, n = $n)";
        return -1;
    }
    $x = $n - $x if $x > $n / 2.0;  # reverse x
    my $prob = pbinom($x, $n, 0.5);
    # $prob = 1 if $prob > 1;  # unnecessary with 1-tailed test
    $BT{$n}{$x} = $BT{$n}{$n - $x} = $prob;
    return $prob;
}

sub binomial_test2($$) {
    my ($x, $n) = @_;
    # x  number of successes
    # n  number of trials
    if (exists $BT2{$n}{$x}) {
        #print STDERR "found [ $x, $n ]\n";
        return $BT2{$n}{$x};
    }
    if ($x > $n or $n <= 0 or $x < 0) {
        carp "invalid number of successes or trials (x = $x, n = $n)";
        return -1;
    }
    $x = $n - $x if $x > $n / 2.0;  # reverse x
    my $prob = pbinom($x, $n, 0.5) * 2;  # double prob, the "2" bit
    $prob = 1 if $prob > 1;
    $BT2{$n}{$x} = $BT2{$n}{$n - $x} = $prob;
    return $prob;
}

sub test {
    my @t = (
        [3, 6],
        [3, 8], 
        [3, 10], 
        [3, 20], 
        [13, 20], 
        [7, 20], 
    );

    foreach my $m (@t) {
        my ($x, $n) = @$m;
        print "binom.test($x, $n, 0.5) = ", binomial_test($x, $n), "\n";
        print "binom.test2($x, $n, 0.5) = ", binomial_test2($x, $n), "\n";
    }

}

1;
