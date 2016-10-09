#!/usr/bin/env perl

use strict;
use warnings;
use feature 'say';
use Carp;
use Getopt::Long;
use List::Util;
use POSIX; # for POSIX::floor()
use FindBin;
use lib $FindBin::RealBin;  # add script directory to @INC
use GameteUtils qw/ fill_ref_lengths open_possibly_gzipped rounddot /;


my $o_windowsize = 50000;
my $o_windowoverlap = 0;
my $o_windowminsize = 10000;
my $o_max_windows = 0;
my $o_mincov = 40;
my $N_mincov = 0;
my $o_minsites = 5;
my $N_minsites = 0;
my $N_zerosites = 0;
my $o_fai_file;
my $o_input_file;

my $usage = "$0 [options] --fai fasta.fa.fai --input FILE

OPTIONS:

    --fai FILE        Fasta index file (samtools faidx FASTA.fa) [default $o_fai_file]
    --input FILE      hetPoolLikelihoods.pl output file [default $o_input_file]
    --size INT        Window size in bp [default $o_windowsize]
    --overlap INT     Window overlap in bp [default $o_windowoverlap]
    --minsize INT     Minimum window size [default $o_windowminsize]
    --max-windows INT Maximum number of windows to output [debugging, default $o_max_windows]
    --mincov INT      Minimum total coverage required to keep a site [default $o_mincov]
    --minsites INT    Minimum number of valid sites required to keep a window [default $o_minsites]
";

GetOptions(
    "fai=s"         => \$o_fai_file,
    "input=s"       => \$o_input_file,
    "size=i"        => \$o_windowsize,
    "overlap=i"     => \$o_windowoverlap,
    "minsize=i"     => \$o_windowminsize,
    "max-windows=i" => \$o_max_windows,
    "mincov=i"      => \$o_mincov,
    "minsites=i"    => \$o_minsites,
) or do { say STDERR $usage; exit(1); };

# fill reference sequence order hash from fai file
my %REF_LENGTHS = fill_ref_lengths($o_fai_file);
say STDERR "Found ".scalar(keys(%REF_LENGTHS))." reference sequences in $o_fai_file";

croak "Please use sensible window and overlap sizes" if $o_windowoverlap < 0 or $o_windowoverlap >= $o_windowsize or $o_windowminsize <= 0;

my $input_fh = open_possibly_gzipped($o_input_file);

my $config = "Starting windowify
fai file                  : $o_fai_file
input                     : $o_input_file
window size               : $o_windowsize
window overlap            : $o_windowoverlap
window min size           : $o_windowminsize
max-windows               : $o_max_windows
minimum total coverage    : $o_mincov
minimum sites per window  : $o_minsites

";
my $outconfig = $config;
$outconfig =~ s/^/#/mg;
print STDERR $config;
print STDOUT $outconfig;

my $l = <$input_fh>;
while ($l and $l =~ /^#/) { $l = <$input_fh>; }  # skip initial comments
chomp $l; # we stopped on the header
my @header = split /\t/, $l;  # array mapping column to header name
my %header = map { $header[$_] => $_ } 0..$#header;  # hash mapping name to header column

my $N_windows = 0;

sub in_window($$) {
    # is this ref and pos within this window?
    my ($l, $w) = @_;
    return 1 if %$w and $l->[0] eq $w->{ref} and $l->[1] >= $w->{start} and $l->[1] <= $w->{end};
    return 0;
}
sub next_window_on_ref($$) {
    # arg1: reference name   arg2: hash of current window
    # return a hash containing the next window relative to this one, empty if none on this ref
    # if arg2 undef, then initiate window on reference
    my ($ref, $w) = @_;
    my %n;
    croak "ref empty" if !$ref;
    croak "this-window hash ref empty" if defined($w) and !%$w;
    if (not defined($w)) { # start first window
        #say STDERR "next_window_on_ref: initiating first window on ref $ref";
        $n{ref} = $ref;
        $n{start} = 1;
    } else {  # continue from window in %$w
        if ($ref ne $w->{ref}) {
            # the next site (the source of $ref) is beyond the reference sequence of %$w
            # just return this window and let the control logic advance to the next reference
            # croak "\$ref ($ref) does not match \$w->{ref} ($w->{ref})";
            say STDERR "site-based \$ref ($ref) does not match \$w->{ref} ($w->{ref}), returning next window anyway";
        }
        #say STDERR "next_window_on_ref: next window computed from \%w: ".(%$w ? "ref=$w->{ref} start=$w->{start} end=$w->{end} size=$w->{size}" : "empty");
        $n{ref} = $w->{ref};
        $n{start} = $w->{end} + 1 - $o_windowoverlap;
    }
    $n{end} = List::Util::min($n{start} + $o_windowsize - 1, $REF_LENGTHS{$n{ref}});
    $n{size} = $n{end} - $n{start} + 1;
    if ($n{start} > $REF_LENGTHS{$n{ref}} or $n{size} < $o_windowminsize) {
        %n = ();
    }
    #say STDERR "next_window_on_ref: *** new window \%n: ".(%n ? "ref=$n{ref} start=$n{start} end=$n{end} size=$n{size}" : "empty");
    return wantarray ? %n : \%n;
}

sub invalid_site($) {
    my $l = shift;
    return 1 if $l->[$header{LRT_FU}] eq '.' or $l->[$header{LRT_FU}] eq 'nan' or
                $l->[$header{LRT_UC}] eq '.' or $l->[$header{LRT_UC}] eq 'nan' or
                $l->[$header{LRT_CO}] eq '.' or $l->[$header{LRT_CO}] eq 'nan' or
                $l->[$header{Maf_F}]  eq '.' or $l->[$header{Maf_F}]  < 0 or
                $l->[$header{Maf_U}]  eq '.' or $l->[$header{Maf_U}]  < 0 or
                $l->[$header{Maf_C}]  eq '.' or $l->[$header{Maf_C}]  < 0 or
                $l->[$header{Maf_O}]  eq '.' or $l->[$header{Maf_O}]  < 0 or
                $l->[$header{lnL_Maf_F}] eq '.' or
                $l->[$header{lnL_Maf_U}] eq '.' or
                $l->[$header{lnL_Maf_C}] eq '.' or
                $l->[$header{lnL_Maf_O}] eq '.';
    return 0;
}

my @cols = qw/ ref start end size Nval Ninv Nv_1k Ni_1k Tcov SDcov Maf_F Maf_U Maf_C Maf_O LRT_FU LRT_UC LRT_CO LRT_Maf_FU LRT_Maf_UC LRT_Maf_CO /;

sub finish_window($) {
    my $w = shift;
    if (not defined($w->{Nval}) or $w->{Nval} == 0) {
        ++$N_zerosites;
        return,
    }
    if ($o_minsites and $w->{Nval} < $o_minsites) {
        ++$N_minsites;
        return,
    }
    $w->{Ninv} = 0 if not exists($w->{Ninv});
    $w->{Nt_1k} = $w->{Ntot} / ($w->{size} / 1000.0);
    $w->{Nv_1k} = $w->{Nval} / ($w->{size} / 1000.0);
    $w->{Ni_1k} = $w->{Ninv} / ($w->{size} / 1000.0);
    $w->{Tcov} /= $w->{Nval};
    $w->{SDcov} = sqrt($w->{VM2} / ($w->{Nval} - 1));
    $w->{Maf_F} /= $w->{Nval};
    $w->{Maf_U} /= $w->{Nval};
    $w->{Maf_C} /= $w->{Nval};
    $w->{Maf_O} /= $w->{Nval};
    $w->{LRT_FU} /= $w->{Nval};
    $w->{LRT_UC} /= $w->{Nval};
    $w->{LRT_CO} /= $w->{Nval};
    $w->{LRT_Maf_FU} = 2*($w->{lnL_Maf_U} - $w->{lnL_Maf_F});
    $w->{LRT_Maf_UC} = 2*($w->{lnL_Maf_C} - $w->{lnL_Maf_U});
    $w->{LRT_Maf_CO} = 2*($w->{lnL_Maf_O} - $w->{lnL_Maf_C});

    $w->{Nt_1k} = rounddot($w->{Nt_1k}, 3);
    $w->{Nv_1k} = rounddot($w->{Nv_1k}, 3);
    $w->{Ni_1k} = rounddot($w->{Ni_1k}, 3);
    $w->{Tcov}  = rounddot($w->{Tcov},  2);
    $w->{SDcov} = rounddot($w->{SDcov}, 2);
    $w->{Maf_F} = rounddot($w->{Maf_F}, 4);
    $w->{Maf_U} = rounddot($w->{Maf_U}, 4);
    $w->{Maf_C} = rounddot($w->{Maf_C}, 4);
    $w->{Maf_O} = rounddot($w->{Maf_O}, 4);
    $w->{LRT_FU} = rounddot($w->{LRT_FU}, 5);
    $w->{LRT_UC} = rounddot($w->{LRT_UC}, 5);
    $w->{LRT_CO} = rounddot($w->{LRT_CO}, 5);
    $w->{LRT_Maf_FU} = rounddot($w->{LRT_Maf_FU}, 5);
    $w->{LRT_Maf_UC} = rounddot($w->{LRT_Maf_UC}, 5);
    $w->{LRT_Maf_CO} = rounddot($w->{LRT_Maf_CO}, 5);
    ++$N_windows;
    return join("\t", @{$w}{@cols});
}

sub add_to_window($$) {
    my ($l, $w) = @_;
    $w->{Ntot}++;
    if (invalid_site($l)) {
        #say STDERR "$l->[0] : $l->[1] is invalid";
        $w->{Ninv}++;
        return;
    }
    $w->{Nval}++;  # valid heterozygous sites
    $w->{Tcov}  += $l->[$header{Tcov}];
    # variance in coverage using Welford's algorithm https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    $w->{Vmean} = 0 if not exists($w->{Vmean});
    my $delta = $l->[$header{Tcov}] - $w->{Vmean};
    $w->{Vmean} += $delta / $w->{Nval};
    $w->{VM2}   += $delta * ($l->[$header{Tcov}] - $w->{Vmean});;
    # frequency and likelihood estimators
    $w->{Maf_F} += $l->[$header{Maf_F}];
    $w->{Maf_U} += $l->[$header{Maf_U}];
    $w->{Maf_C} += $l->[$header{Maf_C}];
    $w->{Maf_O} += $l->[$header{Maf_O}];
    $w->{e_hat} += $l->[$header{e_hat}];
    $w->{lnL_Maf_F} += $l->[$header{lnL_Maf_F}];
    $w->{lnL_Maf_U} += $l->[$header{lnL_Maf_U}];
    $w->{lnL_Maf_C} += $l->[$header{lnL_Maf_C}];
    $w->{lnL_Maf_O} += $l->[$header{lnL_Maf_O}];
    $w->{LRT_FU} += $l->[$header{LRT_FU}];
    $w->{LRT_UC} += $l->[$header{LRT_UC}];
    $w->{LRT_CO} += $l->[$header{LRT_CO}];
}
my %wcurr;
my %wnext;
my $oldref;
my $thisref;
my $out;

say join("\t", @cols);
while ($l = <$input_fh>) {
    next if $l =~ /^#/; # skip interstitial comments
    chomp $l;
    my @l = split /\t/, $l;
    NEWLINE:
    if ($l[$header{Tcov}] < $o_mincov) {
        ++$N_mincov;
        next;
    }
    if (! %wcurr) {
        for (%wcurr = next_window_on_ref($l[0], undef);
            %wcurr and ! in_window(\@l, \%wcurr);
            %wcurr = next_window_on_ref($l[0], \%wcurr)) {
            #say STDERR "*** advancing over empty window with empty wcurr";
        }
        if (! %wcurr) {
            # we can't get a suitable window for this ref
            # it could be that we have a ref with a site but the ref is not long enough for a window
            # if this is the case, we need to advance to the next ref
            $thisref = $l[0];
            while ($l[0] eq $thisref) {
                $l = <$input_fh>;
                $l = <$input_fh> until $l !~ /^#/; # skip interstitial comments
                chomp $l;
                @l = split /\t/, $l;
            }  # at end of loop, this site is at the next reference but windows are not
            last if ! $l;
            goto NEWLINE;
        }
        #last if ! %wcurr;
        %wnext = next_window_on_ref($l[0], \%wcurr);
    }
    if (in_window(\@l, \%wnext)) {
        add_to_window(\@l, \%wnext);
    }
    if (in_window(\@l, \%wcurr)) {
        add_to_window(\@l, \%wcurr);
        $thisref = $l[0];
    } else {
        # we must advance so that %wcurr is the window in which this site belongs
        # first, finish the current one
        $out = finish_window(\%wcurr);
        say $out if $out;
        # switch to new %wcurr
        for (%wcurr = %wnext;
            %wcurr and ! in_window(\@l, \%wcurr);
            %wcurr = next_window_on_ref($l[0], \%wcurr)) {
            #say STDERR "*** advancing over empty window to find site's wcurr";
        }
        if (! %wcurr) {
            # result of ! %wnext, because we could not find another window on this ref
            # or ! in_window(), because this site is in the next ref and not in the current window
            # advance sites so we are off this ref.  we lose the current site.
            while ($l[0] eq $thisref) {
                $l = <$input_fh>;
                $l = <$input_fh> until $l !~ /^#/; # skip interstitial comments
                chomp $l;
                @l = split /\t/, $l;
            }  # at end of loop, this site is at the next reference but windows are not
            goto NEWLINE;
        }
        last if ! %wcurr;
        %wnext = next_window_on_ref($wcurr{ref}, \%wcurr) if %wcurr;
    }
    last if $o_max_windows and $N_windows >= $o_max_windows;
}
$out = finish_window(\%wcurr);
say $out if $out;

say STDERR "Skipped $N_mincov sites with coverage < $o_mincov";
say STDERR "Output $N_windows windows";
say STDERR "Skipped $N_zerosites windows with zero sites";
say STDERR "Skipped $N_minsites windows with number of valid sites < $o_minsites";
