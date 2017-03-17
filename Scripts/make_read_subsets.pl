#!/usr/bin/env perl
use strict;
use warnings;

@ARGV == 6 or die "usage: $0 <bam> <reads1> <reads2> <samtools_range> <outprefix> <samtools>";

my $bam = $ARGV[0];
my $reads1 = $ARGV[1];
my $reads2 = $ARGV[2];
my $samtools_range = $ARGV[3];
my $outprefix = $ARGV[4];
my $samtools = $ARGV[5];
my $seed = 1;

-e $bam or die "bam $bam not found";
-e $reads1 or die "reads1 $reads1 not found";
-e $reads2 or die "reads2 $reads2 not found";

my $original_depth = `$samtools depth -a -r $samtools_range $bam | awk '{s+=\$3} END{print s/NR}'`;
chomp $original_depth;
print "Mean depth on region: $original_depth\n";
my @coverages = (1..25, 30, 35, 40, 45, 50, 75, 100);

foreach my $coverage (@coverages) {
    if ($coverage > $original_depth) {
        print "Skipping coverage $coverage because too low\n";
        next;
    }

    my $percent = 100 * $coverage / $original_depth;
    print "Making reads for coverage $coverage by taking $percent\% of the reads\n";
    my $p = "$outprefix.$coverage.reads";
    my $cmd = "fastaq to_random_subset --seed $seed --mate_file $reads2 $reads1 - $percent | fastaq deinterleave - $p\_1.fq.gz $p\_2.fq.gz";
    print "$cmd\n";
    system($cmd) and die "$!";
    $seed++;
}

