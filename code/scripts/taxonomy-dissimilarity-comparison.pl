#!/usr/bin/perl 
# taxonomy-dissimilarity-comparison.pl by nick bokulich 
# updated 7-14-13
# creates concatenated otu table for taxonomy dissimilarity comparisons in qiime from summarize_taxonomy files in multiple_assign_taxonomy.py subdirectories 
use strict; use warnings;

die "usage: short-taxa-reads-bdiv-table.pl <stem directory fp> <comma-separated list of assigners> <taxonomy level>\n" unless @ARGV == 3;
#Open input files

my $dir = $ARGV[0];
my $level = $ARGV[2];
my $list = $ARGV[1];
my @assignment_methods = split(",",$list);

open(my $out, ">", "$dir/dissimilarity_otu_table.txt") or die "error reading $dir/dissimilarity_otu_table.txt for reading"; 
open(my $map, ">", "$dir/dissimilarity_otu_map.txt") or die "error reading $dir/dissimilarity_otu_map.txt for reading"; 

print $out "# QIIME v1.5.0 OTU table\n#OTU ID\t";

print $map "#SampleID\tassignment_method\toriginal_sample\tDescription\n";

#my @assignment_methods = ('blast_100.0', 'blast_0.0001', 'blast_1e-06', 'blast_1.0', 'blast_1e-10', 'blast_1e-30', 'mothur_0.0', 'mothur_0.1', 'mothur_0.2', 'mothur_0.3', 'mothur_0.4', 'mothur_0.5', 'mothur_0.6', 'mothur_0.7', 'mothur_0.8', 'mothur_0.9', 'mothur_1.0', 'rdp_0.0', 'rdp_0.1', 'rdp_0.2', 'rdp_0.3', 'rdp_0.4', 'rdp_0.5', 'rdp_0.6', 'rdp_0.7', 'rdp_0.8', 'rdp_0.9', 'rdp_1.0', 'rtax_single', 'Expected');

my %taxa;
my @samples;
my $path;

foreach my $assignment_method (@assignment_methods) { #read all subdirectories for sum_taxa files at specified level
	if(-e "$dir/$assignment_method/otu_table_mc2_w_taxa_L$level.txt"){$path = "$dir/$assignment_method/otu_table_mc2_w_taxa_L$level.txt"}
	else{$path = "$dir/$assignment_method/otu_table_mc2_no_pynast_failures_w_taxa_L$level.txt"}
	my ($method,$junk) = split("_",$assignment_method);
	open(my $in, "<", "$path") or die "can't read open $path for reading: $!";
	my $rounds = 0;
	$rounds++;
	my $header_line = <$in>;
	if ($rounds == 1) { # collect headers (sample names) from first file, use downstream to generate taxa listings for each sample
		chomp $header_line;
		@samples = split("\t",$header_line);
		shift @samples;
		foreach my $samples (@samples) {
			print $map "$assignment_method.$samples\t$method\t$samples\t$assignment_method\n";
			print $out "$assignment_method.$samples\t";
		}
	}
	while (<$in>) { # read taxonomy lines, store abundance values for each sample/taxonomy
		chomp;
		my @abundances = split("\t",$_);
		my $taxonomy_string = shift (@abundances);
		if ($taxonomy_string =~ m/NO.*HIT/i) {$taxonomy_string = "NO_HIT"}
		for (my $i = 0; $i < @samples; $i++) { # for loop to loop both @samples and @abundances, add to hash matrix
			$taxa{$taxonomy_string}{$assignment_method}{$samples[$i]} = $abundances[$i];
		}
	}
	close $in;
}
close $map;

print $out "taxonomy\n";

my $otu_id = 0;
foreach my $taxonomy_string (keys %taxa) {
	$otu_id++;
	print $out "$otu_id\t";
#	foreach my $assignment_method (keys %{$taxa{$taxonomy_string}}) {
	foreach my $assignment_method (@assignment_methods) {
#		foreach my $sample (keys %{$taxa{$taxonomy_string}{$assignment_method}}) {
		foreach my $sample (@samples) {
			if (exists $taxa{$taxonomy_string}{$assignment_method}{$sample}) {
				print $out "$taxa{$taxonomy_string}{$assignment_method}{$sample}\t";
			}
			else {print $out "0\t"}
		}
	}
	print $out "$taxonomy_string\n";
}
close $out;
