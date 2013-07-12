#!/usr/bin/perl 
# short-taxa-reads-bdiv-table.pl by nick bokulich 
# updated 6-19-13
# creates concatenated otu table for short-read-taxa-assignment comparisons from summarize_taxonomy files in multiple_assign_taxonomy.py subdirectories 
use strict; use warnings;

die "usage: short-taxa-reads-bdiv-table.pl <stem directory fp> <comma-separated list of assigners> <taxonomy level>\n" unless @ARGV == 3;
#Open input files

my $dir = $ARGV[0];
my $level = $ARGV[2];
my $list = $ARGV[1];
my @assignment_methods = split(",",$list);
my $methods_list = join("\t",@assignment_methods);

open(my $out, ">", "$dir/concatenated_otu_table.txt") or die "error reading $dir/concatenated_otu_table.txt for reading"; 
open(my $map, ">", "$dir/concatenated_otu_map.txt") or die "error reading $dir/concatenated_otu_map.txt for reading"; 
open(my $taxa_sum, ">", "$dir/concatenated_taxa_sum_L$level.txt") or die "error reading $dir/concatenated_taxa_sum_L$level.txt for reading"; 

print $out "# QIIME v1.5.0 OTU table\n#OTU ID\t$methods_list\ttaxonomy\n";

print $map "#SampleID\tassignment_method\tDescription\n";

print $taxa_sum "Taxon\t$methods_list\n";

#my @assignment_methods = ('blast_100.0', 'blast_0.0001', 'blast_1e-06', 'blast_1.0', 'blast_1e-10', 'blast_1e-30', 'mothur_0.0', 'mothur_0.1', 'mothur_0.2', 'mothur_0.3', 'mothur_0.4', 'mothur_0.5', 'mothur_0.6', 'mothur_0.7', 'mothur_0.8', 'mothur_0.9', 'mothur_1.0', 'rdp_0.0', 'rdp_0.1', 'rdp_0.2', 'rdp_0.3', 'rdp_0.4', 'rdp_0.5', 'rdp_0.6', 'rdp_0.7', 'rdp_0.8', 'rdp_0.9', 'rdp_1.0', 'rtax_single', 'Expected');

my %taxa;
my $rounds = 0;
my @samples;
my $path;

foreach my $assignment_method (@assignment_methods) { #read all subdirectories for sum_taxa files at specified level
	if(-e "$dir/$assignment_method/otu_table_mc2_w_taxa_L$level.txt"){$path = "$dir/$assignment_method/otu_table_mc2_w_taxa_L$level.txt"}
	else{$path = "$dir/$assignment_method/otu_table_mc2_no_pynast_failures_w_taxa_L$level.txt"}
	my ($method,$junk) = split("_",$assignment_method);
	print $map "$assignment_method\t$method\t$assignment_method\n";
	open(my $in, "<", "$path") or die "can't read open $path for reading: $!";
	$rounds++;
	my $header_line = <$in>;
	if ($rounds == 1) { # collect headers (sample names) from first file, use downstream to generate taxa listings for each sample
		chomp $header_line;
		@samples = split("\t",$header_line);
		shift @samples;
	}
	while (<$in>) { # read taxonomy lines, store abundance values for each sample/taxonomy
		chomp;
		my @abundances = split("\t",$_);
		my $taxonomy_string = shift (@abundances);
		if ($taxonomy_string =~ m/NO.*HIT/i) {$taxonomy_string = "NO_HIT"}
		for (my $i = 0; $i < @samples; $i++) { # for loop to loop both @samples and @abundances, add to hash matrix
			$taxa{$samples[$i]}{$taxonomy_string}{$assignment_method} = $abundances[$i];
		}
	}
	close $in;
}
close $map;

my $otu_id = 0;
foreach my $sample (keys %taxa) {
	foreach my $taxonomy_string (keys %{$taxa{$sample}}) {
		$otu_id++;
		print $out "$otu_id\t";
		print $taxa_sum "$sample;$taxonomy_string\t";
		foreach my $assignment_method (@assignment_methods) {
			if (exists $taxa{$sample}{$taxonomy_string}{$assignment_method}) {
				print $out "$taxa{$sample}{$taxonomy_string}{$assignment_method}\t";
				print $taxa_sum "$taxa{$sample}{$taxonomy_string}{$assignment_method}\t"
			}
			else {print $out "0\t";print $taxa_sum "0\t"}
		}
		print $out "$sample;$taxonomy_string\n";
		print $taxa_sum "\n";
	}
}
close $out;
