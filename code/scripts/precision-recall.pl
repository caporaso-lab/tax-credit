#!/usr/bin/perl 
# precision-recall.pl by nick bokulich 
# updated 7-19-13
# Computes precision/recall scores for series of input taxonomy summary files of a mock community
# Input consists of stem directory containing series of directories for analysis (multiple mock communities)
# These directies should each contain identical directories for comparison (different taxonomy assignment 
# configurations) AND an "Expected" directory containing a key of the expected taxonomic structure.

# Stem directory should be output of multiple_assign_taxonomy.py or follow that output directory structure

# Input directories should be all directories in that folder that should be read

# Assigners should be list of subdirectories containing taxonomy assignments in ALL input directories that should be read

# Output consists of a table of precision, recall, and f values for each input taxonomy assignment fp at each taxonomic level

use strict; use warnings;

die "usage: precision-recall.pl <stem directory fp> <comma-separated list of input directory names> <comma-separated list of assigners> <taxonomic levels>\n" unless @ARGV == 4;

#Set global parameters from inputs

my $dir = $ARGV[0];
my @input_directories = split(",",$ARGV[1]);
my @assignment_methods = split(",",$ARGV[2]);
my @levels = split(",",$ARGV[3]);
my $methods_list = join("\t",@assignment_methods);
my $top_level = pop(@levels);

#Open and print header lines for output table for this $level
open(my $out, ">", "$dir/precision_recall_table.txt") or die "error reading $dir/precision_recall_table.txt for reading"; 
#print $out "Precison-Recall Output Level $level\nValues: precision,recall,f,false positive,false negative\nInput\t$methods_list\n";
print $out "Precison-Recall Output \nValues: precision,recall,f-measure,true positive, false positive,false negative\nSet\tLevel_Method\tLevel\tMethod\tSample\tP\tR\tF\tTP\tFP\tFN\n";

foreach my $level (@levels) { #iterate through each assignment level

	my %taxa;
	my $path;
	my %key_taxa;
	
	foreach my $input_directory (@input_directories) { #read all subdirectories within each comma-separated input directory
		
		#Identify local variables for this input directory
		my @taxonomy_strings;

		#generate new expected values key for each level for each input directory except top level	
		#Open top level expected taxonomy key, extract taxonomy strings and generate new strings
		#***LIMITATION: Assumes that taxonomy string is found in all samples in expected taxonomy key
		open(my $expected_key, "<", "$dir/$input_directory/Expected/otu_table_mc2_w_taxa_L$top_level.txt") or die "error reading $dir/$input_directory/Expected/otu_table_mc2_w_taxa_L$top_level.txt: $!";
		my $head = <$expected_key>;
		while (<$expected_key>) {
			my @abundances = split("\t",$_);
			my $original_taxonomy_string = shift (@abundances);
			my @taxonomy_levels = split(";",$original_taxonomy_string);
			my @new_taxonomy_levels = @taxonomy_levels[0..$level-1];
			my $taxonomy_string = join (";",@new_taxonomy_levels);
			if ($taxonomy_string ~~ @taxonomy_strings) {}
			else {push (@taxonomy_strings, $taxonomy_string)}
		}
		close $expected_key;

		foreach my $assignment_method (@assignment_methods) { #read all subdirectories for sum_taxa files at specified level
			
			#Identify local variables for each assignment method
			my @samples;
			my %true_positive;
			my %false_positive;
			my %false_negative;
			
			if(-e "$dir/$input_directory/$assignment_method/otu_table_mc2_w_taxa_L$level.txt"){$path = "$dir/$input_directory/$assignment_method/otu_table_mc2_w_taxa_L$level.txt"}
			else{$path = "$dir/$input_directory/$assignment_method/otu_table_mc2_no_pynast_failures_w_taxa_L$level.txt"}
			open(my $in, "<", "$path") or die "can't read open $path for reading: $!";			
			my $header_line = <$in>; # collect headers (sample names) from first file, use downstream to generate taxa listings for each sample
			chomp $header_line;
			@samples = split("\t",$header_line);
			shift @samples;
			foreach my $sample (@samples) {
				$false_negative{$sample}=@taxonomy_strings;
				$false_positive{$sample}=0;
				$true_positive{$sample}=0;
			}
			while (<$in>) { # read taxonomy lines, store abundance values for each sample/taxonomy
				chomp;
				my @abundances = split("\t",$_);
				my $taxonomy_string = shift (@abundances);
				#set false negatives = size of expected taxonomy array
				if ($taxonomy_string ~~ @taxonomy_strings) { # Does $taxonomy_string match Expected? If no, add to $false_positive count
					for (my $i = 0; $i < @samples; $i++) { #  loop @samples and @abundances, match against Expected key and add to hash matrix
						if ($abundances[$i] == 0) { #if the string matches but abundance is 0, that's a false negative
							if (exists $false_negative{$samples[$i]}) {$false_negative{$samples[$i]}++}
							else					             {$false_negative{$samples[$i]} = 1}
						}
						else { #if it matches and is not 0, that's true positive, also subtract from false negative count
							if (exists $true_positive{$samples[$i]}) {$true_positive{$samples[$i]}++}
							else					            {$true_positive{$samples[$i]} = 1}
							if (exists $false_negative{$samples[$i]}) {$false_negative{$samples[$i]}--}
							else					             {$false_negative{$samples[$i]} = -1}
						}
					}
				}
				else { #if it didn't match, you have a false positive
					foreach my $sample (@samples) {
						if (exists $false_positive{$sample}) {$false_positive{$sample}++}
						else					             {$false_positive{$sample} = 1}
					}
				}
			}
			# Calculate P, R, F for each sample, print P, R, F, TP, FP, FN
			# p =  # true positive / (# true positive + false positive)
			# r = # positive / (# true positive + # false negative)
			# f = (2*p*r) / (p + r) 
			foreach my $sample (@samples) {
				my $precision;
				my $recall;
				my $f_measure;
				if ($true_positive{$sample} == 0) {$precision = 0;$recall = 0;$f_measure = 0}
				else {
					$precision = $true_positive{$sample} / ($true_positive{$sample} + $false_positive{$sample});
					$recall = $true_positive{$sample} / ($true_positive{$sample} + $false_negative{$sample});
					$f_measure = 2 * $precision * $recall / ($precision + $recall);
				}
				print $out "$level $assignment_method\t$level $assignment_method\t$level\t$assignment_method\t$sample\t$precision\t$recall\t$f_measure\t$true_positive{$sample}\t$false_positive{$sample}\t$false_negative{$sample}\n";
			}
			close $in;
		}
	}
}
close $out;
