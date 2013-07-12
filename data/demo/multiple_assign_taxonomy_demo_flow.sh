wget ftp://thebeast.colorado.edu/pub/QIIME_DB_Public_Studies/study_1688_split_library_seqs_and_mapping.tgz

tar -zxf $HOME/study_1688_split_library_seqs_and_mapping.tgz
gunzip $HOME/study_1688_split_library_seqs_and_mapping/study_1688_split_library_seqs.fna.gz > $HOME/study_1688_split_library_seqs_and_mapping/study_1688_split_library_seqs.fna

### Uncomment the following TWO commands to run test analysis on subsampled dataset ***Note: will delete original raw sequences file. Back it up if you care
#subsample_fasta.py -i $HOME/study_1688_split_library_seqs_and_mapping/study_1688_split_library_seqs.fna -p 0.00005 -o $HOME/study_1688_split_library_seqs_and_mapping/study_1688_split_library_seqs_subsampled.fna
#mv $HOME/study_1688_split_library_seqs_and_mapping/study_1688_split_library_seqs_subsampled.fna $HOME/study_1688_split_library_seqs_and_mapping/study_1688_split_library_seqs.fna

#make sure correct parameters file is in path
#make sure to download gg_13_5 reference db, ensure paths are correct
pick_open_reference_otus.py --suppress_taxonomy_assignment -i $HOME/study_1688_split_library_seqs_and_mapping/study_103_split_library_seqs.fna -o $HOME/study_1688_split_library_seqs_and_mapping/ -r $HOME/gg_13_5/gg_13_5_otus/rep_set/97_otus.fasta -p $HOME/params/demo_params.txt  -f
echo "picked otus"

#note: need to modify ref db to remove spaces, otherwise mothur fails
sed 's/; /;/g' $HOME/gg_13_5/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt > $HOME/gg_13_5/gg_13_5_otus/taxonomy/97_otu_taxonomy_MOD.txt

multiple_assign_taxonomy.py --read_1_seqs_filename $HOME/study_1688_split_library_seqs_and_mapping/study_103_split_library_seqs.fna --rtax_amplicon_id_regexes '(\S+)\s(\S+?)\s' --rtax_header_id_regexes '\S+\s(\S+?)\s' -x single -f -t $HOME/gg_13_5/gg_13_5_otus/taxonomy/97_otu_taxonomy_MOD.txt -e 100,0.0001,0.000001,1,1e-10,1e-30 -c 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 --clean_otu_table_filename otu_table_mc2_no_pynast_failures.biom --rdp_max_memory 16000 -i $HOME/study_1688_split_library_seqs_and_mapping/ -o $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/ -m blast,rdp,mothur,rtax -r $HOME/gg_13_5/gg_13_5_otus/rep_set/97_otus.fasta
echo "multiple assign taxonomy"

#need key $HOME/mock-community-compositions/study_1688_split_library_seqs_and_mapping_key.txt
generate_taxa_compare_table.py -r $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/ -k $HOME/mock-community-compositions/ -o $HOME/study_1688_split_library_seqs_and_mapping/taxa_compare_table-13_5/ 
echo "taxa compare"

mkdir $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/Expected/
cp $HOME/mock-community-compositions/S16S_key.txt $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/Expected/otu_table_mc2_w_taxa_L6.txt

#make sure short-taxa-reads-bdiv-table.pl is in path
perl short-taxa-reads-bdiv-table.pl $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/  

convert_biom.py --biom_table_type="otu table" --process_obs_metadata taxonomy -i $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/concatenated_otu_table.txt -o $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/concatenated_otu_table.biom

beta_diversity_through_plots.py -m $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/concatenated_otu_map.txt -p $HOME/params/demo_params.txt --color_by_all_fields -f -i $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/concatenated_otu_table.biom -o $HOME/study_1688_split_library_seqs_and_mapping/concatenated_bdiv_plots/
echo "dbiv"

sed 's/^.*;//g' $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/concatenated_taxa_sum_L6.txt > $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/concatenated_taxa_sum_L6-genus-only.txt
make_emperor.py --biplot_fp $HOME/study_1688_split_library_seqs_and_mapping/concatenated_bdiv_plots/bray_curtis_3d_discrete_taxa/taxa_coordinates.txt -b "assignment_method,Description" -m $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/concatenated_otu_map.txt -n -1 -i $HOME/study_1688_split_library_seqs_and_mapping/concatenated_bdiv_plots/bray_curtis_pc.txt -o $HOME/study_1688_split_library_seqs_and_mapping/concatenated_bdiv_plots/bray_curtis_3d_discrete_taxa/ -t $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/concatenated_taxa_sum_L6-genus-only.txt
echo "emperor" 

make_bipartite_network.py --osizes 'Abundance' --ssizes 'Abundance' --scolors 'assignment_method' --ocolors 'p' --sshapes 'Description' --oshapes 'f,g' -k taxonomy --md_fields mock,k,p,c,o,f,g -m $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/concatenated_otu_map.txt -i $HOME/study_1688_split_library_seqs_and_mapping/multiple_assign_taxonomy_13_5/psrotot_gg_13_5/concatenated_otu_table.biom -o $HOME/study_1688_split_library_seqs_and_mapping/bipartite_network/
echo "network" 
