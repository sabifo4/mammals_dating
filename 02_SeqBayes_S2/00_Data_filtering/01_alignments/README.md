# Alignments 
Here, you have the R scripts that have been used to filter the different data subsets and 
generate the final alignments and a directory with the log files that were output during 
this procedure.

In addition, you have a directory with the corresponding 
"dummy alignments" that were generated during the filtering step for each data subset, which 
can be used to save disk memory when running `MCMCtree` with the approximate likelihood.

Last, you have another directory with different files that summarise the missing data in each 
partition of each data subset. Below, you will find the code you need to run in order to calculate 
the missing data in each data subset, which uses the perl script
[`count_missingdat.pl`](/02_SeqBayes_S2/00_Data_filtering/01_alignments/count_missingdat.pl). 
Note this code will go through each individual partition of each data subset to generate the 
summary files with missing data:

```sh
# First, create dir to keep all output files 
# Run from `01_alignments`
mkdir count_missing_dat

# Now, run perl script and code to summarise results 
# Run from `01_alignments`
curr_dir=$( pwd )
for j in afrotheria chiroptera_subt1 chiroptera_subt2 euarchonta lagomorpha laurasiatheria_cetartiodactyla laurasiatheria_therest marsupialia rodentia_ctenohystrica rodentia_squirrel rodentia_subt1 rodentia_subt2 xenarthra 
do
cd $curr_dir/00_mammal_alns/$j
for i in *mt*cp.aln *nt*cp.aln
do 
name_aln=$i
name=$( echo $i | sed 's/\.aln//' )
mkdir -p ../../count_missing_dat/$name
printf "Parsing alignment "$name_aln"... ...\n"
printf "Parsing alignment "$name_aln"... ...\n" > ../../count_missing_dat/$name/log_count_missdat_$name".txt"
perl ../../count_missingdat.pl $name_aln ../../count_missing_dat/$name | tee ../../count_missing_dat/$name/log_count_missdat_$name".txt"
# Get a summary
name_dir=$( pwd | sed 's/..*\///' )
printf "<< DIR "$name_dir"/"$i" >>\n" >> ../../count_missing_dat/$name/$name_dir"_countNA.tsv"
sed -n '1,2p' ../../count_missing_dat/$name/out_count_NA/*avgmissdata.txt >> ../../count_missing_dat/$name/$name_dir"_countNA.tsv"
sed -n '7,7p' ../../count_missing_dat/$name/out_count_NA/*avgmissdata.txt >> ../../count_missing_dat/$name/$name_dir"_countNA.tsv"
sed -n '9,9p' ../../count_missing_dat/$name/out_count_NA/*avgmissdata.txt >> ../../count_missing_dat/$name/$name_dir"_countNA.tsv"
printf "\n" >> ../../count_missing_dat/$name/$name_dir"_countNA.tsv"
printf "\n\n"
done
done
```