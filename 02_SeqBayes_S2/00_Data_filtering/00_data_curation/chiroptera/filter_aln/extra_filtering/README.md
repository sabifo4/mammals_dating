# Subsetting "Chiroptera"

## 1. Obtain subtrees
First, we used a copy of the 885sp rodents tree, which is saved under the 
[`00_R_parsing`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/00_R_parsing)
directory as `main_bats_uncalib.tree`, to find out how we had to divide it into two data subsets.
The first tree subset will contain "Megachiroptera", while the second subset 
will contain "Microchiroptera". 

The R script 
[`Generate_subtrees_bats.R`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/00_R_parsing/Generate_subtrees_bats.R), 
which is located in the
[`00_R_parsing`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/00_R_parsing)
directory, is used to find out which taxa have to be allocated to each data subset.
The input file is
[`main_bats_uncalib.tree`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/00_R_parsing/main_bats_uncalib.tree)
The output files are
[`ChiroSubt1_megachiroptera.tree`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/00_R_parsing/ChiroSubt1_megachiroptera.tree)
and
[`ChiroSubt2_microchiroptera.tree`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/00_R_parsing/ChiroSubt2_microchiroptera.tree).
We visually inspected 
so we could then manually generate the two calibrated subtrees that are saved in this same directory.

## 2. Generating alignments 
Two of the output files generated in the previous step by the R script, 
[`255sp_Lchiro_megachiroptera_list.txt`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/00_R_parsing/255sp_Lchiro_megachiroptera_list.txt)
and
[`633sp_Lchiro_microchiroptera_list.txt`](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/00_R_parsing/633sp_Lchiro_microchiroptera_list.txt)
files, are used as input files in this step to generate the data alignments for these two 
data subsets.

We used the perl script
[`Subset_alignments.pl`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/01_perl_parsing/Subset_alignments.pl)
saved in the
[`01_perl_parsing`](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/01_perl_parsing)
directory, for that purpose. Specifically, it reads a text file with a list of taxa names (first arg) and 
an alignment in PHYLIP format (second arg). The third argument is the number of taxa that 
are included in this file. 
The script will generate an output file with only those 
taxa that were listed in the text file, a summary file with the length of the 
alignment parsed, and the alignment in fasta file.

We first ran it for the first subset (Megachiroptera, 255 sp):

```sh
# Run this code from `01_perl_parsing` directory
for i in ../../../../../01_alignments/00_mammal_alns/chiroptera/*aln
do
if [[ ! $i =~ "5parts" ]]
then
perl Subset_alignments.pl ../00_R_parsing/255sp_Lchiro_megachiroptera_list.txt $i 255
fi
done
```

Please note that the alignment with 5 partitions needs to be generated again having 
the partitions in the correct order: mt12cp, mt3cp, mtrna, nt12cp, nt3cp.

```sh
# Run this code from `01_perl_parsing` directory
cat chiroptera_mt12cp_out.aln >> chiroptera_5parts.aln 
printf "\n\n" >> chiroptera_5parts.aln 
cat chiroptera_mt3cp_out.aln >> chiroptera_5parts.aln 
printf "\n\n" >> chiroptera_5parts.aln 
cat chiroptera_mtrna_out.aln >> chiroptera_5parts.aln 
printf "\n\n" >> chiroptera_5parts.aln 
cat chiroptera_nt12cp_out.aln >> chiroptera_5parts.aln 
printf "\n\n" >> chiroptera_5parts.aln 
printf "\n\n" >> chiroptera_5parts.aln 
cat chiroptera_nt12cp_out.aln >> chiroptera_5parts.aln
mv chiroptera_5parts.aln chiroptera_5parts_out.aln
```

Then, we moved the alignment subsets to the corresponding directory:

```sh
# Run this code from `01_perl_parsing` directory
for i in *aln
do
fasta=$( echo $i | sed 's/\.aln/\.fasta/' )
name=$( echo $i | sed 's/\_out..*//' | sed 's/chiroptera/chiroptera\_subt1/' )
name=$( echo $name | sed 's/\.aln//' )
mv $i $name".aln"
if [[ ! -d aln_subt1 ]]
then
mkdir aln_subt1
fi
mv $name".aln" aln_subt1
done 

## Move the rest of files 
if [[ ! -d fasta_subt1 ]] 
then 
mkdir fasta_subt1 
fi
mv *fasta fasta_subt1
if [[ ! -d log_subt1 ]] 
then 
mkdir log_subt1 
fi
mv *txt *csv log_subt1/
```

We repeated the same for the second subtree:

```sh
# Run this code from `01_perl_parsing` directory
# Get alignments 
for i in ../../../../../01_alignments/00_mammal_alns/chiroptera/*aln
do
if [[ ! $i =~ "5parts" ]]
then
perl Subset_alignments.pl ../00_R_parsing/633sp_Lchiro_microchiroptera_list.txt $i 633
fi
done

# Get 5 parts aln
cat chiroptera_mt12cp_out.aln >> chiroptera_5parts.aln 
printf "\n\n" >> chiroptera_5parts.aln 
cat chiroptera_mt3cp_out.aln >> chiroptera_5parts.aln 
printf "\n\n" >> chiroptera_5parts.aln 
cat chiroptera_mtrna_out.aln >> chiroptera_5parts.aln 
printf "\n\n" >> chiroptera_5parts.aln 
cat chiroptera_nt12cp_out.aln >> chiroptera_5parts.aln 
printf "\n\n" >> chiroptera_5parts.aln 
printf "\n\n" >> chiroptera_5parts.aln 
cat chiroptera_nt12cp_out.aln >> chiroptera_5parts.aln
mv chiroptera_5parts.aln chiroptera_5parts_out.aln

# Move alignments
for i in *aln
do
fasta=$( echo $i | sed 's/\.aln/\.fasta/' )
name=$( echo $i | sed 's/\_out..*//' | sed 's/chiroptera/chiroptera\_subt2/' )
name=$( echo $name | sed 's/\.aln//' )
mv $i $name".aln"
if [[ ! -d aln_subt2 ]]
then
mkdir aln_subt2
fi
mv $name".aln" aln_subt2
done

## Move the rest of files 
if [[ ! -d fasta_subt2 ]] 
then 
mkdir fasta_subt2 
fi
mv *fasta fasta_subt2
if [[ ! -d log_subt2 ]] 
then 
mkdir log_subt2 
fi
mv *txt *csv log_subt2/
```

## 3. Add new taxa
Now, we need to add extra taxa to both subtrees so we can include extra calibrations that will later 
ease grafting the subtrees to the backbone tree. The steps followed for the first rodentia subtree 
can be found 
[here](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/02_MAFFT_subt1),
while those to generate the second data subset are 
[here](/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/02_MAFFT_subt2).
