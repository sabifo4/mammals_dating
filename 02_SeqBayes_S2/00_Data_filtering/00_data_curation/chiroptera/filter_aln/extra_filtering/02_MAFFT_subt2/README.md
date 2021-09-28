# 1. Generate fasta files with new taxa
We use the perl script
[`phylip_to_fasta_and_extraseq.pl`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/02_MAFFT_subt2/01_perl_parsing/phylip_to_fasta_and_extraseq.pl)
to convert the files for the second chiroptera subtree from phylip format to fasta format. In addition, this
script appends the sequences from the taxa within the first chiroptera subtree specified in the file 
[`taxa_to_add.txt`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/chiroptera/filter_aln/extra_filtering/02_MAFFT_subt1/01_perl_parsing/taxa_to_add.txt).

The following code can be run to format the corresponding alignments for each 
of the partitions:

```sh
# Run from `01_perl_parsing`
# Generate directories
for i in `seq 2 6`
do
mkdir $i 
done 

# mt12cp
cd 2
perl ../phylip_to_fasta_and_extraseq.pl  ../../../01_perl_parsing/aln_subt2/chiroptera_subt2_mt12cp.aln \
../../../01_perl_parsing/aln_subt1/chiroptera_subt1_mt12cp.aln ../taxa_to_add.txt

# mt3cp
cd ../3 
perl ../phylip_to_fasta_and_extraseq.pl  ../../../01_perl_parsing/aln_subt2/chiroptera_subt2_mt3cp.aln \
../../../01_perl_parsing/aln_subt1/chiroptera_subt1_mt3cp.aln ../taxa_to_add.txt

# mtrna
cd ../4
perl ../phylip_to_fasta_and_extraseq.pl  ../../../01_perl_parsing/aln_subt2/chiroptera_subt2_mtrna.aln \
../../../01_perl_parsing/aln_subt1/chiroptera_subt1_mtrna.aln ../taxa_to_add.txt

# nt12cp
cd ../5
perl ../phylip_to_fasta_and_extraseq.pl  ../../../01_perl_parsing/aln_subt2/chiroptera_subt2_nt12cp.aln \
../../../01_perl_parsing/aln_subt1/chiroptera_subt1_nt12cp.aln ../taxa_to_add.txt

# nt3cp
cd ../6
perl ../phylip_to_fasta_and_extraseq.pl  ../../../01_perl_parsing/aln_subt2/chiroptera_subt2_nt3cp.aln \
../../../01_perl_parsing/aln_subt1/chiroptera_subt1_nt3cp.aln ../taxa_to_add.txt

```

After checking that all the `species.txt` files are the same, which means that 
the sequences will be in the same order in all fasta files, then 
we can generate the concatenated file in dir `1`. Also, make sure that all 
the `sequences.txt` files are in UNIX format, or there will be issues with 
the `\r` or other trail characters (all files had to be converted!):

```sh
# Run from `01_perl_parsing`
mkdir 1
cd 1
paste -d "" ../2/sequences.txt ../3/sequences.txt > chiroptera_subt2.txt 
paste -d "" chiroptera_subt2.txt ../4/sequences.txt > chiroptera_subt2_2.txt    
paste -d "" chiroptera_subt2_2.txt ../5/sequences.txt > chiroptera_subt2_3.txt 
paste -d "" chiroptera_subt2_3.txt ../6/sequences.txt > chiroptera_subt2_4.txt
rm chiroptera_subt2.txt chiroptera_subt2_2.txt chiroptera_subt2_3.txt
mv chiroptera_subt2_4.txt chiroptera_subt2.txt

# Now, add the species in fasta format
perl ../concatenated_format.pl chiroptera_subt2.txt ../2/species.txt
```

# 2. Getting alignment with new taxa 
We copied the files in the directories `01_perl_parsing/[2-6]/forMAFFT` to directories 
labelled from 2 to 6 to run MAFFT in the Apocrita HPC.
Then, we downloaded the alignments generated by MAFFT from the cluster.

**NOTE**: You will see that the files now start with "lchiroptera_subt1". This is because 
we changed the file names when we were running `MAFFT` analyses in the cluster, but the 
files correspond to the first chiroptera subtree.

```sh
mafft --add $new_seq $aln_name >$out_name
```

To get the fasta sequences in one line, I ran the following code:

```sh
## Run from the directory `mafft_lchiro2`
for i in `seq 2 6`
do

cd $i 
fasta=$( echo *outforMAFFT_out.fasta )
printf "Parsing "$fasta"...\n"
../../../../../../../../../../src/one_line_fasta.pl $fasta
grep '>' $fasta | sed 's/>//' > species.txt
mkdir out_mafft 
mv $fasta out_mafft
fa=$( ls *fa )
name=$( echo $fa | sed 's/\_out..*//' )
mv $fa $name".fasta"
cd ../

done
```

Then, we used the file `species.txt` output for each of the alignments (one partition per alignment) to check 
that the order of the taxa was the same:

```sh
# Run from the directory `02_mafft/mafft_lchiro2`
diff 2/species.txt 3/species.txt 
diff 2/species.txt 4/species.txt 
diff 2/species.txt 5/species.txt 
diff 2/species.txt 6/species.txt 
```

Once we made sure that the order was the same,
we had to generate the `sequences.txt` file alone without the `>species` header, otherwise the header 
would have been duplicated:

```sh
# Run from the directory `02_mafft/mafft_lchiro2`
for i in `seq 2 6`
do 

cd $i 
perl ../../only_sequences.pl *fasta
cd ..

done
```

Now, we can concatenate the sequences in the same order that had been previously done 
(i.e., mt12cp > mt3cp > mtrna > nt12cp > nt3cp):

```sh
# Run `02_mafft/mafft_lchiro2`
mkdir 1
cd 1
paste -d "" ../2/*fasta ../3/lchiro*seq*txt > lchiroptera_subt2.fasta 
paste -d "" lchiroptera_subt2.fasta ../4/lchiro*seq*txt > lchiroptera_subt2_2.fasta 
paste -d "" lchiroptera_subt2_2.fasta ../5/lchiro*seq*txt > lchiroptera_subt2_3.fasta 
paste -d "" lchiroptera_subt2_3.fasta ../6/lchiro*seq*txt > lchiroptera_subt2_4.fasta
rm lchiroptera_subt2.fasta lchiroptera_subt2_2.fasta lchiroptera_subt2_3.fasta
mv lchiroptera_subt2_4.fasta lchiroptera_subt2.fasta
```

Now, the last step is to convert the FASTA files into PHYLIP format:

```sh
# Run `02_mafft/mafft_lchiro2`
for i in `seq 1 6`
do 

cd $i
num=$( grep '>' *fasta | wc -l )
len=$( sed -n '2,2p' *fasta | sed 's/\r//' | sed 's/\n//' | wc --m )
perl ../../FASTAtoPHYL.pl *fasta $num $len 
cd ..

done
```

Last, generate the file with the five partitions concatenated:

```sh
# Run `02_mafft/mafft_lchiro2`
mkdir 7 
cd 7 

for i in `seq 2 6`
do 

cat ../$i/*aln >> lchiroptera_subt2_5parts.aln
printf "\n\n" >> lchiroptera_subt2_5parts.aln

done 
```
