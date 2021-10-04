# 1. Getting data with _Fukomys damarensis_
Now, we have manually extracted `fukomys_damarensis` from the unfiltered alignments of 
rodentia ctenohystrica for each partition. Then, the same procedure followed in the previous 
steps will be followed now: we used MAFFT to align the new sequence to the previously generated 
alignment (which already includes the sequence for the squirrel). After that, we can get
the fasta sequences in one line:

```sh
# Run from `03_mafft/aln`
for i in `seq 2 6`
do

cd $i 
fasta=$( echo *_out.fasta )
printf "Parsing "$fasta"...\n"
../../../../../../../../../src/one_line_fasta.pl $fasta
grep '>' $fasta | sed 's/>//' > species.txt
mkdir out_mafft 
mv $fasta out_mafft
fa=$( ls *fa )
name=$( echo $fa | sed 's/\_out..*//' )
mv $fa $name".fasta"
cd ../

done
```

Then, check that the order is correct:

```sh
# Run from `03_mafft/aln`
diff 2/species.txt 3/species.txt 
diff 2/species.txt 4/species.txt 
diff 2/species.txt 5/species.txt 
diff 2/species.txt 6/species.txt 
```

As everything is correct, we can generate the `sequences.txt` file before concatenating the partitions:

```sh
# Run from `03_mafft/aln`
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
# Run from `03_mafft/aln`
mkdir 1
cd 1
paste -d "" ../2/cteno*fasta ../3/cteno*seq*txt > ctenohystrica.fasta 
paste -d "" ctenohystrica.fasta ../4/cteno*seq*txt > ctenohystrica2.fasta 
paste -d "" ctenohystrica2.fasta ../5/cteno*seq*txt > ctenohystrica3.fasta 
paste -d "" ctenohystrica3.fasta ../6/cteno*seq*txt > ctenohystrica4.fasta
rm ctenohystrica.fasta ctenohystrica2.fasta ctenohystrica3.fasta
mv ctenohystrica4.fasta ctenohystrica.fasta
```

Now, we convert the FASTA files into PHYLIP format:

```sh
# Run from `03_mafft/aln`
for i in `seq 1 6`
do 

cd $i
num=$( grep '>' cteno*fasta | wc -l )
len=$( sed -n '2,2p' *fasta | sed 's/\r//' | sed 's/\n//' | wc --m )
perl ../../FASTAtoPHYL.pl cteno*fasta $num $len 
cd ..

done
```

Last, generate the file with the five partitions concatenated 

```sh
mkdir 7 
cd 7 

for i in `seq 2 6`
do 

cat ../$i/*aln >> ctenohystrica_5parts.aln
printf "\n\n" >> ctenohystrica_5parts.aln

done 
```

**NOTE**: [Test 3](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_aln/extra_filtering/02_tests_fukomys/test3)
shows that there are issues with the CYTB sequence with *F. damarensis*. Therefore, the mitochondrial
alignments (mit-12CP and mit-3CP) generated at this stage are not going to be used. We need an extra 
filtering step, which we detail in
[`04_mafft_mitaln`](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation/rodentia_ctenohystrica/filter_aln/extra_filtering/04_mafft_mitaln) 

At this stage, we keep the nuclear and the mit-rna alignments only.