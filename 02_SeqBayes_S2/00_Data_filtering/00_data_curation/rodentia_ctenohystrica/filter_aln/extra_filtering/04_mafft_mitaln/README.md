# 1. Getting data
This time, we only ran MAFFT for the mitochondrial data subsets so we could re-align the 
CYTB sequence for _F. damarensis_ to the previously generated alignment (i.e., with the 
squirrel). Once MAFFT finished, we ran the following code to obtain the fasta sequences in one line:

```sh
# Run from 04_mafft/aln
for i in `seq 2 3`
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

Then, we used the file `species.txt` output for each of the alignments (one partition per alignment) to check 
that the order of the taxa was the same:

```sh
# Run from 04_mafft/aln
diff 2/species.txt 3/species.txt 
```

Now, we converted the FASTA files into PHYLIP format:

```sh
# Run from 04_mafft/aln
for i in `seq 2 3`
do 

cd $i
num=$( grep '>' cteno*fasta | wc -l )
len=$( sed -n '2,2p' cteno*fasta | sed 's/\r//' | sed 's/\n//' | wc --m )
perl ../../FASTAtoPHYL.pl cteno*fasta $num $len 
cd ..

done
```

Now, we manually change the name of the sequence for *F. damarensis* so the tag 
is "fukomys_damarensis" instead of "C.damarensis_AY425857".

Last, we generated the new concatenated alignment with these two mit new data as well 
as the alignment with 5 partitions:

```sh 
# Run from 04/mafft_mitaln/aln
for i in `seq 2 3`
do 

cd $i 
perl ../../only_sequences.pl *fasta
cd ..

done

# Now generate concatenated alignment 
mkdir 1
cd 1
paste -d "" ../2/cteno*fasta ../3/cteno*seq*txt > ctenohystrica.fasta 
paste -d "" ctenohystrica.fasta ../../../03_mafft/aln/4/cteno*seq*txt > ctenohystrica2.fasta 
paste -d "" ctenohystrica2.fasta ../../../03_mafft/aln/5/cteno*seq*txt > ctenohystrica3.fasta 
paste -d "" ctenohystrica3.fasta ../../../03_mafft/aln/6/cteno*seq*txt > ctenohystrica4.fasta
rm ctenohystrica.fasta ctenohystrica2.fasta ctenohystrica3.fasta
mv ctenohystrica4.fasta ctenohystrica.fasta

# Now get 5 partitions alignment
cd ..
mkdir 7
cd 7
for i in `seq 2 3`
do 
cat ../$i/*aln >> ctenohystrica_5parts.aln
printf "\n\n" >> ctenohystrica_5parts.aln
done 
for i in `seq 4 6` 
do
cat ../../../03_mafft/aln/$i/*aln >> ctenohystrica_5parts.aln
printf "\n\n" >> ctenohystrica_5parts.aln
done
```