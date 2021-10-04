# 1. Get mit data
We downloaded the mit-tRNA sequences from the NCBI (click [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_027742.1)). The NCBI reference sequence 
was retrieved from the main site for the *Fukomys damarensis* genome (see [here](https://www.ncbi.nlm.nih.gov/genome/?term=fukomys)).

These sequences have been saved in separate files and are saved in the directory `NCBI_mtrna_fukomys_damarensis`.

# 2. Generate concatenated file with all mitrna sequences 
First, we concatenated the separate files in a unique file. Then, we used a perl script to generate 
a FASTA file with all these sequences concatenated:

```sh
# Run from `02_tests_fukomys/test2`
# 1. Get everything in one line
for i in NCBI_mtrna_fukomys_damarensis/*fasta
do
../../../../../../../../src/one_line_fasta.pl $i
done 
# 2. Generate file with all FASTA sequences
for i in NCBI_mtrna_fukomys_damarensis/*fa
do
cat $i >> fukomys_damarensis_mtrna.fasta 
done
# 3. Generate a concatenated file
perl concatenate_FASTA.pl fukomys_damarensis_mtrna.fasta 
```

Then, we manually removed "mitochondrion" from the tag name, so it reads `fukomys_damarensis`.

# 3. Re-align
Now, we need to re-align this sequence to the rest of fukomys and _Heterocephalus glaber_. For that purpose,
we will use MAFFT. 

We got the FASTA file from
`01_mafft/aln/4/` (you need to have downloaded this data before in a previous step).
Then, extracted only those sequences 
for all the Fukomys taxa and the mole rat. This was the input alignment for MAFFT. The mitrna sequence downloaded 
from the NCBI (see steps above) is the sequence that will be re-aligned to the input alignment.

Then, we need to reformat the output alignment by MAFFT. We can do this by running the following commands:

```sh
# Run from test2/MAFFT
# 1. Get everything in one line 
cd 4
fasta=$( echo *outMAFFT.fasta )
printf "Parsing "$fasta"...\n"
../../../../../../../../../../src/one_line_fasta.pl $fasta
grep '>' $fasta | sed 's/>//' > species.txt
mkdir out_mafft 
mv $fasta out_mafft
fa=$( ls *fa )
name=$( echo $fa | sed 's/\_out..*//' )
mv $fa $name".fasta"

# 2. Get PHYLIP format 
num=$( grep '>' *fasta | wc -l )
len=$( sed -n '2,2p' *fasta | sed 's/\r//' | sed 's/\n//' | wc --m )
perl ../../../../01_mafft/FASTAtoPHYL.pl *fasta $num $len 
```

# 4. Run RAxML 
Now, we can use this alignment to infer the best-scoring ML tree with RAxML:

```sh
cd <dir>
seq=$( ls *aln )
# Note: `-x` turns on rapid bootstrapping (I set it for `-# 500` iterations)
raxmlHPC -f a -m GTRGAMMA -p 12345 -# 100 -x 12345 -# 500 -s $seq -o heterocephalus_glaber -n mt12CP
```

# 5. Phylogeny with _Fukomys damarensis_
According to the `mtrna` partition, _Fukomys damarensis_ is part of the cluster with the rest of the 
Fukomys taxa. Nevertheless, the latter form a separate cluster from _Fukomys damarensis_.
