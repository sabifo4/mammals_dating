# 1. Generate fasta files to generate alignments with the squirrel
We use the perl script `phylip_to_fasta_and_extraseq.pl` to convert 
the ctenohystrica files from phylip format to fasta format. In addition, it appends 
the sequences from the taxa within the squirrel subtree specified in the file 
`taxa_to_add.txt`.

The following code can be run to generate the corresponding alignments for each 
of the partitions:

```sh
# Ron from this directory, `00_perl_parsing`
# Generate directories 
for i in `seq 2 6`
do
mkdir $i 
done 

# mt12cp
cd 2
perl ../phylip_to_fasta_and_extraseq.pl  ../../../../../../01_alignments/00_mammal_alns/rodentia_ctenohystrica/rodentia_ctenohystrica_mt12cp.aln \
../../../../../../01_alignments/00_mammal_alns/rodentia_squirrel/rodentia_squirrel_mt12cp.aln ../taxa_to_add.txt

# mt3cp
cd ../3 
perl ../phylip_to_fasta_and_extraseq.pl  ../../../../../../01_alignments/00_mammal_alns/rodentia_ctenohystrica/rodentia_ctenohystrica_mt3cp.aln \
../../../../../../01_alignments/00_mammal_alns/rodentia_squirrel/rodentia_squirrel_mt3cp.aln ../taxa_to_add.txt

# mtrna
cd ../4
perl ../phylip_to_fasta_and_extraseq.pl  ../../../../../../01_alignments/00_mammal_alns/rodentia_ctenohystrica/rodentia_ctenohystrica_mtrna.aln \
../../../../../../01_alignments/00_mammal_alns/rodentia_squirrel/rodentia_squirrel_mtrna.aln ../taxa_to_add.txt

# nt12cp
cd ../5
perl ../phylip_to_fasta_and_extraseq.pl  ../../../../../../01_alignments/00_mammal_alns/rodentia_ctenohystrica/rodentia_ctenohystrica_nt12cp.aln \
../../../../../../01_alignments/00_mammal_alns/rodentia_squirrel/rodentia_squirrel_nt12cp.aln ../taxa_to_add.txt

# nt3cp
cd ../6
perl ../phylip_to_fasta_and_extraseq.pl  ../../../../../../01_alignments/00_mammal_alns/rodentia_ctenohystrica/rodentia_ctenohystrica_nt3cp.aln \
../../../../../../01_alignments/00_mammal_alns/rodentia_squirrel/rodentia_squirrel_nt3cp.aln ../taxa_to_add.txt

```

After checking that all the `species.txt` files are the same, which means that 
the sequences will be in the same order in all fasta files, then 
we can generate the concatenated file in dir `1`:

```sh
# Concatenated -- NOTE: "sequences.txt" files need to be in UNIX format,
# otherwise `paste` does not work!
# Run from `00_perl_parsing`
mkdir 1
cd 1
paste -d "" ../2/sequences.txt ../3/sequences.txt > rodentia_ctenohystrica.txt 
paste -d "" rodentia_ctenohystrica.txt ../4/sequences.txt > rodentia_ctenohystrica2.txt 
paste -d "" rodentia_ctenohystrica2.txt ../5/sequences.txt > rodentia_ctenohystrica3.txt 
paste -d "" rodentia_ctenohystrica3.txt ../6/sequences.txt > rodentia_ctenohystrica4.txt
rm rodentia_ctenohystrica.txt rodentia_ctenohystrica2.txt rodentia_ctenohystrica3.txt
mv rodentia_ctenohystrica4.txt rodentia_ctenohystrica.txt

# Now, add the species in fasta format
perl ../concatenated_format.pl rodentia_ctenohystrica.txt ../2/species.txt
```