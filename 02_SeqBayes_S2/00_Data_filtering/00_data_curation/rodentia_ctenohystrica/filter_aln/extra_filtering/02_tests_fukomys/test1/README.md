# 1. Get mit data
We have extracted the sequences for all Fukomys species and _Heterocephalus glaber_ from the 
sequences that are generated before applying any filtering step so we can identify 
what the issue with the sequences are (download from [here](https://www.dropbox.com/s/2lh9kg8i1xhhds5/SeqBayesS2_Raln_ctenohystrica_unfiltered.zip?dl=0)).
Now, we are going to run `RAxML` using all 14 species (using the mole rat as an outgroup).
The command used for this purpose will be the following:

```
# Run from `02_tests_fukomys` in both directories `mt3cp` and `mt12cp`
# HPC Apocrita used for this purpose
cd <dir>
seq=$( ls *aln )
# Note: `-x` turns on rapid bootstrapping (I set it for `-# 500` iterations)
raxmlHPC -f a -m GTRGAMMA -p 12345 -# 100 -x 12345 -# 500 -s $seq -o heterocephalus_glaber -n mt12CP
```

>> NOTE: WE had to delete four sequences, because they were all gaps (no mitochondrial genes for those!) 
>>       and then re-run again RAxML fotr `mt12cp` and `mt3cp`. These taxa were:
>>         `fukomys_anselli`, `fukomys_foxi`, `fukomys_kafuensis`, `fukomys_ochraceocinereus`

# 2. Get data from the rest of the partitions 
As the previous test showed _Fukomys damarensis_ clustering outside the cluster with the rest of Fukomys taxa, 
then we decided to test the same approach with the rest of the partitions: `nt12cp`, `nt3cp`, and `mtrna`.
Unfortunately, for the nuclear partitions there are only 3 taxa with available genomes. Therefore, RAxML cannot 
run with such a few amount of taxa and we only ran an additional test with the mit-tRNA partition.

# 3. Phylogeny with _Fukomys damarensis_
According to the `mtrna` partition, _Fukomys damarensis_ is part of the cluster with the rest of the 
Fukomys taxa. Nevertheless, the latter form a separate cluster from _Fukomys damarensis_. Besides, the 
branch length is so long that it seems that this species is more similar to the mole rat than to the 
rest of the Fukomys taxa. Therefore, we are going to download again the mit-tRNA sequences from the NCBI 
and re-run everything. This will discard any possible contamination error or any other issues during the
first filtering steps when collecting and filtering the sequences. See `test2` for more information on 
the subsequent test carried out.

