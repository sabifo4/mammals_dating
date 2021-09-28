# 1. Run RAxML
Command used to test coding mitochondrial sequences from PC:

```
# 1. First, fix tag names for RAxML to run
sed -i 's/\.1\_..*\;\_mitochondrial//g' fukomys_spp_cytb.phy

# 2. Run RAxML
raxmlHPC.exe -f a -m GTRGAMMA -p 12345 -# 100 -x 12345 -# 500 -s fukomys_spp_cytb.phy -o Heterocephalus -n cytb
```

# 2. Check phylogeny inferred
Then, to match the sequence names with the tags we had to use, we ran the following:

```
# 1. First, keep a copy of the tree output by RAxML 
#    with the fixed tag names and generate a file 
#    to go over with the long names 
cp RAxML_bestTree.mit RAxML_bestTree.cytb.tags
awk '{printf $1"\n"}' fukomys_spp_cytb_wholetags.phy > species.txt
sed -i 's/\_(..*//' species.txt 
sed -i 's/\:..*//' species.txt

# 2. Go over the `species.txt` file, extract the tag name, 
#    find it in the newick tree, and replace with long name
while read line
do
#echo $line
name=$( echo $line | sed 's/\.1..*//g' )
echo $name
sed -i 's/'"${name}"'/'"${line}"'/g' RAxML_bestTree.cytb.tags
done < species.txt
```

This last test has proven that the RefSeq available for _Fukomys damarensis_ was contaminated.
This is the reason for having shown a wrong placement within the `Rodentia ctenohystrica` subtree. 

The next steps will be to use one of the other CYTB sequences available for _F. damarensis_ (and used to  
infer the tree discussed above), re-align it to the `Rodentia ctenohystrica` alignment (the one to which
Ictidomys had already been re-aligned to), and 
then proceed with the Bayesian dating analysis.


# 3. Partition into 12CP and 3CP
We decided to include in the alignment the sequence tagged as `AY425857.1_Cryptomys_damarensis_isolate_165_cytochrome_b_(Cytb)_gene,_partial_cds;_mitochondrial`.
We extracted the sequence and generated a FASTA file. To partition this sequence into 12CP and 3CP,
we ran the follwing commands:

```
# Run from test3
perl partition12CPand3CP_v2.pl Cryptomis_damarensis_CYTB.fasta
```
