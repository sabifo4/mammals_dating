# Data collection 
A summary of the process followed to collect the data for this analysis is given below. 
We also describe the steps followed to filter these data. Please note that further details 
to filter the first data set can be found 
[here](https://github.com/sabifo4/mammals_dating/blob/main/01_SeqBayes_S1/00_Gene_filtering/README.md),
while additional filterings for the second data set can be found [here](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2/00_Data_filtering/00_data_curation) (you 
will need to access each directory in this link to go through the details followed to 
filter each data subset).

## Data Set 1: 72-genome alignment
We downloaded the set of one-to-one protein-coding orthologs for the 72 mammal genomes
available in Ensembl release 98 (accessed 2019/11/15) using `EnsemblBioMarts` ([Kinsella et al., 2011](https://pubmed.ncbi.nlm.nih.gov/21785142/))
Sequences that did not meet the following requirements were removed from further analysis:   

   * Present in both human and mouse.   
   * Not containing stop codons or gene/transcript mismatches   
   * Present in at least 10 species and at least 100 codons in length.   
   
This left a total of 15,569 orthologs, which we partitioned into two data blocks: 
(i) first and second codon positions (12CP) and third codon positions (3CP). For each
ortholog, an alignment was built with `PRANK` v140603 ([Loytynoja 2013](https://link.springer.com/protocol/10.1007%2F978-1-62703-646-7_10))
and the best-scoring maximum-likelihood (ML) trees were inferred with `RAxML` v8.2.10 ([Stamatakis 2014](https://pubmed.ncbi.nlm.nih.gov/24451623/)).
Note we used only the alignments with the 12CP-partition in the subsequent Bayesian molecular clock analyses.

We further filtered the data set using the estimated best-scoring ML trees for each gene
to identify those having a branch whose length was larger than 60% of the total tree length
(the sum of all branch lengths). The relative branch length test is useful to detect
misaligned or misidentified orthologs in the alignments ([Springer & Gatesy, 2017](https://www.tandfonline.com/doi/full/10.1080/14772000.2017.1401016),
[dos Reis et al., 2012](https://pubmed.ncbi.nlm.nih.gov/22628470/)),  which may result in unusually long branch lengths.
Let $b_{ij}$ be the $i$-th branch length for gene tree $j$ , and let $n$ be the number of
branches in the tree, then the relative branch length is

$r_{ij}=b_{ij}/\sum_{i=j}^{n}b_{ij}$ (1).

We identified 133 ortholog alignments associated with at least one relative branch length larger
than 60%. These ortholog alignments were removed from further analysis.

Then, we estimated the pairwise distance between each ortholog in *Mus musculus*
and *Homo sapiens* using the R function `ape::dist.dna` ([Paradis et al., 2004](https://pubmed.ncbi.nlm.nih.gov/14734327/)). 
There were 4 genes (ENSG00000132185, ENSG00000204544, ENSG00000120937, and ENSG00000236699)
for which the distances were returned as `NaN` or were larger than 0.75 for at least one of
the substitution models used (i.e., TN93, JC69, and raw). Furthermore, when we plotted the
percentage of the tree length inferred for each ortholog alignment versus the corresponding
largest branch length (also in percentage), we found an outlier (ENSG00000176973,
see figure below). 

<p align="center">
  <img width="400" height="400" src="https://github.com/sabifo4/mammals/blob/master/figs/check_logtreelengthVSlargestbl.png">
</p>

>>**Fig S1: Quality control check of the estimated relative branch lengths (%) for each gene**
>>**tree in comparison to the corresponding total tree length (log-scale)**. One of the criteria
>>we had to keep a gene alignment was that the relative branch lengths (x-axis) could not be
>>larger or equal to 60% of the total length of the tree (y-axis). This plot shows how most
>>of the gene trees have relative branch lengths that are 10-30% of the length of the tree,
>>while few gene trees have relative branch lengths very close to the threshold (60%).
>>Nevertheless, there is one outlier whose tree length in the log-scale is 4 units despite
>>the relative branch length not accounting for more than 60% of the tree; there are at
>>least two very long branch lengths in this gene tree. This plot is useful to visually
>>identify which gene trees should be deleted according to the total log tree length that
>>the relative branch lengths of each gene tree account for with regards to a specific threshold.

We removed these 5 ortholog alignments, resulting in 15,431 orthologous gene alignments.
Of those, 163 further orthologs were removed after construction of data set 2 (see below).
This resulted in 15,268 ortholog alignments (Table S1). 

>>**Table S1**. Summary statistics for the 72-taxa gene alignments after each filtering step.

| Filtering step:          | Raw    | Initial filtering | Relative branch test | Pairwise distances | HMMER (data set 2) |
|--------------------------|--------|-------------------|----------------------|--------------------|--------------------|
| Number of genes left     | 15,904 | 15,569            | 15,436               | 15,431             | 15,268             |
| Genes removed            | 0      | 335               | 133                  | 5                  | 163                |
| % Removed                | 0      | 2.11              | 0.85                 | 0.03               | 1.06               |
| Cumulative genes removed | 0      | 335               | 468                  | 473                | 636                |
| % Cumulative removed     | 0      | 2.11              | 2.94                 | 2.97               | 4.00               |

The 15,268 ortholog alignments were sorted from slow- to fast-evolving according to the pairwise-distance
estimates and grouped into four partitions with the same number of genes.
Each of the four partitions contained the concatenated 12CP of the orthologs for the partition (Table S2).

>>**Table S2**. Number of taxa and orthologs, alignment length, corresponding site pattern counts, and missing data for each data subset
>>and partition scheme.

| Data subset | No. taxa | No. orthologs | Alignment length (base pairs) | Site pattern counts | Missing taxa |
|-------------|----------|---------------|-------------------------------|---------------------|--------------|
| Partition 1 | 72       | 3,817         | 8,926,316                     | 3,613,711           | 60.17        |
| Partition 2 | 72       | 3,817         | 8,339,196                     | 2,941,508           | 50.78        |
| Partition 3 | 72       | 3,817         | 8,605,264                     | 2,416,624           | 48.49        |
| Partition 4 | 72       | 3,817         | 7,302,398                     | 1,521,894           | 45.92        |

The rationale for this partitioning strategy is as follows. In our previous works ([dos Reis et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29342307/),
[dos Reis et al., 2015](https://www.sciencedirect.com/science/article/pii/S096098221501177X)), 
we tested phylogenomic data partitioning according to locus rate, principal component analysis of relative
branch lengths, and amino acid composition at loci. However, those analyses showed no noticeable differences
in posterior time estimates across the partitioning strategies ([dos Reis et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29342307/),
[dos Reis et al., 2015](https://www.sciencedirect.com/science/article/pii/S096098221501177X)).
On the other hand, we have shown that uncertainty in time estimates is sensitive to the number of partitions
used ([dos Reis et al., 2012](https://pubmed.ncbi.nlm.nih.gov/22628470/),
[dos Reis et al., 2014](https://www.jstor.org/stable/43699777?seq=1#metadata_info_tab_contents)) 
with more partitions producing more precise estimates, but at the cost of additional computation time.
It appears four partitions give a reasonable trade-off between computational speed and precision of estimates.
For example, twenty partitions would produce slightly more precise estimates but at 5 times the
computational cost.

## Data Set 2: Alignments of 4,705 taxa
We downloaded 832 complete mammal mitochondrial genomes from NCBI RefSeq
(accessed: 2016/01/14). Twelve extinct and two redundant entries were removed, leaving 818
genomes. Twelve protein-coding genes (all but *ND6*) and two rRNA non-coding genes 
(*RNR1* & *RNR2*) were extracted from each genome. The overlapping region in *ATP8* 
(position 95 to end) and overlapping codons at the end of *ND4L* were deleted.

To increase the nuclear and mitochondrial data sets, we mined sequences deposited 
in the European Nucleotide Archive (ENA, https://www.ebi.ac.uk/ena). The GenBank taxonomy 
(ftp://ftp.ncbi.nih.gov/pub/taxonomy/) was used to search for non-Ensembl mammalia species 
(this taxonomy is only used for ENA searches and not for any other analyses).
A total of 7,188 taxa were found, 83 of which were extinct. The GenBank identifiers were 
used to reference the corresponding taxa in the ENA, from which all matching coding and 
non-coding sequences for non-Ensembl mammal taxa were downloaded (accessed: 2016/01/17): we 
found 6,453 taxa with coding sequences and 3,239 taxa with non-coding sequences 
(6,606 distinct taxa).

This project started in early 2016. At the time, we downloaded 15K nuclear orthologous gene
alignments for 43 mammal taxa from Ensembl 83 and used these orthologs to create HMM sequence
profiles with `HMMER` ([Eddy 1998](https://pubmed.ncbi.nlm.nih.gov/9918945/)). The HMM profiles 
were then used to identify orthologs for additional taxa from GenBank, bypassing unreliable GenBank
homology annotations, and thus allowing reliable construction of large mammal subtrees.
In late 2019, we updated the 15K orthologs to 72 genomes available in Ensembl 98 (see above),
but the HMM profiles and corresponding homology searches are based on the 2016 mining of Ensembl.
HMM profiles were also created for mitochondrial protein-coding and non-coding genes and used for
taxa extension of the corresponding alignments. DNA homology searches were performed with 
`nhmmer` ([Wheeler & Eddy, 2013](https://pubmed.ncbi.nlm.nih.gov/23842809/)). 
using the following match criteria:   

   * Sequences with E-value $<1\exp^{-100}$ for a single gene were collected (i.e. sequences
   with multiple low E-values for different genes were removed).   
   * Matched sequences had to be at least 70% as long as the shortest Ensembl sequence in the
   alignment because many deposited sequences are partial sequences.   
   * Matches from hybrid/cross species were removed.   
   * Unspecified species (sp.) were excluded, unless no other member of the genus was represented
   (4 taxa included).   
   * Unconfirmed species (cf.) were excluded, unless no other member of the genus was represented
   (1 taxa included).   
   * Coding sequences were checked for correct open reading frame and translation.   
   
Nuclear genes resulting in an expanded set of at least 50 taxa were selected, resulting in 
a set of 168 nuclear genes. These 168 genes correspond to 163 genes in the 2019 Ensembl 
mining (5 genes did not pass filtering criteria for the 72 taxa, but they did pass the 
criteria with the 43 taxa). Thus, data set 1, based on the 2019 Ensembl mining, was reduced 
from 15,431 genes to 15,268 (Table S1) to avoid data duplication in the sequential dating 
approach. For new mined taxa, sequence annotations were extracted, sorted, and visually 
inspected to help verify homology. Alignments were then extended with homology-matched 
sequences using `PAGAN` v0.619 ([Loytynoja et al., 2012](https://pubmed.ncbi.nlm.nih.gov/22531217/)).
Sequences were added in order of decreasing length (i.e., 
longest sequences were added to the alignment first). Table S3 gives summaries of the numbers of
taxa and alignment lengths for data sets 1 and 2.

>>**Table S3**. Content of the mammal data set after bioinformatics filtering.

| Data set   | Data source            | Data type                     | Data set extended?   | Range of taxa | No. genes |
|------------|------------------------|-------------------------------|----------------------|---------------|-----------|
| Data set 1 | Ensembl                | Nuclear protein-coding        | No                   | 10-72         | 15,268    |
| Data set 2 | ENA/HMM-profile RefSeq | Nuclear protein-coding        | Yes, homology search | 50-1,266      | 168       |
| Data set 2 | RefSeq                 | Mitochondrial protein-coding  | No                   | 818           | 3         |
| Data set 2 | ENA / HMM-profile      | Mitochondrial protein-coding  | Yes, homology search | 1,047-4,131   | 9         |
| Data set 2 | ENA / HMM-profile      | Mitochondrial non-coding rRNA | Yes, homology search | 1,573-2,175   | 2         |

We then used `RAxML` to estimate the topology for each one of the 182 loci under the GTR+Gamma 
model. We then manually inspected the trees and further filtered taxa following these 
criteria:   

   * Remove taxa that did not share genes with their order, family, and genus. This is done to
   avoid unidentifiable positioning of taxa in the subtrees: if a species does not share genes with
   its close relatives, then several positionings of the species within the subtree will have the
   same likelihood (a.k.a. “likelihood terraces”).   
   * Keep only one member of each species while maintaining maximum locus coverage.
   That is, remove redundant subspecies. Many subspecies slow the analysis down and are not
   informative about deep divergences (e.g., *Rangifer tarandus tarandus*). Also, subspecies
   annotations are missing for many loci, leading to integrity problems when resolving tips.   
   * Outdated taxonomic names according to the literature were removed.   
   * Remove taxonomically mismatched or mislabelled taxa.   
   * Flag taxa with topological placement in mismatch with the literature.   
   * Outliers with unusually long branches in estimated trees were removed (three sequences
   in two genes).   

Taxa were separated according to the following taxonomic groups: Afrotheria, Xenarthra,
Marsupialia, Euarchonta, Lagomorpha, Laurasiatheria, and Rodentia. Laurasiatheria,
Rodentia and Chiroptera, which are species-rich, were further divided into additional
subsets to speed-up the dating analysis. Monotremata was added as an outgroup to all subtrees.
Thus, the final dataset has 4,705 taxa and 182 loci divided into 13 subtree alignments.
Each alignment was divided into five-partitions: (i) mitochondrial 12CP, (ii) mitochondrial 3CP,
(iii) mitochondrial RNA, (iv) nuclear 12CP, (v) nuclear 3CP (Tables S4-I-III, and S5).

>>**Table S4 (I)**. Number of taxa (alignment length | site pattern counts) for
>>mitochondrial subtrees and partitions. Monotremata is an outgroup in all subtrees.

| Data subset            | mit-12CP             | mit-3CP              | mit-RNA              |
|------------------------|----------------------|----------------------|----------------------|
| Afrotheria             | 34 (7,226 \| 2,298)  | 34 (3,613 \| 3,352)  | 32 (2,786 \| 1,486)  |
| Euarchonta             | 452 (7,260 \| 4,549) | 452 (3,630 \| 3,639) | 289 (3,087 \| 2,274) |
| Lagomorpha             | 82 (7,198 \| 1,326)  | 82 (3,599 \| 3,034)  | 46 (2,634 \| 928)    |
| Artiodactyla           | 419 (7,240 \| 3,430) | 419 (3,620 \| 3,614) | 293 (2,817 \| 1,828) |
| Chiroptera (I)         | 216 (7,202 \| 1,932) | 216 (3,601 \| 3,457) | 133 (3,287 \| 1,993) |
| Chiroptera (II)        | 566 (7,240 \| 3,020) | 566 (3,620 \| 3,590) | 328 (4,117 \| 3,033) |
| Rest of Laurasiatheria | 598 (7,320 \| 3,710) | 598 (3,660 \| 3,640) | 386 (3,402 \| 2,436) |
| Marsupialia            | 260 (7,278 \| 3,131) | 260 (3,639 \| 3,610) | 215 (3,366 \| 2,342) |
| Ctenohystrica          | 174 (7,224 \| 2,535) | 174 (3,612 \| 3,375) | 125 (2,897 \| 1,709) |
| Sciuridae and related  | 215 (7,200 \| 2,068) | 215 (3,600 \| 3,398) | 100 (2,958 \| 1,727) |
| Rest of Rodentia (I)   | 602 (7,366 \| 3,712) | 602 (3,694 \| 3,663) | 167 (3,072 \| 2053)  |
| Rest of Rodentia (II)  | 636 (7,320 \| 2,826) | 636 (3,662 \| 3,568) | 127 (2,927 \| 1,581) |
| Xenarthra              | 32 (7,210 \| 1,980)  | 32 (3,605 \| 3,404)  | 33 (2,655 \| 1,190)  |


>>**Table S4 (II)**. Number of taxa (alignment length | site pattern counts) for
>>nuclear subtrees and partitions. Monotremata is an outgroup in all subtrees.

| Data subset            | nuc-12CP                | nuc-3CP                 |
|------------------------|-------------------------|-------------------------|
| Afrotheria             | 52 (166,190 \| 6,853)   | 52 (83,095 \| 6,786)    |
| Euarchonta             | 253 (193,708 \| 40,599) | 253 (96,854 \| 42,084)  |
| Lagomorpha             | 43 (130,746 \| 1,834)   | 43 (65,373 \| ,2070)    |
| Artiodactyla           | 189 (191,394 \| 19,451) | 189 (95,697 \| 20,472)  |
| Chiroptera (I)         | 163 (222,499 \| 8,057)  | 163 (80,524 \| 9,316)   |
| Chiroptera (II)        | 448 (234,986 \| 9,407)  | 448 (117,248 \| 10,205) |
| Rest of Laurasiatheria | 453 (198,006 \| 45,479) | 453 (177,439 \| 43,565) |
| Marsupialia            | 249 (171,898 \| 10,417) | 249 (85,949 \| 9,322)   |
| Ctenohystrica          | 116 (3,327 \| 8,980)    | 116 (16,778 \| 7,043)   |
| Sciuridae and related  | 118 (24,890 \| 2,967)   | 118 (12,445 \| 3,115)   |
| Rest of Rodentia (I)   | 379 (423,700 \| 11,881) | 379 (94,825 \| 12,671)  |
| Rest of Rodentia (II)  | 505 (206,736 \| 8,285)  | 505 (95,950 \| 7,878)   |
| Xenarthra              | 20 (102,226 \| 2,513)   | 20 (51,113 \| 2,487)    |


>>**Table S4 (III)**. Total number taxa and alignment length (across all partitions) for each subtree.
>>Monotremata is an outgroup in all subtrees.

| Data subset            | Number of taxa | Total alignment length |
|------------------------|----------------|------------------------|
| Afrotheria             | 60             | 262,910                |
| Euarchonta             | 486            | 304,539                |
| Lagomorpha             | 88             | 209,550                |
| Artiodactyla           | 431            | 300,768                |
| Chiroptera (I)         | 256            | 317,113                |
| Chiroptera (II)        | 634            | 367,211                |
| Rest of Laurasiatheria | 659            | 389,827                |
| Marsupialia            | 307            | 272,130                |
| Ctenohystrica          | 210            | 63,798                 |
| Sciuridae and related  | 267            | 51,093                 |
| Rest of Rodentia (I)   | 630            | 532,657                |
| Rest of Rodentia (II)  | 691            | 316,595                |
| Xenarthra              | 33             | 166,809                |


>>**Table S5**. Missing data (%) for each subtree partition before (top) and after (bottom) removing
>>missing taxa. The percentage of missing data is calculated by dividing the number of gaps in
>>the alignment by the number of taxa times the alignment length. The number of taxa is shown within
>>brackets. We note missing taxa are not used by MCMCtree during likelihood calculation.

| Data subset                      | mit-12CP    | mit-3CP     | mit-RNA     | nuc-12CP    | nuc-3CP     |
|----------------------------------|-------------|-------------|-------------|-------------|-------------|
| Afrotheria (before)              | 70.66 (60)  | 70.66 (60)  | 55.44 (60)  | 93.16 (60)  | 93.16 (60)  |
| Afrotheria (after)               | 48.22 (34)  | 48.22 (34)  | 16.46 (32)  | 92.11 (52)  | 92.11 (52)  |
| Euarchonta (before)              | 51.51 (486) | 51.51 (486) | 58.28 (486) | 96.02 (486) | 96.02 (486) |
| Euarchonta (after)               | 47.86 (452) | 47.86 (452) | 29.84 (289) | 92.35 (253) | 92.35 (253) |
| Lagomorpha (before)              | 65.24 (88)  | 65.24 (88)  | 68.55 (88)  | 96.57 (88)  | 96.57 (88)  |
| Lagomorpha (after)               | 62.69 (82)  | 62.69 (82)  | 39.84 (46)  | 92.98 (43)  | 92.98 (43)  |
| Artiodactyla (before)            | 35.83 (431) | 35.83 (431) | 43.86 (431) | 97.85 (431) | 97.85 (431) |
| Artiodactyla (after)             | 33.99 (419) | 33.99 (419) | 17.42 (293) | 95.10 (189) | 95.10 (189) |
| Chiroptera I (before)            | 82.48 (256) | 82.48 (256) | 69.59 (256) | 98.33 (256) | 97.69 (256) |
| Chiroptera I (after)             | 79.23 (216) | 79.23 (216) | 41.47 (133) | 97.37 (163) | 96.37 (163) |
| Chiroptera II (before)           | 82.84 (634) | 82.84 (634) | 71.68 (634) | 99.16 (634) | 99.16 (634) |
| Chiroptera II (after)            | 80.78 (566) | 80.78 (566) | 45.26 (328) | 98.81 (448) | 98.81 (448) |
| Rest of Laurasiatheria (before)  | 59.35 (659) | 59.35 (659) | 66.98 (659) | 97.84 (659) | 98.80 (659) |
| Rest of Laurasiatheria (after)   | 55.20 (598) | 55.20 (598) | 43.62 (386) | 96.86 (453) | 98.25 (453) |
| Marsupialia (before)             | 61.52 (307) | 61.52 (307) | 57.51 (307) | 97.39 (307) | 97.39 (307) |
| Marsupialia (after)              | 54.56 (260) | 54.56 (260) | 39.33 (215) | 96.79 (249) | 96.79 (249) |
| Ctenohystrica (before)           | 83.33 (210) | 83.33 (210) | 76.03 (210) | 93.96 (210) | 94.01 (210) |
| Ctenohystrica (after)            | 79.88 (174) | 79.88 (174) | 59.74 (125) | 89.06 (116) | 89.15 (116) |
| Sciuridae and related (before)   | 84.83 (267) | 84.83 (267) | 85.37 (267) | 96.57 (267) | 96.57 (267) |
| Sciuridae and related (after)    | 81.16 (215) | 81.16 (215) | 60.95 (100) | 92.24 (118) | 92.24 (118) |
| Rest of Rodentia I (before)      | 81.15 (630) | 81.21 (630) | 87.69 (630) | 99.51 (630) | 98.90 (630) |
| Rest of Rodentia I (after)       | 80.27 (602) | 80.33 (602) | 53.54 (167) | 99.18 (379) | 98.16 (379) |
| Rest of Rodentia II (before)     | 84.02 (691) | 84.03 (691) | 90.69 (691) | 98.96 (691) | 98.88 (691) |
| Rest of Rodentia II (after)      | 82.64 (636) | 82.65 (636) | 49.35 (127) | 98.57 (505) | 98.46 (505) |
| Xenarthra (before)               | 3.30 (33)   | 3.30 (33)   | 4.67 (33)   | 90.04 (33)  | 90.04 (33)  |
| Xenarthra (after)                | 0.28 (32)   | 0.28 (32)   | 4.67 (33)   | 83.57 (20)  | 83.57 (20)  |