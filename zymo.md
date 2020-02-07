# Zymo Mock Community Analysis 

The Zymo mock community was analyzed using three methods:
- [Mapping LoopSeq 16S contigs and Illumina V3V4 short reads to provided references](https://github.com/nico-chung/loopseq/blob/master/zymo.md#mapping-to-provided-references)
- [Mapping LoopSeq 16S contigs to NCBI RefSeq references](https://github.com/nico-chung/loopseq/blob/master/zymo.md#mapping-to-ncbi-refseq-references)
- [Classifying LoopSeq 16S contigs with LCA classifier](https://github.com/nico-chung/loopseq/blob/master/zymo.md#classification-with-lca-classifier)

## Experimental Protocol

[Fill in experimental details here: Sample prep, kits used, platforms used, yield per sample, etc.]


## Data Preparation

### Mock Community Description

The ZymoBIOMICS Microbial Community DNA Standard is composed of ten microbial species: eight bacteria and two fungi, the latter two of which were not amplified by the LoopSeq 16S rRNA primers used. Genomic DNA quantity from the eight species are evenly distributed and vary in Gram stain, genome size, GC content, 16S copy number and phylogenetic relatedness. 

| **Species** | **16S copies** | **Expected contig %** |
|---------|------------|-------------------|
Bacillus subtilis |	10 |	17.4
Enterococcus faecalis	| 4 | 9.9
Escherichia coli | 7	| 10.1
Lactobacillus fermentum	| 5 |	18.4
Listeria monocytogenes	| 6 |	14.1
Pseudomonas aeruginosa	| 4 |	4.2
Salmonella enterica	| 7 |	10.4
Staphylococcus aureus | 6 |	15.5

### Computing Environment 

#### R 
- R v3.6.1
- Bioconductor v3.10
- DADA2 v1.14.0 
- ShortRead v1.44.0
- ggplot2 v3.2.1

#### CentOS6 Linux
- Bowtie v2.2.9
- samtools v1.3.1
- Alfred v0.1.17
- Kraken v2.0.8-beta
- Bracken v2.5
- trimmomatic v0.36
- FLASH v1.2.11
- cutadapt v1.16


### LoopSeq 16S Contig QC & Filtering

Load fastq file containing raw contigs as downloaded from the Loop Genomics cloud platform.

[Raw contigs fastq file](https://github.com/nico-chung/loopseq/raw/master/data/zymo.rar)

```
#R
library(dada2)
library(ShortRead)

F_prim <- "AGAGTTTGATCMTGGCTCAG" #forward primer
R_prim <- "TACCTTGTTACGACTT" #reverse primer

zymo <- file.path("zymo.fastq")

#Basic statistics of raw contigs
summary(nchar(getSequences(zymo)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    500    1494    1508    1430    1510    1572 

#Distribution of raw contigs
hist(nchar(getSequences(zymo)), 100) 
```
![](https://github.com/nico-chung/loopseq/blob/master/figures/zymo_rawcontigs.png)

```
#Remove forward and reverse primers. This step also ensures that the contigs that remain after filtering are the full length 16S gene. 
removePrimers(zymo, paste0("zymo_noprimers.fastq"), primer.fwd=F_prim, primer.rev=rc(R_prim), compress=FALSE, max.mismatch=2, allow.indels=FALSE, verbose=TRUE)
Read in 38629, output 27123 (70.2%) filtered sequences.

#Basic statistics of filtered contigs (full length 16S, primers trimmed) 
summary(nchar(getSequences(zymo_noprimers)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1426    1464    1472    1474    1474    1509 

#Distribution of filtered contigs (full length 16S, primers trimmed)
zymo_noprimers <- file.path("zymo_noprimers.fastq")
hist(nchar(getSequences(zymo_noprimers)), 10) 
```
[Filtered contigs fastq file](https://github.com/nico-chung/loopseq/raw/master/data/zymo_noprimers.rar)

![](https://github.com/nico-chung/loopseq/blob/master/figures/zymo_noprimers.png) 


## Mapping to Provided References

Zymo replaced five strains in the Zymo standard with similar strains beginning with Lot ZRC190633. [The sample we sequenced was from LOT NUMBER ###]. As has been noted elsewhere (Callahan Benjamin C., et al. "High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide resolution", *Nucleic Acids Research*, 47.18 (2019): e103), the exact strains used are not disclosed and the [reference sequences](https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip) provided by Zymo Research do not exactly correspond to the strains provided. There are two versions of the references: version 1 and version 2. In both versions, there are systematic nucleotide substitutions that cannot be explained by sequencing or contig assembly error. The references include whole genomes and 16S genes only. In version 1 of the references, the whole genome references are in multiple contigs and in the 16S gene references, only one 16S copy is given for each species. Version 2 of references aren’t correct either but is still more correct than version 1. The new version has complete whole genome references and all 16S gene copies are given for each species. 

Concatenate all bacterial whole genome fasta files together and build a Bowtie2 index with version 2 of the whole genome references and map the filtered Zymo contigs. 

```
cat *.fasta > zymo_bacteria_genome_ref.fasta

bowtie2-build zymo_bacteria_genome_ref.fasta zymo_bacteria_genome_ref

bowtie2 --very-sensitive -p 40 --un zymo.un.fastq -x zymo_bacteria_genome_ref/zymo_bacteria_genome_ref -q zymo_noprimers.fastq -S zymo.sam
27123 reads; of these:
  27123 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    27123 (100.00%) aligned >1 times
100.00% overall alignment rate 
```

Many of the 16S genes are highly similar/identical and multiple valid alignments may occur. Even though the provided reference sequences are incorrect, given that the contigs are full-length 16S sequences and the contigs are accurate enough, all eight bacterial species are sufficiently phylogenetically distant that the multiple mappings occur within a single genome reference. It is the case generally that accurate full-length 16S sequences are able to discriminate between species (Johnson, Jethro S., et al. "Evaluation of 16S rRNA gene sequencing for species and strain-level microbiome analysis." *Nature communications* 10.1 (2019): 1-11.). There are exceptions (eg. Mycobacterium genus) where inter-species 16S copies are so similar (~99.9%), that with even minor sequence errors, it may be difficult to unambiguously map the sequences to the correct references.

```
#Extract names of hits and count number of each unique hit
samtools view -F 4 zymo.sam | cut -f3 | sort | uniq -c

   5692 BS.pilon.polished.v3.ST170922
   2627 Enterococcus_faecalis_complete_genome
   2962 Escherichia_coli_chromosome
   3066 Lactobacillus_fermentum_complete_genome
   4722 Listeria_monocytogenes_complete_genome
   1151 Pseudomonas_aeruginosa_complete_genome
   3032 Salmonella_enterica_complete_genome
   3871 Staphylococcus_aureus_chromosome

```

| **Species** | **Expected contig %** | **Mapped contig %** |
|-------------|-----------------------|---------------------|
Bacillus subtilis | 17.4 |	21.0
Enterococcus faecalis | 9.9 | 9.7
Escherichia coli | 10.1 | 10.9
Lactobacillus fermentum | 18.4 | 11.3
Listeria monocytogenes | 14.1 | 17.4
Pseudomonas aeruginosa | 4.2 | 4.2
Salmonella enterica | 10.4 | 11.2
Staphylococcus aureus | 15.5 | 14.3

Species ratios are as expected. 

```
#R
library(ggplot2)
zymo <- read.table("zymo_exp_vs_map.txt", header=TRUE, sep="\t")
ggplot(zymo, aes(fill=Species, y=Percentage, x=Contigs)) + geom_bar(position="stack", stat="identity")
```
![](https://github.com/nico-chung/loopseq/blob/master/figures/zymo_exp_vs_map.png) 


### Intergenomic 16S Gene Copies & Contig Accuracy

*B. subtilis* has the most 16S gene copies (ten) of the species in the Zymo standard so it can provide a fine-grain view of relative intergenomic 16S gene ratios. It also appears to be one of the correct references for the strain provided (no systematic mismatches between contig and reference) so it would allow for an unbiased assessment of LoopSeq contig accuracy.   

Extract the full length 16S contigs that orginally mapped to *B. subtilis*.
```
samtools view zymo.sam | cut -f1,3 | awk -F "\t" '{ if($2 == "BS.pilon.polished.v3.ST170922") { print $1}}' > bsubtilis_contigs_list
seqtk subseq zymo_noprimers.fastq bsubtilis_contigs_list > bsubtilis_contigs.fastq
```
Sequences for each of the ten 16S genes is given in the reference fasta file. [Multiple sequence alignment](https://github.com/nico-chung/loopseq/blob/master/data/zymo_bsubtilis_clustal_msa.txt) shows that there are five non-redundant copies . Copy #1, #4, #8, #9, #10 are identical. #2, #3 are identical. #5, #6, #7 are all unique. To prevent multi-mapping, build a bowtie2 index only using #1, #2, #5, #6, #7. Thus a mapping ratio of 5:2:1:1:1 is expected.  
```
bowtie2 --very-sensitive -p 40 --un bsubtilis_contigs_trimmed.un.fastq -x zymo_bsubtilis_ref/Bacillus_subtilis_16S_nr -q bsubtilis_contigs.fastq -S bsubtilis.sam
5692 reads; of these:
  5692 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    5692 (100.00%) aligned >1 times
100.00% overall alignment rate #all mapped to 16S copies

#Extract reference hits
samtools view -F 4 bsubtilis.sam | cut -f3 | sort | uniq -c > bsubtilis_hits_counts

   2466 Bacillus_subtilis_16S_1
   1482 Bacillus_subtilis_16S_2
    604 Bacillus_subtilis_16S_5
    617 Bacillus_subtilis_16S_6
    523 Bacillus_subtilis_16S_7
```
| **16S Gene Copy** | **Expected contig %** | **Mapped contig %** |
|-------------------|-----------------------|---------------------|
Bacillus_subtilis_16S_1 |	50 | 43.3
Bacillus_subtilis_16S_2	| 20 | 26
Bacillus_subtilis_16S_5	| 10 | 10.6
Bacillus_subtilis_16S_6	| 10 | 10.8
Bacillus_subtilis_16S_7	| 10 | 9.1

Intergenomic 16S gene ratios are as expected.

```
zymo_bsubtilis <- read.table("zymo_bsubtilis_16S_copies.txt", header=TRUE, sep="\t")
ggplot(zymo_bsubtilis, aes(fill=Contigs, y=Percentage, x=Gene.Copy)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```
![](https://github.com/nico-chung/loopseq/blob/master/figures/zymo_bsubtilis_16S_copies.png)

Visually inspect the alignments to check if there is anything unusual about the alignments. Sort and index the SAM file and load in IGV.
```
samtools view -b -S bsubtilis.sam | samtools sort -o bsubtilis.sorted.bam - && samtools index bsubtilis.sorted.bam
```
**Bacillus_subtilis_16S_1**
![](https://github.com/nico-chung/loopseq/blob/master/figures/zymo_bsubtilis_16S_1_igv.png)

**Bacillus_subtilis_16S_2**
![](https://github.com/nico-chung/loopseq/blob/master/figures/zymo_bsubtilis_16S_2_igv.png)

**Bacillus_subtilis_16S_5**
![](https://github.com/nico-chung/loopseq/blob/master/figures/zymo_bsubtilis_16S_5_igv.png)

**Bacillus_subtilis_16S_6**
![](https://github.com/nico-chung/loopseq/blob/master/figures/zymo_bsubtilis_16S_6_igv.png)

**Bacillus_subtilis_16S_7**
![](https://github.com/nico-chung/loopseq/blob/master/figures/zymo_bsubtilis_16S_7_igv.png)

Quantify the quality of the alignments in Alfred. 

```
alfred qc -r zymo_bsubtilis_ref/Bacillus_subtilis_16S_nr.fasta -o alignment_qc.tsv.gz bsubtilis.sorted.bam
zgrep ^ME alignment_qc.tsv.gz | cut -f 2- | datamash transpose | column -t > zymo_bsubtilis_16S_alignment_qc.tsv
```
[Result file](https://github.com/nico-chung/loopseq/blob/master/data/zymo_bsubtilis_16S_alignment_qc.tsv)

##### QC Summary
- MatchRate                       0.999798
- MismatchRate                    0.000201824
- DeletionRate                    2.38704e-06
- InsertionRate                   2.14834e-06

The error rates are extremely low considering these contigs have only been filtered to cover the full-length of the 16S gene. They have not been quality filtered, denoised or any other type of algorithmic error correction. Indels are negligible. Since consensus sequence accuracy is contigent on read depth, higher per base accuracy and a higher proportion of full-length 16S contigs can be expected with increased sequencing yield. But since there is a limited number of UMI barcodes per 16S fragment, there will be an upper limit to gains in accuracy. 

Examine the distribution of mismatch errors in the contigs. Use a custom script to quantify the position and frequencies of variants.
```
samtools mpileup -B -L 1 -Q 0 -F 0 -f zymo_bsubtilis_ref/Bacillus_subtilis_16S_nr.fasta bsubtilis.sorted.bam > bsubtilis.pileup
pileup2SNVbaseCounts.pl -i bsubtilis.pileup -snvDepth 0 -indelDepth 0

```
[Result file](https://github.com/nico-chung/loopseq/blob/master/data/zymo_bsubtilis_mpileup_baseCount.csv)

The five non-redundant copies of the 16S gene are the same length, so they can be overlapped on the same plot. Only base mismatches are shown since there are practically no indels to visualize. 
```
#R
bsubtilis_variants <- read.csv("zymo_bsubtilis_mpileup_baseCount.csv",header=TRUE)
ggplot(data=bsubtilis_variants,aes(x=bsubtilis_variants$Position,fill=Reference,y=bsubtilis_variants$Variant_Ratio)) + geom_bar(stat="identity")
```
![](https://github.com/nico-chung/loopseq/blob/master/figures/zymo_bsubtilis_snp_var.png)

Mismatch errors are randomly distributed. The mismatches at the two peaks (gene copy #2 at pos: 285 [G->A] and gene copy #1 pos: 268 [G->A]) likely represent real variants since they are too high to be explained by sequencing or assembly error: the sequencing library consists of short randomly fragmented segments of the 16S gene. Consequently, the true mismatch error rate may be even lower: the bases in the two peaks combined accounts for 454/1691 total mismatch errors.


### Comparison with Illumina Short Paired-End Reads

Short read 16S barcoding studies frequently sequence the V3V4 hypervariable region in 300bp paired-end mode. These paired-end reads are then filtered, trimmed and merged into ~440bp fragments which are then used for denoising and/or OTU clustering and classification. To make a fair comparison with LoopSeq contigs, the short reads are processed in a way as might be seen in a typical 16S barcoding study. The denoising step is omitted because it would potentially produce error-free reads; as would be the case if LoopSeq contigs were to undergo denoising. This would not permit a meaningful comparison of accuracies between the two methodologies.

The fastq files are too large to host on Github so have not been provided.

```
#Check how many pairs of reads in raw fastq file
echo $(zcat zymo-v3v4_R1.fastq.gz|wc -l)/4|bc
3278936

#Map against the PhiX reference to remove spike-in sequences
bowtie2 -p 24 --un-conc zymo-v3v4.fastq -x phiX/phiX -1 zymo-v3v4_R1.fastq.gz -2 zymo-v3v4_R2.fastq.gz -S zymo-v3v4_phiX.sam

#Filter out reads with low average quality and clip low quality bases in the tails
java -jar trimmomatic-0.36.jar PE zymo-v3v4.1.fastq.gz zymo-v3v4.2.fastq.gz zymo-3V4_clean_1.fastq.gz zymo-v3v4_unclean_1.fastq.gz zymo-v3v4_clean_2.fastq.gz zymo-v3v4_unclean_2.fastq.gz AVGQUAL:20
java -jar trimmomatic-0.36.jar PE zymo-v3v4_clean_1.fastq.gz zymo-v3v4_clean_2.fastq.gz zymo-v3v4_good_1.fastq.gz zymo-v3v4_bad_1.fastq.gz zymo-v3v4_good_2.fastq.gz zymo-v3v4_bad_2.fastq.gz TRAILING:20 

#Merge reads
flash -t 24 -m 20 -M 250 -d merged -o zymo-v3v4 zymo-v3v4_good_1.fastq.gz zymo-v3v4_good_2.fastq.gz 

#Trim V3V4 primers. 3’ primer sequence needs to be given in the reverse complement.
cutadapt -e 0.05 -g CCTACGGGNGGCWGCAG -a GGATTAGATACCCBDGTAGT -o zymo-v3v4.noprimers.fastq zymo-v3v4.extendedFrags.fastq --discard-untrimmed --no-indels
```

Check merged read distribution.
```
#R
zymo_v3V4_noprimers <- file.path("zymo-v3v4.noprimers.fastq")
summary(nchar(getSequences(zymo_v3V4_noprimers)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   17.0   444.0   444.0   444.1   444.0   574.0 
   
#Some merged reads are too short/long to be correct so these need to be filtered out
fastqFilter(zymo_v3V4_noprimers, "zymo-v3v4_noprimers_length.fastq", minLen=440, maxLen=450, verbose=TRUE, compress=FALSE)
Read in 2761523, output 2239559 (81.1%) filtered sequences.

```

```
#Map to Zymo whole genome bacteria references
bowtie2 --very-sensitive -p 40 --un zymo-v3v4.un.fastq -x zymo_bacteria_genome_ref/zymo_bacteria_genome_ref -q zymo-v3v4_noprimers_length.fastq -S zymo-v3v4.sam
2239559 reads; of these:
  2239559 (100.00%) were unpaired; of these:
    2051 (0.09%) aligned 0 times
    51 (0.00%) aligned exactly 1 time
    2237457 (99.91%) aligned >1 times
99.91% overall alignment rate

#Extract names of hits and count number of each unique hit
samtools view -F 4 zymo-v3v4.sam | cut -f3 | sort | uniq -c
 348730 BS.pilon.polished.v3.ST170922
 202215 Enterococcus_faecalis_complete_genome
 266642 Escherichia_coli_chromosome
 351894 Lactobacillus_fermentum_complete_genome
 282619 Listeria_monocytogenes_complete_genome
 177842 Pseudomonas_aeruginosa_complete_genome
 268401 Salmonella_enterica_complete_genome
 339165 Staphylococcus_aureus_chromosome

```

| **Species** | **Expected %** | **Mapped %** |
|-------------|----------------|------------------------------|
Bacillus subtilis | 17.4 |	15.6  
Enterococcus faecalis | 9.9 | 9.0 
Escherichia coli | 10.1 | 11.9 
Lactobacillus fermentum | 18.4 | 15.7 
Listeria monocytogenes | 14.1 | 12.6  
Pseudomonas aeruginosa | 4.2 | 8.0 
Salmonella enterica | 10.4 | 12.0 
Staphylococcus aureus | 15.5 | 15.2

The species ratios are as expected. It does not appear that short read methodologies significanlty skew the relative abundances of the species examined.

To compare read accuracies with the LoopSeq contigs, extract the merged reads that mapped to *B. subtilis* and map the reads to the five non-redundant *B. subtilis* 16S gene copy references. 
```
samtools view zymo-v3v4.sam | cut -f1,3 | awk -F "\t" '{ if($2 == "BS.pilon.polished.v3.ST170922") { print $1}}' > bsubtilis-v3v4_reads_list
seqtk subseq zymo-v3v4.noprimers.fastq bsubtilis-v3v4_reads_list > bsubtilis-v3v4_reads.fastq
bowtie2 --very-sensitive -p 40 --un bsubtilis-v3v4_reads.un.fastq -x zymo_bsubtilis_ref/Bacillus_subtilis_16S_nr -q bsubtilis-v3v4_reads.fastq -S bsubtilis-v3v4_reads.sam
348730 reads; of these:
  348730 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    4 (0.00%) aligned exactly 1 time
    348726 (100.00%) aligned >1 times
100.00% overall alignment rate

#Extract names of hits and count number of each unique hit
samtools view -F 4 bsubtilis-v3v4_reads.sam | cut -f3 | sort | uniq -c
  78053 Bacillus_subtilis_16S_1
  77557 Bacillus_subtilis_16S_2
  38240 Bacillus_subtilis_16S_5
  77371 Bacillus_subtilis_16S_6
  77509 Bacillus_subtilis_16S_7

```
| **16S Gene Copy** | **Expected contig %** | **Mapped reads %**  |
|-------------------|-----------------------|---------------------|
Bacillus_subtilis_16S_1 | 50 | 22.4
Bacillus_subtilis_16S_2	| 20 | 22.3
Bacillus_subtilis_16S_5	| 10 | 11.0
Bacillus_subtilis_16S_6	| 10 | 22.2
Bacillus_subtilis_16S_7	| 10 | 22.2

Many of the nucleotide differences between the five non-redundant *B. subtilis* 16S copies are not within the V3V4 region (pos: 348-792). Within the V3V4 region of *B. subtilis*, only 16S copy #5 is different (pos: 373 C->A) while copies #1, #2, #6 and #7 are all identical in the V3V4 region. Reads corresponding to these identical sequences would map equally well to any of the 16S copies so are randomly placed, resulting in the approximately equal ratios observed for copy #1, #2, #6 and #7. This shows that sequencing short hypervariable regions may not permit an accurate assessment of intergenomic 16S gene ratios. This limits the use of short reads in inferring bacterial subtypes which may be distinguished on the basis of their intergenomic 16S copies. In general, the use of short hypervariable regions (V1-V2, V1-V3, V4, etc) has limited utility in identifying prokaryotes to the species level (Johnson, Jethro S., et al. "Evaluation of 16S rRNA gene sequencing for species and strain-level microbiome analysis." *Nature Communications* 10.1 (2019): 1-11.).

Quantify the quality of the alignments in Alfred. 

```
samtools view -b -S bsubtilis-v3v4_reads.sam | samtools sort -o bsubtilis-v3v4_reads.sorted.bam - && samtools index bsubtilis-v3v4_reads.sorted.bam
alfred qc -r zymo_bsubtilis_ref/Bacillus_subtilis_16S_nr.fasta -o bsubtilis-v3v4_alignment_qc.tsv.gz bsubtilis-v3v4_reads.sorted.bam
zgrep ^ME bsubtilis-v3v4_alignment_qc.tsv.gz | cut -f 2- | datamash transpose | column -t > zymo_bsubtilis-v3v4_alignment_qc.tsv
```
[Result file](https://github.com/nico-chung/loopseq/blob/master/data/zymo_bsubtilis-v3v4_alignment_qc.tsv)

##### QC Summary
- MatchRate                       0.992513
- MismatchRate                    0.00748692
- DeletionRate                    8.36867e-05
- InsertionRate                   3.11562e-05       

The mismatch, insertion and deletion error rates for short reads are at least an order of magnitude larger than the LoopSeq contigs. Optimized experimental protocols and stricter quality filtering criteria may result in reads with even lower mean error rates. There are also denoising algorithms that can correct PCR and sequencing errors, but such corrected reads would not accurately reflect the underlying error profile. Conversely, not performaing any quality filtering would result in reads with a higher error rate, thus not meeting basic quality standards for analysis. The filtering utilized above simply represents a typical pre-processing workflow which allows us to make a fair comparison with LoopSeq contigs which have only been filtered to full-length.

## Mapping to NCBI RefSeq References

The eight bacterial species in the Zymo mock community are sufficiently phylogenetically distant that if given reads/contigs are long and accurate enough, reads/contigs are unlikely to me errorneously mapped to the wrong species. However, it is not always the case that species within a sample are known and/or reference sequences for the species are available. It it more often the case that samples with unknown composition are compared to a more general database. To determine if full-length 16S sequences can be used to accurately identify the correct species using a general database, the full-length LoopSeq 16S Zymo contigs were mapped to all bacterial genomes (35,494 genomes as of October 2019) in the NCBI RefSeq CG (complete genomes) database. The database contains genomes of the same species as those in the Zymo sample (although it is unclear if the strains are the same since this information is lacking). There are also many genomes of closely related species (species from the same genus) that could potentially counfound a direct mapping approach to species identification. 

Bacterial genomes were download from NCBI RefSeq using Kraken2. A Bowtie2 index was built (this takes some time) and full-length LoopSeq 16S contigs were mapped to the database.
```
bowtie2 --very-sensitive -p 40 --un zymo_refseq.un.fastq -x bt2_refseq_bacteria_genomes/refseq_bacteria -q zymo_noprimers.fastq -S zymo_refseq.sam
27123 reads; of these:
  27123 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    1 (0.00%) aligned exactly 1 time
    27122 (100.00%) aligned >1 times
100.00% overall alignment rate

#Extract reference hits
samtools view -F 4 zymo_refseq.sam | cut -f3 > zymo_refseq_hits
```
The bacterial references and taxonomy files were downloaded with Kraken2. The header does not contain the original species names. It does have taxids which can be converted and reformatted using taxonkit. 
```
#Cut out the taxids from the header and convert to a 7-rank taxonomy
cat zymo_refseq_hits | cut -f2 -d"|" \
| taxonkit lineage --data-dir Kraken2/standard_DB/taxonomy \
| taxonkit reformat --data-dir Kraken2/standard_DB/taxonomy \
| cut -f3 | sort | uniq -c
      2 Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus halotolerans
      1 Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus sp. JS
      1 Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus sp. SJ-10
   5689 Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus subtilis
   4722 Bacteria;Firmicutes;Bacilli;Bacillales;Listeriaceae;Listeria;Listeria monocytogenes
     81 Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus argenteus
   3789 Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus
   2627 Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecalis
   3066 Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus fermentum
   2952 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli
      1 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Klebsiella;Klebsiella aerogenes
   3032 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella enterica
      7 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Shigella;Shigella dysenteriae
      2 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Shigella;Shigella flexneri
   1142 Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas aeruginosa
      7 Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas fluorescens
      2 Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas sp. AK6U
```

| **Species** | **Expected %** | **Mapped %** |
|-------------|----------------|------------------------------|
Bacillus subtilis | 17.4 | 21.0
Enterococcus faecalis | 9.9 | 9.7	
Escherichia coli | 10.1 | 10.9	
Lactobacillus fermentum | 18.4 | 11.3
Listeria monocytogenes | 14.1 | 17.4
Pseudomonas aeruginosa | 4.2 | 4.2
Salmonella enterica | 10.4 | 11.2
Staphylococcus aureus | 15.5 | 14.0	
Other species | 0 | 0.4

```
#R
library(ggplot2)
zymo_refseq <- read.table("zymo_refseq_exp_vs_map.txt", header=TRUE, sep="\t")
ggplot(zymo_refseq, aes(fill=Species, y=Percentage, x=Contigs, group=customOrder)) + geom_bar(position="stack", stat="identity")
```
![](https://github.com/nico-chung/loopseq/blob/master/figures/zymo_refseq_exp_vs_map.png) 

Only 0.4% of contigs were mapped to incorrect species. With the exception of a single contig that mapped to a different genus (Klebsiella), the incorrect mappings were consistent with the genuses in the Zymo sample. This demonstrates that it is possible to accurately classify bacteria to the species level using full-length LoopSeq 16S contigs against a general whole genome database. This also demonstrates that even minor differences may produce incorrect species mappings.


## Classification with LCA Classifier

Although the direct mapping approach may provide higher classification precision when references are available, it is frequently the case that exact references are not availble. One may nontheless still be interested in classifying the sequence to higher taxonomic levels (genus, family etc.). For such cases, lowest common ancestor (LCA) classifiers may be utilized. Commonly used 16S classification tools include mothur and QIIME used in conjunction with the SILVA, GreenGenes or RDP databases. The databases may contain 16S sequences from many representative species but the full intergenomic 16S diversity is usually not represented. Thus, it may be more advantageous to use whole genome references (which contain all the 16S copies) when they are available. The obvious disadvantage is that fewer species are represented as compared to 16S databases.    

Kraken2 is a LCA classifer that searches for exact alignment of k-mers against an index built from whole genome references. The database used is the standard Kraken2 database (NCBI RefSeq CG: bacteria, archaea, virus, human, common cloning vectors).


```
kraken2 --db Kraken2/standard_DB --threads 24 --out zymo.kraken --report zymo.kreport zymo_noprimers.fastq
```

The result columns are interpreted as follows: 
1. Percentage of fragments covered by the clade rooted at this taxon
2. Number of fragments covered by the clade rooted at this taxon
3. Number of fragments assigned directly to this taxon
4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
   (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
   Taxa that are not at any of these 10 ranks have a rank code that is
   formed by using the rank code of the closest ancestor rank with
   a number indicating the distance from that rank.  E.g., "G2" is a
   rank code indicating a taxon is between genus and species and the
   grandparent taxon is at the genus rank.
5. NCBI taxonomic ID number
6. Indented scientific name

```
100.00	27123	0	R	1	root
100.00	27123	0	R1	131567	  cellular organisms
100.00	27123	0	D	2	    Bacteria
 73.66	19978	0	D1	1783272	      Terrabacteria group
 73.66	19978	0	P	1239	        Firmicutes
 73.66	19978	0	C	91061	          Bacilli
 52.66	14284	0	O	1385	            Bacillales
 20.98	5691	0	F	186817	              Bacillaceae
 20.98	5691	21	G	1386	                Bacillus
 20.86	5658	3877	G1	653685	                  Bacillus subtilis group
  6.51	1767	1759	S	1423	                    Bacillus subtilis
  0.03	7	4	S1	135461	                      Bacillus subtilis subsp. subtilis
  0.00	1	1	S2	224308	                        Bacillus subtilis subsp. subtilis str. 168
  0.00	1	1	S2	1052588	                        Bacillus subtilis subsp. subtilis str. RO-NN-1
  0.00	1	1	S2	1302650	                        Bacillus subtilis subsp. subtilis str. BAB-1
  0.00	1	0	S1	96241	                      Bacillus subtilis subsp. spizizenii
  0.00	1	1	S2	1052585	                        Bacillus subtilis subsp. spizizenii TU-B-10
  0.04	11	0	G2	1938374	                    Bacillus amyloliquefaciens group
  0.03	8	8	S	492670	                      Bacillus velezensis
  0.01	3	2	S	1390	                      Bacillus amyloliquefaciens
  0.00	1	1	S1	1034836	                        Bacillus amyloliquefaciens XH7
  0.01	2	2	S	1402	                    Bacillus licheniformis
  0.00	1	1	S	119858	                    Bacillus sonorensis
  0.02	5	1	G1	86661	                  Bacillus cereus group
  0.01	3	1	S	1396	                    Bacillus cereus
  0.00	1	1	S1	1454382	                      Bacillus cereus D17
  0.00	1	0	S1	1179100	                      Bacillus cereus biovar anthracis
  0.00	1	1	S2	637380	                        Bacillus cereus biovar anthracis str. CI
  0.00	1	0	S	1392	                    Bacillus anthracis
  0.00	1	1	S1	673518	                      Bacillus anthracis str. A16R
  0.01	3	3	S	152268	                  Bacillus litoralis
  0.01	2	0	S	1404	                  Bacillus megaterium
  0.01	2	2	S1	592022	                    Bacillus megaterium DSM 319
  0.01	2	0	S	665099	                  Bacillus oceanisediminis
  0.01	2	2	S1	1196031	                    Bacillus oceanisediminis 2691
 17.41	4722	0	F	186820	              Listeriaceae
 17.41	4722	0	G	1637	                Listeria
17.41	4722	4722	S	1639	                  Listeria monocytogenes 
 14.27	3871	0	F	90964	              Staphylococcaceae
 14.27	3871	2725	G	1279	                Staphylococcus
  4.14	1122	1115	S	1280	                  Staphylococcus aureus
  0.03	7	1	S1	46170	                    Staphylococcus aureus subsp. aureus
  0.02	5	5	S2	985006	                      Staphylococcus aureus subsp. aureus LGA251
  0.00	1	1	S2	548473	                      Staphylococcus aureus subsp. aureus TCH60
  0.01	4	4	S	1295	                  Staphylococcus schleiferi
  0.01	3	3	S	29379	                  Staphylococcus auricularis
  0.01	3	3	S	28035	                  Staphylococcus lugdunensis
  0.01	2	2	S	1282	                  Staphylococcus epidermidis
  0.01	2	2	S	29384	                  Staphylococcus kloosii
  0.00	1	1	S	70258	                  Staphylococcus piscifermentans
  0.00	1	1	S	70255	                  Staphylococcus condimenti
  0.00	1	1	S	246432	                  Staphylococcus equorum
  0.00	1	0	S	29388	                  Staphylococcus capitis
  0.00	1	1	S1	72758	                    Staphylococcus capitis subsp. capitis
  0.00	1	1	S	29385	                  Staphylococcus saprophyticus
  0.00	1	1	S	283734	                  Staphylococcus pseudintermedius
  0.00	1	1	S	308354	                  Staphylococcus simiae
  0.00	1	1	S	1286	                  Staphylococcus simulans
  0.00	1	1	S	1283	                  Staphylococcus haemolyticus
  0.00	1	1	S	2044912	                  Staphylococcus sp. SDB 2975
 20.99	5694	0	O	186826	            Lactobacillales
 11.30	3066	0	F	33958	              Lactobacillaceae
 11.30	3066	0	G	1578	                Lactobacillus
 11.30	3066	3064	S	1613	                  Lactobacillus fermentum
  0.01	2	2	S1	767453	                    Lactobacillus fermentum F-6
  9.69	2628	0	F	81852	              Enterococcaceae
  9.69	2628	0	G	1350	                Enterococcus
  9.69	2628	2628	S	1351	                  Enterococcus faecalis
 26.34	7145	0	P	1224	      Proteobacteria
 26.34	7145	0	C	1236	        Gammaproteobacteria
 22.10	5994	1	O	91347	          Enterobacterales
 22.10	5993	2038	F	543	            Enterobacteriaceae
 11.19	3034	0	G	590	              Salmonella
 11.19	3034	3013	S	28901	                Salmonella enterica
  0.08	21	9	S1	59201	                  Salmonella enterica subsp. enterica
  0.01	3	1	S2	149539	                    Salmonella enterica subsp. enterica serovar Enteritidis
  0.01	2	2	S3	1244111	                      Salmonella enterica subsp. enterica serovar Enteritidis str. EC20110357
  0.01	2	0	S2	57046	                    Salmonella enterica subsp. enterica serovar Paratyphi C
  0.01	2	2	S3	476213	                      Salmonella enterica subsp. enterica serovar Paratyphi C str. RKS4594
  0.00	1	1	S2	90370	                    Salmonella enterica subsp. enterica serovar Typhi
  0.00	1	0	S2	90371	                    Salmonella enterica subsp. enterica serovar Typhimurium
  0.00	1	1	S3	1008297	                      Salmonella enterica subsp. enterica serovar Typhimurium str. 798
  0.00	1	1	S2	98360	                    Salmonella enterica subsp. enterica serovar Dublin
  0.00	1	0	S2	108619	                    Salmonella enterica subsp. enterica serovar Newport
  0.00	1	1	S3	1454618	                      Salmonella enterica subsp. enterica serovar Newport str. USDA-ARS-USMARC-1925
  0.00	1	1	S2	605	                    Salmonella enterica subsp. enterica serovar Pullorum
  0.00	1	0	S2	58712	                    Salmonella enterica subsp. enterica serovar Anatum
  0.00	1	1	S3	1454592	                      Salmonella enterica subsp. enterica serovar Anatum str. CDC 06-0532
 0.00	1	1	S2	29474	                    Salmonella enterica subsp. enterica serovar California
  3.33	903	0	G	561	              Escherichia
  3.33	902	899	S	562	                Escherichia coli
  0.01	2	0	S1	83333	                  Escherichia coli K-12
  0.01	2	2	S2	1110693	                    Escherichia coli str. K-12 substr. MDS42
  0.00	1	1	S1	2048777	                  Escherichia coli O15:H11
  0.00	1	1	S	1499973	                Escherichia marmotae
  0.04	10	3	G	620	              Shigella
  0.03	7	3	S	623	                Shigella flexneri
  0.01	3	0	S1	424718	                  Shigella flexneri 5a
  0.01	3	3	S2	1086030	                    Shigella flexneri 5a str. M90T
  0.00	1	1	S1	42897	                  Shigella flexneri 2a
  0.01	4	0	G	570	              Klebsiella
  0.01	4	4	S	573	                Klebsiella pneumoniae
  0.01	3	0	G	547	              Enterobacter
  0.01	2	0	G1	354276	                Enterobacter cloacae complex
  0.00	1	1	S	550	                  Enterobacter cloacae
  0.00	1	1	S	1915310	                  Enterobacter cloacae complex sp. ECNIH7
  0.00	1	1	S	881260	                Enterobacter bugandensis
  0.00	1	0	F1	191675	              unclassified Enterobacteriaceae
  0.00	1	0	F2	84563	                ant, tsetse, mealybug, aphid, etc. endosymbionts
  0.00	1	0	F3	146507	                  aphid secondary symbionts
  0.00	1	1	S	1199245	                    secondary endosymbiont of Ctenarytaina eucalypti
  4.24	1151	0	O	72274	          Pseudomonadales
  4.24	1151	0	F	135621	            Pseudomonadaceae
  4.24	1151	1143	G	286	              Pseudomonas
  0.01	2	0	S	101564	                Pseudomonas alcaliphila
  0.01	2	2	S1	741155	                  Pseudomonas alcaliphila JAB1
  0.01	2	0	G1	136842	                Pseudomonas chlororaphis group
  0.01	2	2	S	296	                  Pseudomonas fragi
  0.00	1	1	S	1981174	                Pseudomonas sp. M30-35
  0.00	1	1	S	1302376	                Candidatus Pseudomonas adelgestsugas
  0.00	1	0	G1	136841	                Pseudomonas aeruginosa group
  0.00	1	1	S	287	                  Pseudomonas aeruginosa
  0.00	1	0	G1	136843	                Pseudomonas fluorescens group
  0.00	1	1	S	47883	                  Pseudomonas synxantha
```

Most contigs were correctly assigned to the species level. The exception was *P. aureginosa*, for which most contigs could only be assigned to the genus level. However, the percentage of contigs (4.24%) that were classified to the Pseudomonas genus does conform to expectations. This lack of species resolution for *P. aureginosa* is due to characteristics of the Pseudomonas genus in general and not a particular shortcoming of LoopSeq contigs. A previous study (Gomila, Margarita, et al. "Phylogenomics and systematics in Pseudomonas." *Frontiers in Microbiology* 6 (2015): 214.) found that: 

> *In conclusion, because the resolution of the 16S rRNA tree was not sufficient to differentiate 63 genomes from other closely    related Pseudomonas species, the classification of these bacteria should follow the phylogeny of the housekeeping genes until the whole genome sequence of the type strains of all Pseudomonas species is known.*

Using Bracken, the companion tool to Kraken that redistributes all the classified reads to a single taxonomic level, it is apparent that at the genus level, the ratios of all eight species are as expected.

```
bracken -d Kraken2/standard_DB -i zymo.kreport -o zymo_genus.bracken -r 1500 -l G -t 1
```
[Result file](https://github.com/nico-chung/loopseq/blob/master/data/zymo_genus.bracken)

| **Genus** | **Expected %** | **Classified %** |
|-----------|----------------|-----------------|
Bacillus | 17.4 | 21.0
Enterococcus | 9.9 | 9.7
Escherichia | 10.1 | 10.0
Lactobacillus | 18.4 | 11.3
Listeria | 14.1 | 17.4
Pseudomonas | 4.2 | 4.2
Salmonella | 10.4 | 11.9
Staphylococcus | 15.5 | 14.3
Other species | 0 | 0.3


Redistributing contigs to the species level results in a higher proportion of false positives. Of the expected 4.2% *P. aureginosa* species, approximately 3/4 were classied to other Pseudomonas species. 
```
bracken -d Kraken2/standard_DB -i zymo.kreport -o zymo_species.bracken -r 1500 -l S -t 1
```
[Result file](https://github.com/nico-chung/loopseq/blob/master/data/zymo_species.bracken)

| **Species** | **Expected %** | **Classified %** |
|-----------|----------------|-----------------|
Bacillus subtilis	| 17.4 | 19.0
Enterococcus faecalis | 9.9 | 9.7
Escherichia coli	| 10.1 | 10.0
Lactobacillus fermentum |	18.4 | 11.3
Listeria monocytogenes | 14.1 | 17.4
Salmonella enterica | 10.4 | 11.9
Staphylococcus aureus | 15.5 | 13.8
**Pseudomonas aeruginosa** | 4.2 | 1.0
**Other Pseudomonas spp.** | 0 | 3.3 		
Other species	| 0 | 2.6
