# Madison Water Samples 

To demonstrate the the application of LoopSeq 16S kit on real world environmental samples, four water samples were collected from Madison (Wisconsin, USA): 
1) Home aquarium (Aquarium) 
2) Feynman Pond outside the Promega office buildings (FeynmanPond)
3) Residential rain barrel collecting rain runoff (RainBarrel) 
4) Riva Pond in a residential neighbourhood (RivaPond)

Samples were classified using an LCA classifer and diversity analysis was performed on the classification results.

## Experimental Protocol 

[Fill in experimental details here: Sample prep, kits used, platforms used, yield per sample, etc.]

The 16S + 18S kit was used but only the 16S contigs are examined. 

## Data Preparation

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

#### Madison-Aquarium

Load fastq file containing raw contigs as downloaded from the Loop Genomics cloud platform.

[Raw contigs fastq file](https://github.com/nico-chung/loopseq/raw/master/data/madison_aquarium.rar)

```
#R
library(dada2)
library(ShortRead)

F_prim <- "AGAGTTTGATCMTGGCTCAG" #forward primer
R_prim <- "TACCTTGTTACGACTT" #reverse primer

madison_aquarium <- file.path("madison_aquarium.fastq")

#Basic statistics of raw contigs
summary(nchar(getSequences(madison_aquarium)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    231    1393    1457    1355    1489    2921 
    
#Remove forward and reverse primers. This step also ensures that the contigs that remain after filtering are the full length 16S gene. 
removePrimers(madison_aquarium, "madison_aquarium_noprimers.fastq", primer.fwd=F_prim, primer.rev=rc(R_prim), compress=FALSE, max.mismatch=2, allow.indels=FALSE, verbose=TRUE)
6 sequences out of 18585 are being reverse-complemented.
Read in 18585, output 6967 (37.5%) filtered sequences.

#Basic statistics of filtered contigs (full length 16S, primers trimmed) 
madison_aquarium_noprimers <- file.path("madison_aquarium_noprimers.fastq")
summary(nchar(getSequences(madison_aquarium_noprimers)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1176    1442    1453    1443    1455    1544 
   
#Still some contigs that are too short. Filter for length.
filterAndTrim(madison_aquarium_noprimers, "madison_aquarium_noprimers_lengthfilter.fastq", minLen=1400, maxLen=1600, verbose=TRUE, compress=FALSE)
Read in 6967, output 6428 (92.3%) filtered sequences.

#Basic statistics of filtered contigs (full length 16S, primers trimmed, length filtered)
madison_aquarium_noprimers_lengthfilter <- file.path("madison_aquarium_noprimers_lengthfilter.fastq")
summary(nchar(getSequences(madison_aquarium_noprimers_lengthfilter)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1402    1451    1453    1447    1455    1544 

```
[Filtered contigs fastq file](https://github.com/nico-chung/loopseq/raw/master/data/madison_aquarium_noprimers_lengthfilter.rar)


#### FeynmanPond

Load fastq file containing raw contigs as downloaded from the Loop Genomics cloud platform.

[Raw contigs fastq file](https://github.com/nico-chung/loopseq/raw/master/data/madison_feynmanpond.rar)

```
#R
library(dada2)
library(ShortRead)

F_prim <- "AGAGTTTGATCMTGGCTCAG" #forward primer
R_prim <- "TACCTTGTTACGACTT" #reverse primer

madison_feynman_pond <- file.path("madison_feynman_pond.fastq")

#Basic statistics of raw contigs
summary(nchar(getSequences(madison_feynman_pond)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    222    1376    1451    1352    1487    2714 
    
#Remove forward and reverse primers. This step also ensures that the contigs that remain after filtering are the full length 16S gene. 
removePrimers(madison_feynman_pond, "madison_feynman_pond_noprimers.fastq", primer.fwd=F_prim, primer.rev=rc(R_prim), compress=FALSE, max.mismatch=2, allow.indels=FALSE, verbose=TRUE)
5 sequences out of 20189 are being reverse-complemented.
Read in 20189, output 8251 (40.9%) filtered sequences.

#Basic statistics of filtered contigs (full length 16S, primers trimmed) 
madison_feynman_pond_noprimers <- file.path("madison_feynman_pond_noprimers.fastq")
summary(nchar(getSequences(madison_feynman_pond_noprimers)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1376    1434    1446    1442    1457    1522 

```
[Filtered contigs fastq file](https://github.com/nico-chung/loopseq/raw/master/data/madison_feynman_pond_noprimers.rar)


#### RainBarrel

Load fastq file containing raw contigs as downloaded from the Loop Genomics cloud platform.

[Raw contigs fastq file](https://github.com/nico-chung/loopseq/raw/master/data/madison_rain_barrel.rar)

```
#R
library(dada2)
library(ShortRead)

F_prim <- "AGAGTTTGATCMTGGCTCAG" #forward primer
R_prim <- "TACCTTGTTACGACTT" #reverse primer

madison_rain_barrel <- file.path("madison_rain_barrel.fastq")
     
#Basic statistics of raw contigs
summary(nchar(getSequences(madison_rain_barrel)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    204    1337    1443    1349    1482    2671 
 
#Remove forward and reverse primers. This step also ensures that the contigs that remain after filtering are the full length 16S gene. 
removePrimers(madison_rain_barrel, "madison_rain_barrel_noprimers.fastq", primer.fwd=F_prim, primer.rev=rc(R_prim), compress=FALSE, max.mismatch=2, allow.indels=FALSE, verbose=TRUE)
2 sequences out of 15140 are being reverse-complemented.
Read in 15140, output 4673 (30.9%) filtered sequences.

#Basic statistics of filtered contigs (full length 16S, primers trimmed) 
madison_rain_barrel_noprimers <- file.path("madison_rain_barrel_noprimers.fastq")
summary(nchar(getSequences(madison_rain_barrel_noprimers)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1375    1436    1446    1443    1456    1496 

```
[Filtered contigs fastq file](https://github.com/nico-chung/loopseq/raw/master/data/madison_rain_barrel_noprimers.rar)

#### RivaPond

Load fastq file containing raw contigs as downloaded from the Loop Genomics cloud platform.

[Raw contigs fastq file](https://github.com/nico-chung/loopseq/raw/master/data/madison_riva_pond.rar)

```
#R
library(dada2)
library(ShortRead)

F_prim <- "AGAGTTTGATCMTGGCTCAG" #forward primer
R_prim <- "TACCTTGTTACGACTT" #reverse primer

madison_riva_pond <- file.path("madison_riva_pond.fastq")

#Basic statistics of raw contigs
summary(nchar(getSequences(madison_riva_pond)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    242    1090    1412    1253    1448    2036 
    
#Remove forward and reverse primers. This step also ensures that the contigs that remain after filtering are the full length 16S gene. 
removePrimers(madison_riva_pond, "madison_riva_pond_noprimers.fastq", primer.fwd=F_prim, primer.rev=rc(R_prim), compress=FALSE, max.mismatch=2, allow.indels=FALSE, verbose=TRUE)
7 sequences out of 14352 are being reverse-complemented.
Read in 14352, output 3808 (26.5%) filtered sequences.

#Basic statistics of filtered contigs (full length 16S, primers trimmed) 
madison_riva_pond_noprimers <- file.path("madison_riva_pond_noprimers.fastq")
summary(nchar(getSequences(madison_riva_pond_noprimers)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1350    1410    1419    1428    1449    2000 
 
#Still some contigs are too long. Filter for length.
filterAndTrim(madison_riva_pond_noprimers, "madison_riva_pond_noprimers_lengthfilter.fastq", minLen=1400, maxLen=1600, verbose=TRUE, compress=FALSE)
Overwriting file:C:\Users\Nico\Documents\workdir\madison_riva_pond_noprimers_lengthfilter.fastq
Read in 3808, output 3564 (93.6%) filtered sequences.

#Basic statistics of filtered contigs (full length 16S, primers trimmed, length filtered)
madison_riva_pond_noprimers_lengthfilter <- file.path("madison_riva_pond_noprimers_lengthfilter.fastq")
summary(nchar(getSequences(madison_riva_pond_noprimers_lengthfilter)))
```
[Filtered contigs fastq file](https://github.com/nico-chung/loopseq/raw/master/data/madison_riva_pond_noprimers_lengthfilter.rar)



## Classification with LCA Classifier

Taxonomic classification with Kraken 2 outputs reports which show a hierarchical LCA scheme. While most full length 16S contigs can be classified to species level, some contigs can only be classifed to higher taxonomic levels. To obtain species abundances, the companion tool Bracken is run to redistribute all contigs to species level. This redistribution is limited to lineages with at least one species identified. While this method may not fully represent the species diversity in the sample, it does allow a comparison of relative abundances of species that can be classified.

#### Aquarium
```
kraken2 --db Kraken2/standard_DB --threads 24 --out madison_aquarium.kraken --report madison_aquarium.kreport madison_aquarium_noprimers_lengthfilter.fastq 
Loading database information... done.
6428 sequences (9.30 Mbp) processed in 0.440s (876.8 Kseq/m, 1268.39 Mbp/m).
  6428 sequences classified (100.00%)
  0 sequences unclassified (0.00%)
```
[Kraken 2 report](https://github.com/nico-chung/loopseq/raw/master/data/madison_aquarium.kreport)
```
bracken -d Kraken2/standard_DB -i madison_aquarium.kreport -o madison_aquarium.bracken -r 1500 -l S -t 1
 >> Checking for Valid Options...
 >> Running Bracken 
      >> python src/est_abundance.py -i madison_aquarium.kreport -o madison_aquarium.bracken -k Kraken2/standard_DB/database1500mers.kmer_distrib -l S -t 1
PROGRAM START TIME: 12-31-2019 05:36:30
BRACKEN SUMMARY (Kraken report: madison_aquarium.kreport)
    >>> Threshold: 1 
    >>> Number of species in sample: 118 
	  >> Number of species with reads > threshold: 118 
	  >> Number of species with reads < threshold: 0 
    >>> Total reads in sample: 6428
	  >> Total reads kept at species level (reads > threshold): 6129
	  >> Total reads discarded (species reads < threshold): 0
	  >> Reads distributed: 68
	  >> Reads not distributed (eg. no species above threshold): 231
	  >> Unclassified reads: 0
BRACKEN OUTPUT PRODUCED: madison_aquarium.bracken
PROGRAM END TIME: 12-31-2019 05:36:31
  Bracken complete.

```
[Bracken report](https://github.com/nico-chung/loopseq/raw/master/data/madison_aquarium.bracken)


#### FeynmanPond
```
kraken2 --db Kraken2/standard_DB --threads 24 --out madison_feynman_pond.kraken --report madison_feynman_pond.kreport madison_feynman_pond_noprimers.fastq 
Loading database information... done.
8251 sequences (11.89 Mbp) processed in 0.516s (959.3 Kseq/m, 1382.90 Mbp/m).
  8251 sequences classified (100.00%)
  0 sequences unclassified (0.00%)
```
[Kraken 2 report](https://github.com/nico-chung/loopseq/raw/master/data/madison_feynman_pond.kreport)
```
bracken -d Kraken2/standard_DB -i madison_feynman_pond.kreport -o madison_feynman_pond.bracken -r 1500 -l S -t 1
 >> Checking for Valid Options...
 >> Running Bracken 
      >> python src/est_abundance.py -i madison_feynman_pond.kreport -o madison_feynman_pond.bracken -k Kraken2/standard_DB/database1500mers.kmer_distrib -l S -t 1
PROGRAM START TIME: 12-31-2019 05:45:10
BRACKEN SUMMARY (Kraken report: madison_feynman_pond.kreport)
    >>> Threshold: 1 
    >>> Number of species in sample: 332 
	  >> Number of species with reads > threshold: 332 
	  >> Number of species with reads < threshold: 0 
    >>> Total reads in sample: 8251
	  >> Total reads kept at species level (reads > threshold): 5870
	  >> Total reads discarded (species reads < threshold): 0
	  >> Reads distributed: 496
	  >> Reads not distributed (eg. no species above threshold): 1885
	  >> Unclassified reads: 0
BRACKEN OUTPUT PRODUCED: madison_feynman_pond.bracken
PROGRAM END TIME: 12-31-2019 05:45:10
  Bracken complete.
```
[Bracken report](https://github.com/nico-chung/loopseq/raw/master/data/madison_feynman_pond.bracken)


#### RainBarrel
```
kraken2 --db Kraken2/standard_DB --threads 24 --out madison_rain_barrel.kraken --report madison_rain_barrel.kreport madison_rain_barrel_noprimers.fastq
Loading database information... done.
4673 sequences (6.74 Mbp) processed in 0.274s (1022.6 Kseq/m, 1475.63 Mbp/m).
  4673 sequences classified (100.00%)
  0 sequences unclassified (0.00%)
```
[Kraken 2 report](https://github.com/nico-chung/loopseq/raw/master/data/madison_rain_barrel.kreport)
```
bracken -d Kraken2/standard_DB -i madison_rain_barrel.kreport -o madison_rain_barrel.bracken -r 1500 -l S -t 1
 >> Checking for Valid Options...
 >> Running Bracken 
      >> python src/est_abundance.py -i madison_rain_barrel.kreport -o madison_rain_barrel.bracken -k Kraken2/standard_DB/database1500mers.kmer_distrib -l S -t 1
PROGRAM START TIME: 12-31-2019 06:28:32
BRACKEN SUMMARY (Kraken report: madison_rain_barrel.kreport)
    >>> Threshold: 1 
    >>> Number of species in sample: 165 
	  >> Number of species with reads > threshold: 165 
	  >> Number of species with reads < threshold: 0 
    >>> Total reads in sample: 4673
	  >> Total reads kept at species level (reads > threshold): 2336
	  >> Total reads discarded (species reads < threshold): 0
	  >> Reads distributed: 1489
	  >> Reads not distributed (eg. no species above threshold): 848
	  >> Unclassified reads: 0
BRACKEN OUTPUT PRODUCED: madison_rain_barrel.bracken
PROGRAM END TIME: 12-31-2019 06:28:32
  Bracken complete.

```
[Bracken report](https://github.com/nico-chung/loopseq/raw/master/data/madison_rain_barrel.bracken)

#### RivaPond
```
kraken2 --db Kraken2/standard_DB --threads 24 --out madison_riva_pond.kraken --report madison_riva_pond.kreport madison_riva_pond_noprimers_lengthfilter.fastq 
Loading database information... done.
3564 sequences (5.10 Mbp) processed in 0.339s (631.5 Kseq/m, 903.08 Mbp/m).
  3564 sequences classified (100.00%)
  0 sequences unclassified (0.00%)
```
[Kraken 2 report](https://github.com/nico-chung/loopseq/raw/master/data/madison_riva_pond.kreport)
```
bracken -d Kraken2/standard_DB -i madison_riva_pond.kreport -o madison_riva_pond.bracken -r 1500 -l S -t 1
 >> Checking for Valid Options...
 >> Running Bracken 
      >> python src/est_abundance.py -i madison_riva_pond.kreport -o madison_riva_pond.bracken -k Kraken2/standard_DB/database1500mers.kmer_distrib -l S -t 1
PROGRAM START TIME: 12-31-2019 06:16:39
BRACKEN SUMMARY (Kraken report: madison_riva_pond.kreport)
    >>> Threshold: 1 
    >>> Number of species in sample: 166 
	  >> Number of species with reads > threshold: 166 
	  >> Number of species with reads < threshold: 0 
    >>> Total reads in sample: 3564
	  >> Total reads kept at species level (reads > threshold): 2507
	  >> Total reads discarded (species reads < threshold): 0
	  >> Reads distributed: 87
	  >> Reads not distributed (eg. no species above threshold): 970
	  >> Unclassified reads: 0
BRACKEN OUTPUT PRODUCED: madison_riva_pond.bracken
PROGRAM END TIME: 12-31-2019 06:16:39
  Bracken complete.
```
[Bracken report](https://github.com/nico-chung/loopseq/raw/master/data/madison_riva_pond.bracken)


## Diversity Analysis 

The number of species identified and top 3 species by abundance.

|Sample     |Species count|First                                   |Second                                 |Third                                   |
|-----------|-------------|----------------------------------------|---------------------------------------|----------------------------------------|
|Aquarium   |118          |Limnohabitans sp. 63ED37-2 (0.55)       |Polynucleobacter necessarius (0.12)    |Novosphingobium ginsenosidimutans (0.05)|
|FeynmanPond|332          |Candidatus Planktophila limnetica (0.09)|Arachidicoccus ginsenosidivorans (0.09)|Limnohabitans sp. 103DPR2 (0.05)        |
|RainBarrel |165          |Steroidobacter denitrificans (0.25)     |Pseudomonas sp. SWI6 (0.11)            |Planctomycetes bacterium Pan189 (0.08)  |
|RivaPond   |166          |Scytonema sp. HK-05 (0.09)              |Niabella soli (0.07)                   |Planctomycetes bacterium Pan97 (0.06)   |

Combine all bracken reports into a contingency table, convert to biom file and run plot rank abundance.

[Madison water combined](https://github.com/nico-chung/loopseq/raw/master/data/madison_water_combined.txt)

```
biom convert -i madison_water_combined.txt -o madison_water_combined.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
plot_rank_abundance_graph.py -i madison_water_combined.biom -s '*' -x -f png -v -o madison_water_combined_rank_abundance.png
```
![](https://github.com/nico-chung/loopseq/blob/master/figures/madison_water_combined_rank_abundance.png)

The Aquarium sample is dominated by fewer high abundance species while the FeynmanPond sample is most diverse and with a more even species distribution.


```
#R
library(dplyr)
library(ggplot2)
library(RColorBrewer)

#Melt and add sample info to madison_water_combined.txt
data <- read.csv("madison_water_combined.csv", header = TRUE, stringsAsFactors = FALSE)
data <- data %>% arrange(Sample, Abundance)

x <- ggplot(data, aes(x=Sample, y=Abundance, fill=factor(Species),group=Abundance)) + 
     geom_bar(stat="identity",position="stack") + 
     theme(legend.position="none",
           axis.text.x = element_text(size=10,angle=90,vjust=0))
	   
cols <- colorRampPalette(brewer.pal(12, "Paired")) 
myPal <- cols(length(unique(data$Species)))

x + scale_fill_manual(values = sample(myPal))

newpath <- paste0("madison_water_combined_barchart.png")
ggsave(filename=newpath,width=3,height=8,units="cm",dpi=600)

```
<img src="https://github.com/nico-chung/loopseq/blob/master/figures/madison_water_combined_barchart.png" width="500">


&nbsp; 

