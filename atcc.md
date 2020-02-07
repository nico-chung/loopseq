# ATCC Mock Human Communities Analysis 

The two ATCC mock human communities were analyzed using two methods: 
- [Mapping LoopSeq 16S contigs to whole genome strain references](https://github.com/nico-chung/loopseq/blob/master/atcc.md#mapping-to-whole-genome-references)
- [Strain detection using LCA classifier](https://github.com/nico-chung/loopseq/blob/master/atcc.md#strain-level-detection-using-lca-classifier)

## Experimental Protocol
 
[Fill in experimental details here: Sample prep, kits used, platforms used, yield per sample, etc.] 
 
## Data Preparation

### Mock Community Description

Although the Zymo mock community is widely used and is arguably the common standard in microbiome studies, many of the references that are provided do not match the strains in the sample. Zymo Research also does not offer the extensive supporting materials that the American Type Culture Association (ATCC) offers through their mock community samples (BEI Resources, USA). Many of the species in the ATCC mock communities have exact strain genome references along with assembly details and genome annotations. Where exact reference strains are not provided, strains names allows researchers to look for identical or similar references in common genome repositories (eg. NCBI RefSeq CG database).   

Two ATCC mock communities representing human-derived samples were analyzed. These two mock communities demonstrates the utility of applying LoopSeq for human site-specific studies. The provided exact strain names also allows for a more fine-grained assessment of the level of taxonomic resolution that can be achieved through accurate full-length 16S sequences. Two mock communities were analyzed: ATCC-gut and ATCC-oral. Six out of the 12 genome references are provided by the ATCC for ATCC-gut, while the remaining six have equivalent NCBI RefSeq assemblies. All six of genome references are provided by the ATCC for ATCC-oral. 

&nbsp; 

**[Gut Microbiome Genomic Mix (MSA-1006)](https://www.atcc.org/products/all/MSA-1006.aspx) ("ATCC-gut")**

|Species                                  |Strain         |Genomic DNA %|ATCC genome reference|Reference                                               |RefSeq assembly| RefSeq assembly strain           |16S copy information                                                                                                                                          |16S copies|
|-----------------------------------------|---------------|-------------|---------------------|-------------------------------------------------------------|---------------|----------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|----------|
|Bacteriodes fragilis                     |ATCC 25285 |8.3          |No                   |[link](https://www.ncbi.nlm.nih.gov/assembly/GCF_000025985.1)|Yes            |NCTC 9343 (ATCC 25285)|Yes                                                                                                                                                           |6         |
|Bacteroides vulgatus                     |ATCC 8482  |8.3          |No                   |[link](https://www.ncbi.nlm.nih.gov/assembly/GCF_000012825.1)|Yes            |ATCC 8482             |Yes                                                                                                                                                           |7         |
|Clostridioides difficile                 |ATCC 9689  |8.3          |No                   |[link](https://www.ncbi.nlm.nih.gov/assembly/GCF_001077535.1)|Yes            |ATCC 9689             |Yes                                                                                                                                                           |12        |
|Enterobacter cloacae                     |ATCC 13047 |8.3          |No                   |[link](https://www.ncbi.nlm.nih.gov/assembly/GCF_000025565.1)|Yes            |ATCC 13047            |[link](https://rrndb.umms.med.umich.edu/genomes/GCF_000025565.1?query=Enterobacter+cloacae&button=Search+taxonomy&name_type=13&per_page=&search_type=taxonomy)|8         |
|Escherichia coli                         |ATCC 700926 |8.3          |No                   |[link](https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2)|Yes            |MG1655 (ATCC 700926)  |[link](https://rrndb.umms.med.umich.edu/genomes/GCF_000005845.2?query=escherichia+coli+&button=Search+taxonomy&name_type=13&per_page=&search_type=taxonomy)   |7         |
|Salmonella enterica subsp. enterica      |ATCC 9150  |8.3          |No                   |[link](https://www.ncbi.nlm.nih.gov/assembly/GCF_000011885.1)|Yes            |ATCC 9150             |Yes                                                                                                                                                           |7         |
|Bifidobacterium adolescentis             |ATCC 15703 |8.3          |Yes                  |[link](https://genomes.atcc.org/genomes/90eb97d11e4b445f)    |               |                      |                                                                                                                                                              |5         |
|Enterococcus faecalis                    |ATCC 700802 |8.3          |Yes                  |[link](https://genomes.atcc.org/genomes/535b47f8d8a142a3)    |               |                      |                                                                                                                                                              |4         |
|Lactobacillus plantarum                  |ATCC BAA-793|8.3          |Yes                  |[link](https://genomes.atcc.org/genomes/daded83eca2b4d03)    |               |                      |                                                                                                                                                              |5         |
|Helicobacter pylori                      |ATCC 700392 |8.3          |Yes                  |[link](https://genomes.atcc.org/genomes/9038b5a9e94245e8)    |               |                      |                                                                                                                                                              |2         |
|Yersinia enterocolitica                  |ATCC 27729  |8.3          |Yes                  |[link](https://genomes.atcc.org/genomes/099e5acebc284d19)    |               |                      |                                                                                                                                                              |7         |
|Fusobacterium nucleatum subsp. nucleatum |ATCC 25586  |8.3          |Yes                  |[link](https://genomes.atcc.org/genomes/d7727bb49ac34e11)    |               |                      |                                                                                                                                                              |5         |

&nbsp; 


**[Oral Microbiome Genomic Mix (MSA-1004)](https://www.atcc.org/products/all/MSA-1004.aspx) ("ATCC-oral")**

|Species                                  |Strain       |Genomic DNA %|ATCC genome reference|Reference                                                |16S copies|
|-----------------------------------------|-------------|-------------|---------------------|---------------------------------------------------------|----------|
|Schaalia odontolytica                    |ATCC 17982|16.7         |Yes                  |[link](https://genomes.atcc.org/genomes/6823ab7475dd4769)|3         |
|Prevotella melaninogenica                |ATCC 25845|16.7         |Yes                  |[link](https://genomes.atcc.org/genomes/7fc50ab8f518490d)|4         |
|Fusobacterium nucleatum subsp. nucleatum |ATCC 25586|16.7         |Yes                  |[link](https://genomes.atcc.org/genomes/d7727bb49ac34e11)|5         |
|Streptococcus mitis                      |ATCC 49456|16.7         |Yes                  |[link](https://genomes.atcc.org/genomes/28282a8a3f2a4e92)|4         |
|Veillonella parvula                      |ATCC 17745|16.7         |Yes                  |[link](https://genomes.atcc.org/genomes/3dc68c1d228c4c81)|4         |
|Haemophilus parainfluenzae               |ATCC 33392|16.7         |Yes                  |[link](https://genomes.atcc.org/genomes/8b417db59c2e4292)|6         |

&nbsp; 


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


#### ATCC-gut 

Load fastq file containing raw contigs as downloaded from the Loop Genomics cloud platform.

[Raw contigs fastq file](https://github.com/nico-chung/loopseq/blob/master/data/atcc_gut.rar)

```
#R
library(dada2)
library(ShortRead)

F_prim <- "AGAGTTTGATCMTGGCTCAG" #forward primer
R_prim <- "TACCTTGTTACGACTT" #reverse primer

atcc_gut <- file.path("atcc_gut.fastq")

#Basic statistics of raw contigs
summary(nchar(getSequences(atcc_gut)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    259    1436    1462    1376    1489    2242 


#Distribution of raw contigs
hist(nchar(getSequences(atcc_gut)), 100) 
```
![](https://github.com/nico-chung/loopseq/blob/master/figures/atcc_gut_rawcontigs.png)
```
#Remove forward and reverse primers. This step also ensures that the contigs that remain after filtering are the full length 16S gene. 
removePrimers(atcc_gut, "atcc_gut_noprimers.fastq", primer.fwd=F_prim, primer.rev=rc(R_prim),  compress=FALSE, max.mismatch=2, allow.indels=FALSE, verbose=TRUE)
Multiple matches to the primer(s) in some sequences. Using the longest possible match.
Multiple matches to the primer(s) in some reverse-complement sequences. Using the longest possible match.
86 sequences out of 20027 are being reverse-complemented.
Read in 20027, output 9087 (45.4%) filtered sequences.

#Basic statistics of filtered contigs (full length 16S, primers trimmed) 
atcc_gut_noprimers <- file.path("atcc_gut_noprimers.fastq")
summary(nchar(getSequences(atcc_gut_noprimers)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1416    1426    1451    1451    1464    1585 

#Distribution of filtered contigs (full length 16S, primers trimmed)
hist(nchar(getSequences(atcc_gut_noprimers)), 10) 
```
[Filtered contigs fastq file](https://github.com/nico-chung/loopseq/blob/master/data/atcc_gut_noprimers.rar)

![](https://github.com/nico-chung/loopseq/blob/master/figures/atcc_gut_noprimers.png)

&nbsp; 


#### ATCC-oral

Load fastq file containing raw contigs as downloaded from the Loop Genomics cloud platform.

[Raw contigs fastq file](https://github.com/nico-chung/loopseq/blob/master/data/atcc_oral.rar)

```
#R
library(dada2)
library(ShortRead)

F_prim <- "AGAGTTTGATCMTGGCTCAG" #forward primer
R_prim <- "TACCTTGTTACGACTT" #reverse primer

atcc_oral <- file.path("atcc_oral.fastq")

#Basic statistics of raw contigs
summary(nchar(getSequences(atcc_oral)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    205    1444    1478    1372    1498    2255   
    
#Distribution of raw contigs
hist(nchar(getSequences(atcc_oral)), 100) 
```
![](https://github.com/nico-chung/loopseq/blob/master/figures/atcc_oral_rawcontigs.png)
```
#Remove forward and reverse primers. This step also ensures that the contigs that remain after filtering are the full length 16S gene. 
removePrimers(atcc_oral, "atcc_oral_noprimers.fastq"), primer.fwd=F_prim, primer.rev=rc(R_prim),  compress=FALSE, max.mismatch=2, allow.indels=FALSE, verbose=TRUE)
Multiple matches to the primer(s) in some sequences. Using the longest possible match.
193 sequences out of 33222 are being reverse-complemented.
Read in 33222, output 16392 (49.3%) filtered sequences.

#Basic statistics of filtered contigs (full length 16S, primers trimmed) 
atcc_oral_noprimers <- file.path("atcc_oral_noprimers.fastq")
summary(nchar(getSequences(atcc_oral_noprimers)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    300    1442    1462    1460    1469    2219 

#The max contig length after primer trimming is 2219 bp. 16S genes shouldn't be that long. There are also contigs that are too short. These may represent experimental or assembly errors. These are filtered out.
filterAndTrim(atcc_oral_noprimers, "atcc_oral_noprimers_lengthfilter.fastq", minLen=1400, maxLen=1600, verbose=TRUE, compress=FALSE)
Read in 16392, output 16173 (98.7%) filtered sequences.

atcc_oral_noprimers_lengthfilter <- file.path("atcc_oral_noprimers_lengthfilter.fastq")
summary(nchar(getSequences(atcc_oral_noprimers_lengthfilter)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1409    1442    1462    1460    1469    1595 
   
#Distribution of filtered contigs (full length 16S, primers trimmed, length filtered)
hist(nchar(getSequences(atcc_oral_noprimers_lengthfilter)), 10) 

```
[Filtered contigs fastq file](https://github.com/nico-chung/loopseq/blob/master/data/atcc_oral_noprimers_lengthfilter.rar)

![](https://github.com/nico-chung/loopseq/blob/master/figures/atcc_oral_noprimers_lengthfilter.png)



## Mapping to Whole Genome References

### ATCC-gut 

Concatenate all bacterial whole genome fasta files together and build a Bowtie2 index with the whole genome references and map the filtered ATCC-gut contigs.

[ATCC-gut genome references](https://github.com/nico-chung/loopseq/blob/master/data/atcc_gut_genome_references.rar)

```
bowtie2 --very-sensitive -p 40 --un atcc_gut.un.fastq -x atcc_gut_ref/atcc_gut_ref -q atcc_gut_noprimers.fastq -S atcc_gut.sam
9087 reads; of these:
  9087 (100.00%) were unpaired; of these:
    3 (0.03%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    9084 (99.97%) aligned >1 times
99.97% overall alignment rate

#Extract names of hits and count number of each unique hit
samtools view -F 4 atcc_gut.sam | cut -f3 | sort | uniq -c 
   1067 Bacteroides_fragilis_NCTC_9343
    968 Bacteroides_vulgatus_ATCC_8482
     20 bifidobacterium_adolescentis_atcc_15703_contig1
     11 bifidobacterium_adolescentis_atcc_15703_contig2
   1966 Clostridioides_difficile_ATCC_9689
    692 Enterobacter_cloacae_ATCC_13047
    582 enterococcus_faecalis_atcc_700802_contig1
    735 Escherichia_coli_MG1655
    894 fusobacterium_nucleatum_subsp_nucleatum_atcc_25586
    531 helicobacter_pylori_atcc_700392
    556 lactobacillus_plantarum_atcc_baa-793_contig1
    391 Salmonella_enterica_ATCC_9150
    671 yersinia_enterocolitica_atcc_27729_contig1
```

|Species                                  |Expected contig %|Mapped contig %|
|-----------------------------------------|-----------------|---------------|
|Bacteriodes fragilis                     |8.0              |11.7           |
|Bacteroides vulgatus                     |9.3              |10.7           |
|Clostridioides difficile                 |16.0             |21.6           |
|Enterobacter cloacae                     |10.7             |7.6            |
|Escherichia coli                         |9.3              |8.1            |
|Salmonella enterica subsp. enterica      |9.3              |4.3            |
|Bifidobacterium adolescentis             |6.7              |0.3            |
|Enterococcus faecalis                    |5.3              |6.4            |
|Lactobacillus plantarum                  |6.7              |6.1            |
|Helicobacter pylori                      |2.7              |5.8            |
|Yersinia enterocolitica                  |9.3              |7.4            |
|Fusobacterium nucleatum subsp. nucleatum |6.7              |9.8            |

The proportion of species 16S ratios are converted from 16S copy numbers. All 12 species are observed in approximately the expected ratios, except *B. adolescentis* which is much lower than expected. 


&nbsp; 

### ATCC-oral

Concatenate all bacterial whole genome fasta files together and build a Bowtie2 index with the whole genome references and map the filtered ATCC-oral contigs.

[ATCC-oral genome references](https://github.com/nico-chung/loopseq/blob/master/data/atcc_oral_genome_references.rar)

```
bowtie2 --very-sensitive -p 40 --un atcc_oral.un.fastq -x atcc_oral_ref/atcc_oral_ref -q atcc_oral_noprimers_lengthfilter.fastq -S atcc_oral.sam
16173 reads; of these:
  16173 (100.00%) were unpaired; of these:
    1 (0.01%) aligned 0 times
    0 (0.00%) aligned exactly 1 time
    16172 (99.99%) aligned >1 times
99.99% overall alignment rate

#Extract names of hits and count number of each unique hit
samtools view -F 4 atcc_oral.sam | cut -f3 | sort | uniq -c 
   4542 fusobacterium_nucleatum_subsp_nucleatum_atcc_25586
   4144 haemophilus_parainfluenzae_atcc_33392_contig1
   1083 prevotella_melaninogenica_atcc_25845_contig1
   1129 prevotella_melaninogenica_atcc_25845_contig2
    623 schaalia_odontolytica_atcc_17982
   2485 streptococcus_mitis_atcc_49456
   2166 veillonella_parvula_atcc_17745
```
|Species                                  |Expected contig %|Mapped contig %|
|-----------------------------------------|-----------------|---------------|
|Schaalia odontolytica                    |11.5             |3.9            |
|Prevotella melaninogenica                |15.4             |13.7           |
|Fusobacterium nucleatum subsp. nucleatum |19.2             |28.1           |
|Streptococcus mitis                      |15.4             |15.4           |
|Veillonella parvula                      |15.4             |13.4           |
|Haemophilus parainfluenzae               |23.1             |25.6           |

&nbsp; 


## Strain-level Detection Using LCA Classifier

### ATCC-gut 

```
kraken2 --db Kraken2/standard_DB --threads 24 --out atcc_gut.kraken --report atcc_gut.kreport atcc_gut_noprimers.fastq 
Loading database information... done.
9087 sequences (13.19 Mbp) processed in 0.384s (1418.5 Kseq/m, 2058.48 Mbp/m).
  9087 sequences classified (100.00%)
  0 sequences unclassified (0.00%)
```

```
100.00	9087	0	R	1	root
100.00	9087	0	R1	131567	  cellular organisms
100.00	9087	0	D	2	    Bacteria
 34.52	3137	0	D1	1783272	      Terrabacteria group
 34.18	3106	0	P	1239	        Firmicutes
 21.65	1967	0	C	186801	          Clostridia
 21.65	1967	0	O	186802	            Clostridiales
 21.64	1966	0	F	186804	              Peptostreptococcaceae
 21.64	1966	0	G	1870884	                Clostridioides
 21.64	1966	1965	S	1496	                  Clostridioides difficile
  0.01	1	1	S1	699037	                    Clostridioides difficile M68
  0.01	1	0	F	31979	              Clostridiaceae
  0.01	1	0	G	1485	                Clostridium
  0.01	1	1	S	1502	                  Clostridium perfringens
 12.52	1138	0	C	91061	          Bacilli
 12.52	1138	0	O	186826	            Lactobacillales
  6.40	582	0	F	81852	              Enterococcaceae
  6.40	582	0	G	1350	                Enterococcus
  6.40	582	582	S	1351	                  Enterococcus faecalis
  6.12	556	0	F	33958	              Lactobacillaceae
  6.12	556	229	G	1578	                Lactobacillus
  3.57	324	324	S	1590	                  Lactobacillus plantarum
  0.01	1	1	S	1613	                  Lactobacillus fermentum
  0.01	1	0	S	109790	                  Lactobacillus jensenii
  0.01	1	1	S1	525329	                    Lactobacillus jensenii JV-V16
  0.01	1	1	S	1610	                  Lactobacillus coryniformis
  0.01	1	0	C	909932	          Negativicutes
  0.01	1	0	O	1843489	            Veillonellales
  0.01	1	0	F	31977	              Veillonellaceae
  0.01	1	0	G	29465	                Veillonella
  0.01	1	1	S	29466	                  Veillonella parvula
  0.34	31	0	P	201174	        Actinobacteria
  0.34	31	0	C	1760	          Actinobacteria
  0.34	31	0	O	85004	            Bifidobacteriales
  0.34	31	0	F	31953	              Bifidobacteriaceae
  0.34	31	0	G	1678	                Bifidobacterium
  0.34	31	29	S	1680	                  Bifidobacterium adolescentis
  0.02	2	2	S1	367928	                    Bifidobacterium adolescentis ATCC 15703
 33.23	3020	0	P	1224	      Proteobacteria
 27.39	2489	0	C	1236	        Gammaproteobacteria
 27.39	2489	0	O	91347	          Enterobacterales
 20.00	1817	1503	F	543	            Enterobacteriaceae
  2.45	223	0	G	561	              Escherichia
  2.45	223	222	S	562	                Escherichia coli
  0.01	1	0	S1	1072458	                  Escherichia coli O7:K1
  0.01	1	1	S2	1072459	                    Escherichia coli O7:K1 str. CE10
  0.77	70	0	G	590	              Salmonella
  0.77	70	56	S	28901	                Salmonella enterica
  0.14	13	13	S1	59201	                  Salmonella enterica subsp. enterica
  0.01	1	0	S1	59202	                  Salmonella enterica subsp. salamae
  0.01	1	1	S2	2500152	                    Salmonella enterica subsp. salamae serovar 42:r:-
  0.14	13	1	G	570	              Klebsiella
  0.09	8	7	S	573	                Klebsiella pneumoniae
  0.01	1	1	S1	72407	                  Klebsiella pneumoniae subsp. pneumoniae
  0.01	1	1	S	244366	                Klebsiella variicola
  0.01	1	1	S	1134687	                Klebsiella michiganensis
  0.01	1	1	S	1463165	                Klebsiella quasipneumoniae
  0.01	1	1	S	2153354	                Klebsiella sp. WCHKl090001
  0.03	3	0	G	620	              Shigella
  0.01	1	1	S	621	                Shigella boydii
  0.01	1	1	S	622	                Shigella dysenteriae
  0.01	1	1	S	623	                Shigella flexneri
  0.02	2	0	F1	191675	              unclassified Enterobacteriaceae
  0.02	2	0	F2	84563	                ant, tsetse, mealybug, aphid, etc. endosymbionts
  0.01	1	0	F3	146507	                  aphid secondary symbionts
  0.01	1	1	S	1199245	                    secondary endosymbiont of Ctenarytaina eucalypti
  0.01	1	0	G	1906660	                  Candidatus Mikella
  0.01	1	1	S	1778264	                    Candidatus Mikella endobia
  0.01	1	0	G	547	              Enterobacter
  0.01	1	0	G1	354276	                Enterobacter cloacae complex
  0.01	1	1	S	158836	                  Enterobacter hormaechei
  0.01	1	0	G	544	              Citrobacter
  0.01	1	0	G1	1344959	                Citrobacter freundii complex
  0.01	1	1	S	2066049	                  Citrobacter freundii complex sp. CFNIH2
  0.01	1	0	G	1780190	              Izhakiella
  0.01	1	1	S	2579935	                Izhakiella sp. KSNA2
  7.40	672	0	F	1903411	            Yersiniaceae
  7.40	672	3	G	629	              Yersinia
  7.36	669	669	S	630	                Yersinia enterocolitica
  5.84	531	0	P1	68525	        delta/epsilon subdivisions
  5.84	531	0	C	29547	          Epsilonproteobacteria
  5.84	531	0	O	213849	            Campylobacterales
  5.84	531	0	F	72293	              Helicobacteraceae
  5.84	531	0	G	209	                Helicobacter
  5.84	531	530	S	210	                  Helicobacter pylori
  0.01	1	1	S1	1163741	                    Helicobacter pylori Shi169
 22.39	2035	0	D1	1783270	      FCB group
 22.39	2035	0	D2	68336	        Bacteroidetes/Chlorobi group
 22.39	2035	0	P	976	          Bacteroidetes
 22.39	2035	0	C	200643	            Bacteroidia
 22.39	2035	0	O	171549	              Bacteroidales
 22.39	2035	0	F	815	                Bacteroidaceae
 22.39	2035	0	G	816	                  Bacteroides
 11.73	1066	1066	S	817	                    Bacteroides fragilis
 10.65	968	271	S	821	                    Bacteroides vulgatus
  7.67	697	697	S1	435590	                      Bacteroides vulgatus ATCC 8482
  0.01	1	1	S	820	                    Bacteroides uniformis
  9.84	894	0	P	32066	      Fusobacteria
  9.84	894	0	C	203490	        Fusobacteriia
  9.84	894	0	O	203491	          Fusobacteriales
  9.84	894	0	F	203492	            Fusobacteriaceae
  9.84	894	0	G	848	              Fusobacterium
  9.84	894	0	S	851	                Fusobacterium nucleatum
  9.84	894	379	S1	76856	                  Fusobacterium nucleatum subsp. nucleatum
  5.67	515	515	S2	190304	                    Fusobacterium nucleatum subsp. nucleatum ATCC 25586
  0.01	1	0	D1	1783257	      PVC group
  0.01	1	0	P	74201	        Verrucomicrobia
  0.01	1	0	C	203494	          Verrucomicrobiae
  0.01	1	0	O	48461	            Verrucomicrobiales
  0.01	1	0	F	1647988	              Akkermansiaceae
  0.01	1	0	G	239934	                Akkermansia
  0.01	1	1	S	239935	                  Akkermansia muciniphila
```

|Species                                  |Strain      |Species detected|Strain detected|Name of strain detected                                                 |
|-----------------------------------------|------------|----------------|---------------|------------------------------------------------------------------------|
|Bacteriodes fragilis                     |ATCC 25285  |Yes             |No             |                                                                        |
|Bacteroides vulgatus                     |ATCC 8482   |Yes             |Yes            |Bacteroides vulgatus ATCC 8482                                          |
|Clostridioides difficile                 |ATCC 9689   |Yes             |Yes            |Clostridioides difficile M68                                            |
|Enterobacter cloacae                     |ATCC 13047  |No              |               |                                                                        |
|Escherichia coli                         |ATCC 700926 |Yes             |Yes            |Escherichia coli O7:K1                                                  |
|Salmonella enterica subsp. enterica      |ATCC 9150   |Yes             |Yes            |Salmonella enterica subsp. Enterica & Salmonella enterica subsp. Salamae|
|Bifidobacterium adolescentis             |ATCC 15703  |Yes             |Yes            |Bifidobacterium adolescentis ATCC 15703                                 |
|Enterococcus faecalis                    |ATCC 700802 |Yes             |No             |                                                                        |
|Lactobacillus plantarum                  |ATCC BAA-793|Yes             |No             |                                                                        |
|Helicobacter pylori                      |ATCC 700392 |Yes             |Yes            |Helicobacter pylori Shi169                                              |
|Yersinia enterocolitica                  |ATCC 27729  |Yes             |No             |                                                                        |
|Fusobacterium nucleatum subsp. nucleatum |ATCC 25586  |Yes             |Yes            |Fusobacterium nucleatum subsp. nucleatum ATCC 25586                     |

Of the 12 species in the sample, 11 could be detected to species level. Of these 11 detected species, seven could be detected to strain level. Of these seven detected strains, three are exact matches to the ATCC strain designation. 

&nbsp; 


### ATCC-oral

```
kraken2 --db Kraken2/standard_DB --threads 24 --out atcc_oral.kraken --report atcc_oral.kreport atcc_oral_noprimers_lengthfilter.fastq 
Loading database information... done.
16173 sequences (23.61 Mbp) processed in 0.481s (2018.2 Kseq/m, 2945.82 Mbp/m).
  16173 sequences classified (100.00%)
  0 sequences unclassified (0.00%)
```
```
100.00	16173	0	R	1	root
100.00	16173	0	R1	131567	  cellular organisms
100.00	16173	0	D	2	    Bacteria
 32.61	5274	0	D1	1783272	      Terrabacteria group
 28.76	4651	0	P	1239	        Firmicutes
 15.37	2485	0	C	91061	          Bacilli
 15.37	2485	0	O	186826	            Lactobacillales
 15.37	2485	0	F	1300	              Streptococcaceae
 15.37	2485	2450	G	1301	                Streptococcus
  0.10	16	16	S	2420308	                  Streptococcus sp. FDAARGOS_520
  0.02	3	3	S	1346	                  Streptococcus iniae
  0.02	3	3	S	1313	                  Streptococcus pneumoniae
  0.01	2	0	G1	119603	                  Streptococcus dysgalactiae group
  0.01	2	0	S	1334	                    Streptococcus dysgalactiae
  0.01	2	0	S1	119602	                      Streptococcus dysgalactiae subsp. equisimilis
  0.01	2	2	S2	486410	                        Streptococcus dysgalactiae subsp. equisimilis GGS_124
  0.01	1	1	S	2420310	                  Streptococcus sp. FDAARGOS_522
  0.01	1	1	S	1917441	                  Streptococcus ruminantium
  0.01	1	1	S	150055	                  Streptococcus lutetiensis
  0.01	1	1	S	28037	                  Streptococcus mitis
  0.01	1	1	S	1319	                  Streptococcus sp. 'group B'
  0.01	1	0	S	1311	                  Streptococcus agalactiae
  0.01	1	1	S1	1309807	                    Streptococcus agalactiae ILRI005
  0.01	1	1	S	1309	                  Streptococcus mutans
  0.01	1	1	S	1308	                  Streptococcus thermophilus
  0.01	1	1	S	1111760	                  Streptococcus troglodytae
  0.01	1	1	S	1305	                  Streptococcus sanguinis
  0.01	1	0	G1	671232	                  Streptococcus anginosus group
  0.01	1	1	S	1338	                    Streptococcus intermedius
 13.39	2166	0	C	909932	          Negativicutes
 13.39	2166	0	O	1843489	            Veillonellales
 13.39	2166	0	F	31977	              Veillonellaceae
 13.39	2166	0	G	29465	                Veillonella
 13.39	2166	2166	S	29466	                  Veillonella parvula
  3.85	623	0	P	201174	        Actinobacteria
  3.85	623	0	C	1760	          Actinobacteria
  3.85	623	0	O	2037	            Actinomycetales
  3.85	623	0	F	2049	              Actinomycetaceae
  3.80	615	0	G	2529408	                Schaalia
  3.80	615	615	S	1660	                  Schaalia odontolytica
  0.05	8	0	G	1654	                Actinomyces
  0.05	8	8	S	1852377	                  Actinomyces pacaensis
 28.08	4542	0	P	32066	      Fusobacteria
 28.08	4542	0	C	203490	        Fusobacteriia
 28.08	4542	0	O	203491	          Fusobacteriales
 28.08	4542	0	F	203492	            Fusobacteriaceae
 28.08	4542	0	G	848	              Fusobacterium
 28.08	4542	0	S	851	                Fusobacterium nucleatum
 28.08	4542	2056	S1	76856	                  Fusobacterium nucleatum subsp. nucleatum
 15.37	2485	2485	S2	190304	                    Fusobacterium nucleatum subsp. nucleatum ATCC 25586
  0.01	1	1	S2	525283	                    Fusobacterium nucleatum subsp. nucleatum ATCC 23726
 25.62	4144	0	P	1224	      Proteobacteria
 25.62	4144	0	C	1236	        Gammaproteobacteria
 25.62	4144	0	O	135625	          Pasteurellales
 25.62	4144	0	F	712	            Pasteurellaceae
 25.62	4144	0	G	724	              Haemophilus
 25.62	4144	4144	S	729	                Haemophilus parainfluenzae
 13.68	2213	0	D1	1783270	      FCB group
 13.68	2213	0	D2	68336	        Bacteroidetes/Chlorobi group
 13.68	2213	0	P	976	          Bacteroidetes
 13.68	2213	0	C	200643	            Bacteroidia
 13.68	2213	0	O	171549	              Bacteroidales
 13.68	2212	0	F	171552	                Prevotellaceae
 13.68	2212	0	G	838	                  Prevotella
 13.68	2212	2212	S	28132	                    Prevotella melaninogenica
  0.01	1	0	F	815	                Bacteroidaceae
  0.01	1	0	G	816	                  Bacteroides
  0.01	1	1	S	817	                    Bacteroides fragilis
```

|Species                                  |Strain      |Species detected|Strain detected|Name of strain detected                                                 |
|-----------------------------------------|------------|----------------|---------------|------------------------------------------------------------------------|
|Schaalia odontolytica                    |ATCC 17982  |Yes             |No             |                                                                        |
|Prevotella melaninogenica                |ATCC 25845  |Yes             |No             |                                                                        |
|Fusobacterium nucleatum subsp. nucleatum |ATCC 25586  |Yes             |Yes            |Fusobacterium nucleatum subsp. nucleatum ATCC 25586                     |
|Streptococcus mitis                      |ATCC 49456  |Yes             |No             |                                                                        |
|Veillonella parvula                      |ATCC 17745  |Yes             |No             |                                                                        |
|Haemophilus parainfluenzae               |ATCC 33392  |Yes             |No             |                                                                        |

Of the six species in the sample, all six could be detected to species level. Of these six species detected, one could be detected to strain level, which was an exact match to the ATCC strain designation. 


&nbsp; 

