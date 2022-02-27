# Metagenome Assembled Genomes from Churince, Cuatro Cienegas

## Step by step workflow used to assemble the Churince metagenomes to recover MAGs


The following pipeline for the assembly and identification of MAGs was perfomed in three different samples:

1) Churince sediment, previously known to be enriched in _Bordetella_. Hereafter referred to as Bordetella.

2) Mesocosms enrichment experiment. Contains sediment and water samples from the enrichment experiment from 2011 in Lagunita.

3) Microbial mats samples collected by Valerie de Anda for her Ph.D. thesis. The samples corresponded to a time lapse and will be referred to as Mats.

The Mats and Bordetella samples where kindly provided by Valerie de Anda from the backup data of the lab.

The samples corresponding to the metagenome of the mesocosm enrichment experiment were retrieved from the NCBI (https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=ERP104517).

The download and rename of the samples were done with the script `DownloadSamplesFromSRA.sh` in the cluster from FZJ.

```bash
cd /mnt/data/yaxal/CC/data

bash DownloadSamplesFromSRA.sh ListSamplesSRA
```

All samples were located in `/mnt/data/yaxal/CC/data`, and were processed with Trimmomatic v0.39.

```bash
cd /mnt/data/yaxal/CC/data
qsub Trimmomatic.sh

```

### Raw reads and data formating

Once the sequences were quality trimmed, the paired end sequences were "merged" as interleaved sequences. To achieve this, the program `reformat.sh` from the package BBMap was used. In the case of the Bordetella sample, as it was only one, BBMap was run manually, while in the case of the other samples, each was run in a script.

It is important to note, that the mesocosm experiment was subdivided in sediment and water samples, thus each environment it was run independently. 

```bash
cd /mnt/data/yaxal/CC/data/Trimmed
mkdir Interleaved
cd ../..
qsub BBMapMats.sh
qsub BBMapSediment.sh
qsub BBMapWater.sh
```

### Metagenome assembly and MetaBAT binning

For the rest of the Water and MicrobialMats, the difference were the number of tasks run `#$ -t 1:15` for water and `#$ -t 1:12` for the mats, as well as the manifest files `manifest.water` and `manifest.mats`, respectively.

For the assembly of the Churince MAGs', the samples were processed with MEGAHIT v.1.2.9

The following steps were run as a single script containing several processes:

* Merging and compress the sequences for the same site.

* Assemble the genome with MEGAHIT and retain only the contigs longer than 2000 bases.

* Map the assembly to the reads with BWA and generate a BAM file.

* Binning with MetaBAT (https://bitbucket.org/berkeleylab/metabat/src/master/)

```bash
cd /mnt/data/yaxal/CC
qsub PipeAssemblyBordetella.sh
qsub PipeAssemblyMats.sh
qsub PipeAssemblySediment.sh
qsub PipeAssemblyWater.sh
```

To filter the contigs based on their lenght, the script `FilterFastaByLength.pl` was used. 

The evaluation of this first assembly step and binning with MetaBAT, returned a total of 88 bins distributed among the four samples, with the mats containing the greatest diversity with 45 bins, and the Bordetella sample containing only 1 bin. We present the evaluation results only for those with completeness values higher than 40%.

* Bordetella - 1 bin

* Lagunita Sediment - 5 bins

|  Bin Id        |         Marker lineage        |  Completeness | Contamination |  Strain heterogeneity  |
|----------------|-------------------------------|---------------|---------------|------------------------|
|  binsMetaBAT.2 |      k__Bacteria (UID203)     |     59.40     |     13.42     |          6.25          |
|  binsMetaBAT.3 |  f__Lachnospiraceae (UID1286) |     41.25     |      0.12     |          0.00          |

* Lagunita Water - 37 bins

|  Bin Id         |           Marker lineage           |  Completeness |  Contamination |  Strain heterogeneity | 
|-----------------|------------------------------------|---------------|----------------|-----------------------|
|  binsMetaBAT.20 |        k__Bacteria (UID203)        |     100.00    |      137.55    |         54.92         | 
|  binsMetaBAT.24 |     p__Bacteroidetes (UID2591)     |     99.49     |       2.02     |          0.00         | 
|  binsMetaBAT.9  |     p__Bacteroidetes (UID2591)     |     98.02     |       0.99     |         50.00         | 
|  binsMetaBAT.4  |       k__Bacteria (UID2570)        |     96.45     |       0.00     |          0.00         | 
|  binsMetaBAT.2  |       k__Bacteria (UID2982)        |     95.75     |       1.40     |         33.33         | 
|  binsMetaBAT.23 |  c__Alphaproteobacteria (UID3422)  |     95.55     |       1.29     |         85.71         | 
|  binsMetaBAT.13 |     p__Bacteroidetes (UID2605)     |     95.24     |       1.67     |         50.00         | 
|  binsMetaBAT.14 |     p__Bacteroidetes (UID2591)     |     94.43     |       3.71     |          0.00         | 
|  binsMetaBAT.1  |     p__Bacteroidetes (UID2605)     |     94.29     |       0.00     |          0.00         | 
|  binsMetaBAT.31 |       k__Bacteria (UID3187)        |     94.20     |       3.74     |         42.86         | 
|  binsMetaBAT.10 |       k__Bacteria (UID2565)        |     93.55     |       0.00     |          0.00         | 
|  binsMetaBAT.25 |   f__Rhodobacteraceae (UID3340)    |     93.02     |       5.02     |         25.00         | 
|  binsMetaBAT.33 |       k__Bacteria (UID3187)        |     92.80     |       1.79     |          0.00         | 
|  binsMetaBAT.29 |  c__Alphaproteobacteria (UID3305)  |     91.52     |       0.43     |         100.00        | 
|  binsMetaBAT.36 |  c__Betaproteobacteria (UID3971)   |     89.16     |       0.96     |         33.33         | 
|  binsMetaBAT.8  |        k__Bacteria (UID203)        |     84.83     |       1.72     |          0.00         | 
|  binsMetaBAT.17 |  c__Alphaproteobacteria (UID3422)  |     84.31     |       6.43     |         18.92         | 
|  binsMetaBAT.19 |    p__Proteobacteria (UID3887)     |     81.94     |       1.83     |         33.33         | 
|  binsMetaBAT.5  |        k__Bacteria (UID203)        |     75.71     |       0.00     |          0.00         | 
|  binsMetaBAT.26 |   o__Sphingomonadales (UID3310)    |     75.56     |       6.61     |         41.03         | 
|  binsMetaBAT.12 |        k__Bacteria (UID203)        |     73.53     |      38.65     |         72.22         | 
|  binsMetaBAT.6  |    p__Proteobacteria (UID3214)     |     58.93     |       0.00     |          0.00         | 
|  binsMetaBAT.28 |        k__Bacteria (UID203)        |     54.31     |       5.09     |         31.82         | 
|  binsMetaBAT.34 |        k__Bacteria (UID203)        |     47.18     |       0.00     |          0.00         | 
|  binsMetaBAT.37 |  c__Alphaproteobacteria (UID3337)  |     39.45     |       8.85     |         11.43         |

* Mats - 45 bins

|  Bin Id          |          Marker lineage           |  Completeness |  Contamination  | Strain heterogeneity |
|------------------|-----------------------------------|---------------|-----------------|----------------------|
|  binsMetaBAT.41  |      k__Bacteria (UID2570)        |     93.46     |       2.30      |        33.33         |
|  binsMetaBAT.29  | c__Deltaproteobacteria (UID3216)  |     93.23     |       6.99      |        17.65         |
|  binsMetaBAT.6   |       k__Bacteria (UID203)        |     92.95     |      37.37      |         0.00         |
|  binsMetaBAT.26  |      k__Bacteria (UID3187)        |     90.60     |       3.25      |        33.33         |
|  binsMetaBAT.27  |      k__Bacteria (UID1453)        |     90.18     |       2.63      |         0.00         |
|  binsMetaBAT.2   |      k__Bacteria (UID1453)        |     89.59     |       3.51      |        40.00         |
|  binsMetaBAT.44  |     p__Euryarchaeota (UID3)       |     89.33     |       1.60      |         0.00         |
|  binsMetaBAT.37  |       k__Bacteria (UID203)        |     87.15     |      56.57      |        12.50         |
|  binsMetaBAT.16  | c__Gammaproteobacteria (UID4274)  |     86.57     |      71.30      |        52.22         |
|  binsMetaBAT.14  |       k__Bacteria (UID203)        |     84.48     |       7.42      |        14.29         |
|  binsMetaBAT.23  | c__Gammaproteobacteria (UID4274)  |     83.58     |       3.15      |        11.11         |
|  binsMetaBAT.1   |       k__Bacteria (UID203)        |     75.86     |      56.24      |         1.43         |
|  binsMetaBAT.34  |      k__Bacteria (UID1452)        |     74.72     |       3.28      |        40.00         |
|  binsMetaBAT.12  |      k__Bacteria (UID3187)        |     65.12     |       2.67      |        20.00         |
|  binsMetaBAT.17  |      k__Bacteria (UID2495)        |     62.98     |       1.10      |         0.00         |
|  binsMetaBAT.5   |    c__Spirochaetia (UID2496)      |     55.62     |      15.49      |         2.56         |
|  binsMetaBAT.30  |       k__Bacteria (UID203)        |     49.92     |      19.07      |         6.25         |
|  binsMetaBAT.19  |    c__Spirochaetia (UID2496)      |     45.92     |       0.85      |         0.00         |
|  binsMetaBAT.22  | c__Alphaproteobacteria (UID3337)  |     44.25     |       7.52      |        15.38         |
|  binsMetaBAT.18  |       k__Bacteria (UID203)        |     43.42     |       0.00      |         0.00         |
|  binsMetaBAT.11  |       k__Bacteria (UID203)        |     43.09     |       7.02      |         0.00         |
|  binsMetaBAT.33  |       k__Bacteria (UID203)        |     41.54     |      10.82      |         6.67         |
|  binsMetaBAT.13  |    p__Bacteroidetes (UID2591)     |     41.28     |       0.51      |         0.00         |
|  binsMetaBAT.40  |      k__Bacteria (UID1452)        |     38.18     |       0.91      |         0.00         |

An evaluation of the taxonomic assignation was performed with GTDB-Tk

After installing the program and downloading the required database (Release 95). The data was set as an environmental variable

```
export GTDBTK_DATA_PATH=/mnt/data/yaxal/CC/data/GTDB-Tk_release95/
```

### 16S Identification and extraction for MetaBAT bins

Once with the assembly of the metagenomes and the selection of the bins with MetaBAT, we proceeded to identify and extract the ribosomal genes. The identification was done with Barrnap for all the sequences greater than 2000 bases, as well as only for the identified bins. Barrnap identifies and returns a gff3 file with the 5S, 16S and 23S ribosomal genes, however, we retained only the 16S. With this list, the coordinates of the 16S and the string, the gene sequence was extracted with Bedtools.

Barrnap was executed three times, _(i)_ one without distinguishing between bacteria and archaea, _(ii)_ the second one only for bacteria, and _(iii)_ only for archaea.

Several 16S were identified, however, there were no big differences between the first run of Barrnap, and those dividing by kingdom, bacteria or archaea. When the 16S prediction was done for bins, they were greatly reduced. 


```bash
qsub RunBarrnap.sh
```

|                          | Bordetella  |  Bordetella IDBA | Water Lagunita | Sediment Lagunita | Mats |
| ------------------------ |:-----------:|:----------------:|:--------------:|:-----------------:|:----:|
| All >2000 bases          | 2           | 0                | 39             | 3                 | 35   |
| All >2000 bases bacteria | 2           | 0                | 37             | 3                 | 32   |
| All >2000 bases archaea  | 2           | 0                | 39             | 3                 | 33   |
| All bins MetaBAT         | 0           | NA               | 10             | 0                 | 3    |
| Bins bacteria MetaBAT    | 0           | NA               | 8              | 0                 | 3    |
| Bins archaea MetaBAT     | 0           | NA               | 8              | 0                 | 3    |

* Bordetella

| query name	| Domain	| Phylum	| Class	| Order	| Family	| Genus |
|:-------------:|:-------------:|:-------------:|:-----:|:-----:|:-------------:|:-----:|
| Bordetella_1	| Bacteria	| Proteobacteria	| Betaproteobacteria	| Burkholderiales	| Alcaligenaceae	| Bordetella |
| Bordetella_2	| Archaea	| Thaumarchaeota	| Nitrososphaerales	| Nitrososphaeraceae	| Nitrososphaera	|            |

* Mats

| query name	| Domain	| Phylum	| Class	| Order	| Family	| Genus |
|:-------------:|:-------------:|:-------------:|:-----:|:-----:|:-------------:|:-----:|
| Mats_1	| Bacteria	| Chloroflexi	| unclassified_Chloroflexi |	|	| 	|
| Mats_10	| Archaea	| Thaumarchaeota	| Nitrososphaerales	| Nitrososphaeraceae	| Nitrososphaera	|
| Mats_11	| Bacteria	| unclassified_Bacteria	|	|	|	|
| Mats_12	| Bacteria	| Cyanobacteria/Chloroplast	| Chloroplast	| Chloroplast	| Streptophyta | 	
| Mats_13	| Bacteria	| Cyanobacteria/Chloroplast	| Cyanobacteria	Family VII	|  GpVII	| 
| Mats_14	| Bacteria	| Proteobacteria	| Alphaproteobacteria	| Sphingomonadales	| Sphingomonadaceae	| Sphingomonas | 
| Mats_15	| Bacteria	| Proteobacteria	| Alphaproteobacteria	| Rhizobiales	| Phyllobacteriaceae	| unclassified_Phyllobacteriaceae | 
| Mats_16	| Bacteria	| Proteobacteria	| Gammaproteobacteria	| unclassified_Gammaproteobacteria | 		
| Mats_17	| Fungi 	| Chytridiomycota	| Chytridiomycetes	| Chytridiales	| Chytridiaceae	| Obelidium | 
| Mats_18	| Bacteria	| Verrucomicrobia	| Opitutae	| Puniceicoccales	| Puniceicoccaceae	| unclassified_Puniceicoccaceae | 
| Mats_19	| Bacteria	| Cyanobacteria/Chloroplast	| Chloroplast	| Chloroplast	| Bacillariophyta | 	
| Mats_2	| Bacteria	| Bacteroidetes	Bacteroidia	| Bacteroidales	| unclassified_Bacteroidales	| 
| Mats_20	| Fungi 	| Basidiomycota	Agaricomycetes	| Sebacinales	| Sebacinaceae	| Sebacinaceae | 
| Mats_21	| Bacteria	| Cyanobacteria/Chloroplast	| Cyanobacteria	Family VII	| GpVII	| 
| Mats_22	| Bacteria	| Chloroflexi	Anaerolineae	| unclassified_Anaerolineae | 		
| Mats_23	| Archaea	| Thaumarchaeota	| Nitrososphaerales	| Nitrososphaeraceae	| Nitrososphaera | 	
| Mats_24	| Bacteria	| Cyanobacteria/Chloroplast	| Cyanobacteria	| unclassified_Cyanobacteria | 	| 	| 
| Mats_25	| Fungi 	| Chytridiomycota	| Chytridiomycetes	| Chytridiales	| Chytridiaceae	| Obelidium | 
| Mats_26	| Bacteria	| Balneolaeota	| Balneolia	| Balneolales	| Balneolaceae	| unclassified_Balneolaceae | 
| Mats_27	| Bacteria	| Gemmatimonadetes	| unclassified_Gemmatimonadetes	| 	| 	| 
| Mats_28	| Bacteria	| Verrucomicrobia	| Opitutae	| unclassified_Opitutae	| 	| 
| Mats_29	| Fungi 	| Chytridiomycota	| Chytridiomycetes	| Chytridiales	| Chytridiaceae	| Obelidium
| Mats_3	| Bacteria	| Cyanobacteria/Chloroplast	| Chloroplast	| Chloroplast	| Bacillariophyta	| 
| Mats_30	| Archaea	| Thaumarchaeota	| Nitrososphaerales	| Nitrososphaeraceae	| Nitrososphaera	| 
| Mats_31	| Bacteria	| Bacteroidetes	| unclassified_Bacteroidetes	| 	| 	| 
| Mats_32	| Archaea	| Thaumarchaeota	| Nitrososphaerales	| Nitrososphaeraceae	| Nitrososphaera	| 
| Mats_33	| Bacteria	| Verrucomicrobia	| Subdivision3	| Subdivision3_genera_incertae_sedis	| 	| 
| Mats_34	| Bacteria	| Bacteroidetes	Cytophagia	| Cytophagales	| Reichenbachiellaceae	| Marinoscillum | 
| Mats_35	| Fungi 	| Basidiomycota	| Agaricomycetes	| Sebacinales	| Sebacinaceae	| Sebacinaceae | 
| Mats_4	| Bacteria	| Proteobacteria	| Deltaproteobacteria	| Syntrophobacterales	| Syntrophobacteraceae	| Syntrophobacter | 
| Mats_5	| Bacteria	| Firmicutes	| Bacilli	| Bacillales	| Bacillaceae 1	| Metabacillus | 
| Mats_6	| Fungi 	| Chytridiomycota	| Chytridiomycetes	| Chytridiales	| Chytridiaceae	| Obelidium | 
| Mats_7	| Bacteria	| Bacteroidetes	| Saprospiria	| Saprospirales	| unclassified_Saprospirales	| 
| Mats_8	| Fungi 	| Blastocladiomycota	| Blastocladiomycetes	| Blastocladiales	| Coelomomycetaceae	| Coelomomyces | 
| Mats_9	| Bacteria	| Cyanobacteria/Chloroplast	| Cyanobacteria	| Family IX	| GpIX	| 

* Sediment Lagunita

| query name	| Domain	| Phylum	| Class	| Order	| Family	| Genus |
|:-------------:|:-------------:|:-------------:|:-----:|:-----:|:-------------:|:-----:|
| SedimentLagunita_1	| Bacteria	| Actinobacteria	| Actinobacteria	| Bifidobacteriales	| Bifidobacteriaceae	| Bifidobacterium |
| SedimentLagunita_2	| Bacteria	| Actinobacteria	| Coriobacteriia	| Coriobacteriales	| Coriobacteriaceae	| Collinsella |
| SedimentLagunita_3	| Bacteria	| Bacteroidetes 	| Saprospiria   	| Saprospirales   	| unclassified_Saprospirales	| |

* Water Lagunita

| query name	| Domain	| Phylum	| Class	| Order	| Family	| Genus |
|:-------------:|:-------------:|:-------------:|:-----:|:-----:|:-------------:|:-----:|
| WaterLagunita_1	| Bacteria	|	Cyanobacteria/Chloroplast	| Chloroplast	| Chloroplast	| Chlorophyta	| 	| 
| WaterLagunita_10	| Bacteria	|	Proteobacteria	| Alphaproteobacteria	| Rhodospirillales	| Acetobacteraceae	| unclassified_Acetobacteraceae	| 
| WaterLagunita_11	| Bacteria	|	Proteobacteria	| Betaproteobacteria	| Burkholderiales	| Alcaligenaceae	| unclassified_Alcaligenaceae	| 
| WaterLagunita_12	| Archaea	|	Thaumarchaeota	| Nitrososphaerales	| Nitrososphaeraceae	| Nitrososphaera	| 	| 
| WaterLagunita_13	| Archaea	|	Thaumarchaeota	| Nitrososphaerales	| Nitrososphaeraceae	| Nitrososphaera	| 	| 
| WaterLagunita_14	| Bacteria	|	Proteobacteria	| Alphaproteobacteria	| Rhodobacterales	| Rhodobacteraceae	| Rubribacterium	| 
| WaterLagunita_15	| Bacteria	|	Proteobacteria	| Alphaproteobacteria	| Caulobacterales	| Hyphomonadaceae	| unclassified_Hyphomonadaceae	| 
| WaterLagunita_16	| Bacteria	|	Proteobacteria	| Alphaproteobacteria	| Kiloniellales	| Kiloniellaceae	| Thalassocola	| 
| WaterLagunita_17	| Bacteria	|	Bacteroidetes	| Saprospiria	| Saprospirales	| Lewinellaceae	| unclassified_Lewinellaceae	| 
| WaterLagunita_18	| Bacteria	|	Proteobacteria	| Deltaproteobacteria	| unclassified_Deltaproteobacteria	| 	| 	| 
| WaterLagunita_19	| Fungi	|	Basidiomycota	| Agaricomycetes	| Sebacinales	| Sebacinaceae	| Sebacinaceae	| 
| WaterLagunita_2	| Bacteria	|	Cyanobacteria/Chloroplast	| Chloroplast	| Chloroplast	| Chlorophyta	| 	| 
| WaterLagunita_20	| Bacteria	|	Bacteroidetes	| unclassified_Bacteroidetes	| 	| 	| 	| 
| WaterLagunita_21	| Fungi	|	Basidiomycota	| Agaricomycetes	| Sebacinales	| Sebacinaceae	| Sebacinaceae	| 
| WaterLagunita_22	| Archaea	|	Thaumarchaeota	| Nitrososphaerales	| Nitrososphaeraceae	| Nitrososphaera	| 	| 
| WaterLagunita_23	| Archaea	|	Thaumarchaeota	| Nitrososphaerales	| Nitrososphaeraceae	| Nitrososphaera	| 	| 
| WaterLagunita_24	| Fungi	|	Chytridiomycota	| Chytridiomycetes	| Chytridiales	| Chytridiaceae	| Obelidium	| 
| WaterLagunita_25	| Bacteria	|	Verrucomicrobia	| Verrucomicrobiae	| Verrucomicrobiales	| Verrucomicrobiaceae	| Luteolibacter	| 
| WaterLagunita_26	| Bacteria	|	Proteobacteria	| Deltaproteobacteria	| unclassified_Deltaproteobacteria	| 	| 	| 
| WaterLagunita_27	| Bacteria	|	Proteobacteria	| Alphaproteobacteria	| Caulobacterales	| Hyphomonadaceae	| Ponticaulis	| 
| WaterLagunita_28	| Bacteria	|	Cyanobacteria/Chloroplast	| Chloroplast	| Chloroplast	| unclassified_Chloroplast	| 	| 
| WaterLagunita_29	| Bacteria	|	Proteobacteria	| Oligoflexia	| Bacteriovoracales	| Bacteriovoracaceae	| Peredibacter	| 
| WaterLagunita_3	| Bacteria	|	SR1	| SR1_genera_incertae_sedis	| 	| 	| 	| 
| WaterLagunita_30	| Bacteria	|	Bacteroidetes	| unclassified_Bacteroidetes	| 	| 	| 	| 
| WaterLagunita_31	| Bacteria	|	Bacteroidetes	| Flavobacteriia	| Flavobacteriales	| Weeksellaceae	| Cloacibacterium	| 
| WaterLagunita_32	| Bacteria	|	Proteobacteria	| Alphaproteobacteria	| Kiloniellales	| Kiloniellaceae	| Thalassocola	| 
| WaterLagunita_33	| Bacteria	|	Firmicutes	| Bacilli	| Lactobacillales	| Streptococcaceae	| Streptococcus	| 
| WaterLagunita_34	| Bacteria	|	Bacteroidetes	| Saprospiria	| Saprospirales	| unclassified_Saprospirales	| 	| 
| WaterLagunita_35	| Fungi	|	Blastocladiomycota	| Blastocladiomycetes	| Blastocladiales	| Coelomomycetaceae	| Coelomomyces	| 
| WaterLagunita_36	| Bacteria	|	unclassified_Bacteria	| 	| 	| 	| 	| 
| WaterLagunita_37	| Bacteria	|	Proteobacteria	| unclassified_Proteobacteria	| 	| 	| 	| 
| WaterLagunita_38	| Bacteria	|	Proteobacteria	| Alphaproteobacteria	| Caulobacterales	| Hyphomonadaceae	| Hyphomonas	| 
| WaterLagunita_39	| Bacteria	|	Bacteroidetes	| unclassified_Bacteroidetes	| 	| 	| 	| 
| WaterLagunita_4	| Archaea	|	Thaumarchaeota	| Nitrososphaerales	| Nitrososphaeraceae	| Nitrososphaera	| 	| 
| WaterLagunita_5	| Bacteria	|	Firmicutes	| Bacilli	| Bacillales	| Staphylococcaceae	| Staphylococcus	| 
| WaterLagunita_6	| Bacteria	|	Proteobacteria	| Alphaproteobacteria	| unclassified_Alphaproteobacteria	| 	| 	| 
| WaterLagunita_7	| Bacteria	|	Proteobacteria	| Betaproteobacteria	| Burkholderiales	| unclassified_Burkholderiales	| 	| 
| WaterLagunita_8	| Archaea	|	Thaumarchaeota	| Nitrososphaerales	| Nitrososphaeraceae	| Nitrososphaera	| 	| 
| WaterLagunita_9	| Fungi	|	Chytridiomycota	| Chytridiomycetes	| Chytridiales	| Chytridiaceae	| Obelidium	| 

### CONCOCT binning

A second approach for binning inference was implemented with CONCOCT (https://concoct.readthedocs.io/en/latest/) In contrast to MetaBAT which does the binning based on tetranucleotide (4-mers) frequency in the samples, CONCOCT does the binning by using nucleotide composition, k-mer frequencies and coverage data.

As it was run with conda, the Concoct conda environment has to be activated and deactivated during the pipeline. 

To activate this environment, use

```bash
conda activate concoct_env
```

To deactivate an active environment, use
```bash
conda deactivate
```

The CONCOCT binning and  evaluation with CheckM were run in the same script.

```bash
qsub RunConcoct.sh*
```

The evaluation of the binning with CONCOCT, returned a total of 192 bins distributed among the four samples, with the mats containing the greatest diversity with 82 bins, and the Bordetella sample containing 4 bins. We only display the values for those bins with a completeness higher than 40%

* Bordetella - 4 bins

* Bordetella_IDBA - 4 bins

* Lagunita Sediment - 27 bins

|  Bin Id  |         Marker lineage          | Completeness | Contamination |  Strain heterogeneity  |
|----------|---------------------------------|--------------|---------------|------------------------|
|  0       | o__Rhodospirillales (UID3754)   |    56.05     |      1.69     |         66.67          |
|  14      |  f__Lachnospiraceae (UID1286)   |    44.13     |      0.65     |         0.00           |

* Lagunita Water - 79 bins

|  Bin Id  |          Marker lineage           | Completeness | Contamination | Strain heterogeneity  |
|----------|-----------------------------------|--------------|---------------|-----------------------|
|  45      |       k__Bacteria (UID203)        |    100.00    |     136.93    |        57.30          |
|  38      |       k__Bacteria (UID203)        |    100.00    |     79.15     |         0.00          |
|  24      |    p__Bacteroidetes (UID2591)     |    99.49     |      2.02     |         0.00          |
|  48      |  f__Rhodobacteraceae (UID3340)    |    98.52     |      1.03     |         0.00          |
|  19      |    p__Bacteroidetes (UID2591)     |    98.51     |      0.99     |        50.00          |
|  36      |       k__Bacteria (UID203)        |    98.28     |     94.83     |         0.00          |
|  11      |    p__Bacteroidetes (UID2605)     |    97.14     |      1.67     |        50.00          |
|  47      | c__Alphaproteobacteria (UID3422)  |    96.71     |      2.38     |        72.73          |
|  15      |      k__Bacteria (UID2982)        |    95.75     |      1.40     |        33.33          |
|  18      |   p__Proteobacteria (UID3887)     |    94.26     |      2.34     |        100.00         |
|  67      |      k__Bacteria (UID3187)        |    93.59     |      1.79     |         0.00          |
|  10      |      k__Bacteria (UID3187)        |    93.30     |      3.57     |        75.00          |
|  60      |       k__Bacteria (UID203)        |    92.84     |     89.09     |        71.70          |
|  41      |    p__Bacteroidetes (UID2591)     |    92.55     |      1.26     |         0.00          |
|  53      | c__Alphaproteobacteria (UID3305)  |    91.52     |      0.43     |        100.00         |
|  3       |  f__Rhodobacteraceae (UID3340)    |    89.45     |      1.99     |        30.00          |
|  56      | c__Alphaproteobacteria (UID3422)  |    89.08     |      8.12     |        16.67          |
|  74      |   p__Proteobacteria (UID3887)     |    82.25     |      2.44     |        42.86          |
|  23      |   p__Proteobacteria (UID3214)     |    63.85     |      0.00     |         0.00          |
|  33      |       k__Bacteria (UID203)        |    58.70     |      7.13     |        26.92          |
|  16      |     o__Rhizobiales (UID3447)      |    54.57     |      0.40     |         0.00          |
|  69      | c__Alphaproteobacteria (UID3337)  |    43.43     |     10.33     |        30.23          |
 
* Mats - 82 bins

|  Bin Id  |          Marker lineage           |  Completeness  | Contamination  | Strain heterogeneity  |
|----------|-----------------------------------|----------------|----------------|-----------------------|
|  0       |       k__Bacteria (UID203)        |     96.55      |     88.64      |        72.63          |
|  71      | c__Deltaproteobacteria (UID3216)  |     94.72      |     68.03      |         3.19          |
|  19      |      k__Bacteria (UID2570)        |     93.46      |      3.93      |        22.22          |
|  21      |      k__Bacteria (UID2569)        |     90.93      |      1.88      |         0.00          |
|  41      |     p__Euryarchaeota (UID3)       |     89.33      |      1.60      |         0.00          |
|  64      |      k__Bacteria (UID3187)        |     89.32      |      3.04      |        40.00          |
|  14      |      k__Bacteria (UID1453)        |     88.71      |      2.63      |        25.00          |
|  48      |      k__Bacteria (UID1453)        |     88.42      |      4.39      |         0.00          |
|  79      |       k__Bacteria (UID203)        |     87.30      |     64.48      |        10.67          |
|  13      | c__Gammaproteobacteria (UID4274)  |     82.22      |      3.48      |         7.41          |
|  39      |  o__Rhodospirillales (UID3754)    |     78.99      |      6.30      |        20.00          |
|  4       |      k__Bacteria (UID1452)        |     75.17      |      3.28      |        40.00          |
|  16      |      k__Bacteria (UID3187)        |     62.64      |      0.96      |        33.33          |
|  40      |      k__Bacteria (UID2495)        |     61.71      |      3.75      |         0.00          |
|  24      |       k__Bacteria (UID203)        |     60.58      |     21.25      |         5.00          |
|  37      |    p__Bacteroidetes (UID2605)     |     53.32      |     14.50      |         5.88          |
|  25      | c__Alphaproteobacteria (UID3337)  |     47.31      |      8.14      |        18.18          |
|  81      |      k__Bacteria (UID2982)        |     46.57      |      2.17      |         0.00          |
|  2       |    p__Bacteroidetes (UID2591)     |     45.82      |      2.36      |        14.29          |
|  29      |    c__Spirochaetia (UID2496)      |     45.65      |      1.67      |         0.00          |
|  7       |      k__Bacteria (UID2570)        |     42.60      |     11.78      |         1.79          |
|  67      |      k__Bacteria (UID1452)        |     42.43      |      6.11      |         0.00          |
|  33      |       k__Bacteria (UID203)        |     41.74      |     35.39      |         1.17          |
|  22      |    p__Bacteroidetes (UID2591)     |     41.63      |      0.51      |         0.00          |
|  27      |       k__Bacteria (UID203)        |     40.83      |      7.84      |         9.09          |
|  1       |      k__Bacteria (UID3187)        |     40.62      |      4.38      |         0.00          |

### DAS_Tool

Using the Megahit metagenome assembly, and the two binning approaches with MetaBat and Concoct, the integration of both results was done with DAS Tool to recover and optimized and non-redundant set of bins.

A fragment of the script to prepare the list and run DAS Tool is shown below.

```bash
qsub RunDasTool.sh
```

**Bins selected:**

* Bordetella - 1 bin selected from 2

* Mats - 17 bins selected from 96

* Sediment - 2 bins selected from 13

* Water - 21 bins selected from 72

To keep a structure and a meaningful nomenclature, the integrated bins by DAS Tool, were renamed to give information regarding the sample origin and bin number. For the case of Water bins, as we have bins assigned to the same number (metabat.016 and concoct.016), the bin concoct.016.fa was renamed to Water.Bin_73.fa to avoid overwriting the files.

```bash
bash RenameDasToolBins.sh
```

After renaming the integrated bins, they were evaluated with CheckM, and the 16S sequences were predicted and extracted with Barrnp.

```bash
qsub RunCheckMBarrnapDASToolBins.sh
```

The DAS_Tool integrated bins evaluation resulted in overall better statistics than individual binning approaches.

* Bordetella - 1 bin

|     Bin Id       |         Marker lineage        |  Completeness |  Contamination  | Strain heterogeneity  |
|:----------------:|:-----------------------------:|:-------------:|:---------------:|:---------------------:|
|  Churince.Bin_0  |  o__Burkholderiales (UID4000) |     26.98     |       0.51      |        50.00          |

* Mats - 17 bins

|   Bin Id    |         Marker lineage        |  Completeness |  Contamination  | Strain heterogeneity  |
|:-----------:|:-----------------------------:|:-------------:|:---------------:|:---------------------:|
|  Mats.Bin_21 |  c__Deltaproteobacteria (UID3216)  | 93.23   |         6.99    |          17.65        |  
|  Mats.Bin_14 |       k__Bacteria (UID2569)        | 90.93   |         1.88    |           0.00        |  
|  Mats.Bin_19 |       k__Bacteria (UID1453)        | 90.18   |         2.63    |           0.00        |  
|  Mats.Bin_22 |       k__Bacteria (UID1453)        | 89.59   |         3.51    |          40.00        |  
|  Mats.Bin_38 |      p__Euryarchaeota (UID3)       | 89.33   |         1.60    |           0.00        |  
|  Mats.Bin_61 |       k__Bacteria (UID3187)        | 89.32   |         3.04    |          40.00        |  
|  Mats.Bin_7  |  c__Gammaproteobacteria (UID4274)  | 86.57   |        71.30    |          52.22        |  
|  Mats.Bin_5  |  c__Gammaproteobacteria (UID4274)  | 82.22   |         3.48    |           7.41        |  
|  Mats.Bin_33 |   o__Rhodospirillales (UID3754)    | 78.99   |         6.30    |          20.00        |  
|  Mats.Bin_45 |       k__Bacteria (UID1452)        | 75.17   |         3.28    |          40.00        |  
|  Mats.Bin_8  |       k__Bacteria (UID2495)        | 62.98   |         1.10    |           0.00        |  
|  Mats.Bin_41 |     c__Spirochaetia (UID2496)      | 55.62   |        15.49    |           2.56        |  
|  Mats.Bin_11 |  c__Alphaproteobacteria (UID3337)  | 54.51   |        15.98    |           1.52        |  
|  Mats.Bin_80 |       k__Bacteria (UID2982)        | 45.89   |         1.35    |           0.00        |  
|  Mats.Bin_23 |     p__Bacteroidetes (UID2591)     | 45.82   |         2.36    |          14.29        |  
|  Mats.Bin_4  |       k__Bacteria (UID2328)        | 33.33   |         0.12    |           0.00        |  
|  Mats.Bin_1  |  c__Gammaproteobacteria (UID4274)  |  8.48   |         0.17    |          100.00       | 

* Sediment - 2 bins

|   Bin Id    |         Marker lineage        |  Completeness |  Contamination  | Strain heterogeneity  |
|:-----------:|:-----------------------------:|:-------------:|:---------------:|:---------------------:|
|  Sediment.Bin_1 |  o__Rhodospirillales (UID3754)  | 56.05   |         1.69    |          66.67        |  
|  Sediment.Bin_6 |   f__Lachnospiraceae (UID1286)  | 44.13   |         0.65    |           0.00        |  

* Water - 21 bins

|   Bin Id    |         Marker lineage        |  Completeness |  Contamination  | Strain heterogeneity  |
|:-----------:|:-----------------------------:|:-------------:|:---------------:|:---------------------:|
|  Water.Bin_16  |    p__Bacteroidetes (UID2591)    |  99.49  |          2.02   |            0.00       |   
|  Water.Bin_43  |  f__Rhodobacteraceae (UID3340)   |  98.52  |          1.03   |            0.00       |   
|  Water.Bin_37  |    p__Bacteroidetes (UID2591)    |  98.02  |          0.99   |           50.00       |   
|  Water.Bin_32  |      k__Bacteria (UID2570)       |  96.45  |          0.00   |            0.00       |   
|  Water.Bin_22  |      k__Bacteria (UID2982)       |  95.75  |          1.40   |           33.33       |   
|  Water.Bin_15  | c__Alphaproteobacteria (UID3422) |  95.55  |          1.29   |           85.71       |   
|  Water.Bin_4   |    p__Bacteroidetes (UID2605)    |  95.24  |          1.67   |           50.00       |   
|  Water.Bin_11  |    p__Bacteroidetes (UID2605)    |  94.29  |          0.00   |            0.00       |   
|  Water.Bin_10  |   p__Proteobacteria (UID3887)    |  94.26  |          2.34   |           100.00      |   
|  Water.Bin_64  |      k__Bacteria (UID3187)       |  93.59  |          1.79   |            0.00       |   
|  Water.Bin_1   |      k__Bacteria (UID2565)       |  93.55  |          0.00   |            0.00       |   
|  Water.Bin_2   |      k__Bacteria (UID3187)       |  93.30  |          3.57   |           75.00       |   
|  Water.Bin_17  |  f__Rhodobacteraceae (UID3340)   |  93.02  |          5.02   |           25.00       |   
|  Water.Bin_36  |    p__Bacteroidetes (UID2591)    |  92.55  |          1.26   |            0.00       |   
|  Water.Bin_49  | c__Alphaproteobacteria (UID3305) |  91.52  |          0.43   |           100.00      |   
|  Water.Bin_52  | c__Alphaproteobacteria (UID3422) |  89.08  |          7.80   |           18.18       |   
|  Water.Bin_72  |   p__Proteobacteria (UID3887)    |  82.25  |          2.44   |           42.86       |   
|  Water.Bin_33  |       k__Bacteria (UID203)       |  75.71  |          0.00   |            0.00       |   
|  Water.Bin_73  |   p__Proteobacteria (UID3214)    |  63.85  |          0.00   |            0.00       |   
|  Water.Bin_8   |     o__Rhizobiales (UID3447)     |  54.57  |          0.40   |            0.00       |   
|  Water.Bin_9   |   o__Burkholderiales (UID4000)   |  34.18  |          1.25   |           25.00       |   

### GTDB-Tk

Taking once again the integrated bins from DAS_Tool, a taxonomic assignation was performed with GTDB-Tk. GTDB-Tk is a software for assigning taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy GTDB, and is designed to work with metagenome-assembled genomes (MAGs).

**Technical note.** While running GTDB-Tk, the process was killed several times on the server. Originally it was assumed to be caused by a high memory requirements, but that was not the case. GTDB-Tk uses pplacer in one of its steps, and the default behaviour is to request as much memory to each thread, as the total memory assigned, because of this, the process dies with a low memory error. To overcome this, two options need to be passed, _--pplacer_cpus 1_ and _--scratch_dir /tmp/dir/_.

```bash
qsub RunGTDB_Tk.sh
```

Classification results.

* Bordetella / Churince - 1 bin

| user_genome	 | Domain | Phylum | Class | Order | Family | Genus | Species |
|:--------------:|:------:|:------:|:-----:|:-----:|:------:|:-----:|:-------:|
| Churince.Bin_0 | Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | Burkholderiaceae | Bordetella | Bordetella hinzii |

* Mats - 17 bins

| user_genome	 | Domain | Phylum | Class | Order | Family | Genus | Species |
|:--------------:|:------:|:------:|:-----:|:-----:|:------:|:-----:|:-------:|
| Mats.Bin_38 | 	Archaea | Thermoplasmatota | E2 | DHVEG-1 | DHVEG-1 | EX4572-165 | 
| Mats.Bin_1	 | Bacteria | Proteobacteria | Gammaproteobacteria | Ectothiorhodospirales | Ectothiorhodospiraceae | Ectothiorhodospira | 
| Mats.Bin_11 | 	Bacteria | Proteobacteria | Alphaproteobacteria | Rhodobacterales | Rhodobacteraceae |  | 
| Mats.Bin_14 | 	Bacteria | Bacteroidota | Bacteroidia | Bacteroidales | UBA12077 | UBA12077 | 
| Mats.Bin_19 | 	Bacteria | Cyanobacteria | Cyanobacteriia | PCC-7336 | JA-3-3Ab |  | 
| Mats.Bin_21 | 	Bacteria | Desulfobacterota | Desulfobacteria | Desulfobacterales | UBA2174 |  | 
| Mats.Bin_22 | 	Bacteria | Cyanobacteria | Cyanobacteriia | PCC-7336 | JA-3-3Ab |  | 
| Mats.Bin_23 | 	Bacteria | Bacteroidota | Bacteroidia | Chitinophagales | Saprospiraceae |  | 
| Mats.Bin_33 | 	Bacteria | Proteobacteria | Alphaproteobacteria | Geminicoccales |  |  | 
| Mats.Bin_4	 | Bacteria | Firmicutes | Bacilli | Izemoplasmatales | Izemoplasmataceae |  | 
| Mats.Bin_41 | 	Bacteria | Acidobacteriota | Mor1 | J045 |  |  | 
| Mats.Bin_45 | 	Bacteria | Chloroflexota | Anaerolineae | 4572-78 | J111 |  | 
| Mats.Bin_5	 | Bacteria | Proteobacteria | Gammaproteobacteria | Ectothiorhodospirales |  |  | 
| Mats.Bin_61 | 	Bacteria | Acidobacteriota | Acidobacteriae | Bryobacterales |  |  | 
| Mats.Bin_7	 | Bacteria | Proteobacteria | Gammaproteobacteria | Ectothiorhodospirales | Ectothiorhodospiraceae | Ectothiorhodospira | 
| Mats.Bin_8	 | Bacteria | Zixibacteria | MSB-5A5 | GN15 | FEB-12 |  | 
| Mats.Bin_80 | 	Bacteria | Verrucomicrobiota | Verrucomicrobiae | Opitutales | Puniceicoccaceae |  | 

* Sediment - 2 bins

| user_genome	 | Domain | Phylum | Class | Order | Family | Genus | Species |
|:--------------:|:------:|:------:|:-----:|:-----:|:------:|:-----:|:-------:|
| Sediment.Bin_1 |	Bacteria | Proteobacteria | Alphaproteobacteria | Geminicoccales |  |  | 
| Sediment.Bin_6 |	Bacteria | Firmicutes_A | Clostridia | Lachnospirales | Lachnospiraceae | Agathobacter | Agathobacter rectalis |
 
* Water - 21 bins

| user_genome	 | Domain | Phylum | Class | Order | Family | Genus | Species |
|:--------------:|:------:|:------:|:-----:|:-----:|:------:|:-----:|:-------:| 
| Water.Bin_1 | Bacteria | Proteobacteria | Alphaproteobacteria | UBA7879 | UBA5542 | UBA5542 | 
| Water.Bin_10 | Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | Burkholderiaceae | UBA2463 | 
| Water.Bin_11 | Bacteria | Bacteroidota | Bacteroidia | NS11-12g | UBA8524 |  | 
| Water.Bin_15 | Bacteria | Proteobacteria | Alphaproteobacteria | Caulobacterales | Maricaulaceae | Oceanicaulis | 
| Water.Bin_16 | Bacteria | Bacteroidota | Bacteroidia | Chitinophagales | Saprospiraceae |  | 
| Water.Bin_17 | Bacteria | Proteobacteria | Alphaproteobacteria | Rhodobacterales | Rhodobacteraceae | Defluviimonas_B | 
| Water.Bin_2 | Bacteria | Bdellovibrionota | Bdellovibrionia | Bdellovibrionales | UBA1609 |  | 
| Water.Bin_22 | Bacteria | Verrucomicrobiota | Verrucomicrobiae | Verrucomicrobiales | Akkermansiaceae | UBA1315 | 
| Water.Bin_32 | Bacteria | Bacteroidota | Bacteroidia | UBA7662 |  |  | 
| Water.Bin_33 | Bacteria | Patescibacteria | Gracilibacteria | Absconditabacterales | X112 | SKIU01 | 
| Water.Bin_36 | Bacteria | Bacteroidota | Bacteroidia | Chitinophagales | Saprospiraceae |  | 
| Water.Bin_37 | Bacteria | Bacteroidota | Bacteroidia | Chitinophagales | Saprospiraceae | M3007 | 
| Water.Bin_4 | Bacteria | Bacteroidota | Bacteroidia | Flavobacteriales | UBA10329 |  | 
| Water.Bin_43 | Bacteria | Proteobacteria | Alphaproteobacteria | Rhodobacterales | Rhodobacteraceae | UBA9554 | 
| Water.Bin_49 | Bacteria | Proteobacteria | Alphaproteobacteria | Micavibrionales_C | TMED2 | TMED2 | 
| Water.Bin_52 | Bacteria | Proteobacteria | Alphaproteobacteria | Caulobacterales | Hyphomonadaceae |  | 
| Water.Bin_64 | Bacteria | Bdellovibrionota | Bdellovibrionia | Bdellovibrionales | UBA1609 |  | 
| Water.Bin_72 | Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | Burkholderiaceae | UBA2463 | 
| Water.Bin_73 | Bacteria | Proteobacteria | Alphaproteobacteria | UBA1280 |  |  | 
| Water.Bin_8 | Bacteria | Proteobacteria | Alphaproteobacteria | Rhizobiales | UBA4765 |  | 
| Water.Bin_9 | Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | Burkholderiaceae | Polaromonas | 

### Cleaning

The cleaning of the bins was accomplished using mmgenome2. This step allows to identify contigs that can be assigned as contamination inside the identified bins. As we previously identified candidate contaminated bins win CheckM, we will focus only in those with levels of contamination  greater than 5% and Completeness greater than 50%.

Briefly, the idea of this cleaning, is to use the coverage of the bins as key point to identify any possible contamination. Different organisms should have different depth values, which will allow us to identify those outlayers. mmgenome2 can load multiple coverage files for each bin, however, it can also work with only one. **Technical note.** In our case, as we have all the fastq file concatenated, at the moment of plotting the distribution and selecting the contigs of interest, we have to use the same coverage file for X and Y axis which gave an error. To overcome this, two options were found, the first one, was to copy the coverage file and rename it, that way, at the moment of loading the axis, we could work. The second option was to use the GC% (gc) as an axis and select the contigs. The second option I believe is the most appropriate. 

A "_1" was added to the name of the cleaned bins, and the "old" bin was renamed from \*.fa to \*.fna. 
 
#### Data files preparation.

As we need to have the coverage for each bin, they were mapped to the clean fastq files

```bash
qsub MappingCoverageIntegratedBins.sh
```

The curation of each bin was performed manually in RStudio and the package `mmgenome2`. Here we present the example for Mats.Bin.21

```R
library(mmgenome2)
setwd("/media/yaxal/Emiru/Genomes/Metagenomes/ChurinceCC/Lagunita/DASTool/Mats")

# Load the bin to clean
mm <- mmload(
  assembly = "Mats_DASTool_bins/Mats.Bin_21.fa",
  coverage = "coverage/",
  verbose = TRUE,
  kmer_pca = FALSE,
  kmer_BH_tSNE = FALSE
)

# Plot the coverage files, if there were multiple files, it would need a vector with all of them
mmplot(mm,
	x = "cov_Mats.Bin_21", 
	y = "gc",
	min_length = 2500,
	color_by = "gc",
	x_scale = "log10",
	x_limits = c(7, 30),
	y_scale = "log10",
	y_limits = c(47, 80),
	locator = TRUE
	)
	
# Manually select the area of interest
shiny_selection <- data.frame(cov_Mats.Bin_21 = c(23.133, 20.177, 17.51, 13.594, 11.932, 10.186, 11.415, 13.525),
           gc = c(64.708, 71.735, 71.917, 69.292, 65.423, 57.198, 55.813, 55.297))
           
mmsub <- mmextract(mm, selection = shiny_selection)

mmplot(mmsub,
	x = "cov_Mats.Bin_21", 
	y = "gc",
	min_length = 2500,
	color_by = "gc",
	x_scale = "log10",
	x_limits = c(10, 25),
	y_scale = "log10",
	y_limits = c(53, 72),
	locator = TRUE
	)

shiny_selection <- data.frame(cov_Mats.Bin_21 = c(22.768, 13.913, 10.833, 12.238, 18.316, 22.551),
           gc = c(69.96, 69.994, 61.974, 55.33, 59.289, 65.002))

# Extract the selected contigs           
mmsubset <- mmextract(mmsub, selection = shiny_selection)

# Export to a fasta file
mmexport(mmsubset, assembly = assembly, file = "Mats_DASTool_bins/Mats.Bin_21_1.fa")
```

In this case, 478 scaffolds (or 97.51% of the scaffolds in mm, weighted by length) remain after 22 of 500 scaffolds have been filtered.

### Phylosift

A phylogenetic and taxonomic analysis of the samples was performed with Phylosift. The installation of PhyloSift was done from GitHub, as the tarball is no longer available.

```bash
cd /home/yaxal/Downloads
git clone git://github.com/gjospin/PhyloSift.git
# Install CPAN dependencies as root
sudo su
cpan
cpan[1]> Bio::Phylo
cpan[2]> LWP::Simple
# Install other dependencies as needed

# Database.
# The PhyloSift database was provided by Valerie de Anda
cd /home/yaxal
mkdir share
cd share 
mv ../phylosift.tar.gz /home/yaxal/share
tar -xvzf phylosift.tar.gz
```

After installing and configuring PhyloSift, it was executed

```bash

cd 
ln -s ../DASTool/Mats/Mats_DASTool_bins/*.fa .
ln -s ../DASTool/Sediment/Sediment_DASTool_bins/*.fa 
ln -s ../DASTool/Water/Water_DASTool_bins/*.fa .

find . -maxdepth 1 -name "*.fa" -exec /home/yaxal/Downloads/PhyloSift/phylosift search --isolate --besthit {} \; 

find . -maxdepth 1 -name "*.fa" -exec /home/yaxal/Downloads/PhyloSift/phylosift align --isolate --besthit --threads 3 {} \;

# Extract the alignments for all bins and samples into a single file
find . -type f -regex '.*alignDir/concat.updated.1.fasta' -exec cat {} \; | sed -r 's/\.1\..*//'  > AllSamples_alignment.fna

# Extract the alignments for each sample
find PS_temp/Churince.Bin_0.fa -type f -regex '.*alignDir/concat.updated.1.fasta' -exec cat {} \; | sed -r 's/\.1\..*//'  > Churince_alignment.fa
find PS_temp/Mats.Bin_* -type f -regex '.*alignDir/concat.updated.1.fasta' -exec cat {} \; | sed -r 's/\.1\..*//'  > Mats_alignment.fa
find PS_temp/Sediment.Bin_* -type f -regex '.*alignDir/concat.updated.1.fasta' -exec cat {} \; | sed -r 's/\.1\..*//'  > Sediment_alignment.fa
find PS_temp/Water.Bin_* -type f -regex '.*alignDir/concat.updated.1.fasta' -exec cat {} \; | sed -r 's/\.1\..*//'  > Water_alignment.fa

```

Valerie de Anda provided an aligment containing a selection of different bacteria an archaea used in a previous study. The PhyloSift alignment from water, sediment and mats were incorporated into the provided alignment. Subsequently, the column of the sequences with gaps in at least 50% of the taxa were removed with trimAl v1.4.rev22 build[2015-05-21].

```bash
muscle -profile -in1 2020_10_04-SIP_hydrocarbon-alignment-checked_alignment_with_manually_add_reference_modify_round9-SF_Bacteria_alignment-realigned.fasta -in2 AllSaplesAlignment.fna -out 2021_07_31_Bacteria_alignment.fna

trimal -in 2021_07_31_Bacteria_alignment.fna -out 2021_07_31_Bacteria_alignment_gaps.phy -phylip -gt 0.5
```

A total of 25,669 characters, 6,305 were retained. The trimmed alignment was used for a phylogenetic reconstruction. First, the best fitting model was selected with SMS v1.8.4. Using the best fitted model, a phylogenetic tree was computed with IQ-Tree version 2.1.4-beta COVID-edition

```bash
cd /mnt/data/yaxal/CC
qsub SMS.sh
qsub IQTree.sh
```

