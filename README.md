![example_output](/images/logov4.png)

### Traposome Sequencing Pipeline for Long Read Technologies


- [Installation](#installation)
- [Usage](#usage)
- [Input Data Analysis](#input-data-analysis)
  - [Read Length Histograms](#read-length-histograms)
- [Alignment](#alignment)
  - [Map](#map)
- [Analysis](#analysis)
  - [Insertion Finder](#insertion-finder)
  - [Manually curated Transposons](#manually-curated-transposons)
  - [Annotation of Insertions](#annotation-of-insertions)
  - [Examining the novel inserts](#examining-the-novel-inserts)
  - [Stats of Transposons](#stats-of-transposons)
  - [Heatmap of Transposons](#heatmap-of-transposons)
  - [Insertion Coordinates](#insertion-coordinates)
  - [DR Finder](#dr-finder)
  - [DR Sequence Logos](#dr-sequence-logos)
  - [Kde Plots](#kde-plots)
  - [Possible Excisions](#possible-excisions)
  - [Mutation and Variant Analysis](#mutation-and-variant-analysis)
  - [Map Distribution](#map-distribution)
- [FAQ](#faq)
  - [How Can I Cite trapo-seq?](#how-can-i-cite-trapo-seq?)
- [License](#license)
- [Flowchart](#flowchart)


**trapo-seq** is a bacterial traposome sequencing pipeline for long read technologies. The biggest advantage of the long read technology for traposome analysis is that it can resolve mobile elements (transposons and insertion sequences are repetitive sequences) as a whole, even just one read can reveal the complete sequence of the mobile elements that have trapped by the plasmid (or other genetic elements). Also, the assembly process is not a must unlike short read technologies. This also helps to obtain high throughput traposome data to compare effects of different conditions on transposition. 


**trapo-seq** was initially developed for trap plasmids, but it can be used for other genetic elements if there are insertions of mobile elements on them.

# Installation

## Requirements
**For all modules:**
- [Linux]()
- [Python]() (>v3.5) (For all modules)

**Libraries:**
- [matplotlib]() (read_histogram, heatmap, dr_logo, in_del_plot, kde_mobile, map_dist)
- [biopython]() (filter, blast_annot_batch, is_stat, heatmap, dr_finder, dr_logo)
- [pandas]() (is_stat, heatmap, dr_finder, dr_logo, in_del_plot)
- [seaborn]() (heatmap, dr_logo, kde_mobile)
- [numpy]() (heatmap, dr_logo)
- [pysam]() (insert_finder_batch, in_del_plot, kde_mobile, map_dist)
- [logomaker]() (dr_logo)
### **Option 1:** Using pip

```
pip install matplotlib biopython pandas seaborn numpy pysam logomaker
```
### **Option 2:** Using conda **(recommended)**
```
conda create -n trapo-seq-env python=3.8 -y
conda activate trapo-seq-env
pip install matplotlib biopython pandas seaborn numpy pysam logomaker
```

**External programs**: for installation please visit the links below, each can have different procedures:
- [minimap2](https://github.com/lh3/minimap2) (map)
- [samtools](https://github.com/samtools/samtools) (map)
- [blastn](https://www.ncbi.nlm.nih.gov/books/NBK569861/) (blast_annot_batch)
- [makeblastdb](https://www.ncbi.nlm.nih.gov/books/NBK569861/) (blast_annot_batch)
- [MAFFT](http://mafft.cbrc.jp/alignment/software/) (dr_logo)

After installation, run the **test** module below, If you don't see any errors in the output, means all programs have been successfully installed. trapo-seq expects external programs to be available in **$PATH**.

**Test** (It will use the test data as input under (/test_data/*) and creates outputs in same folder):
```
conda activate trapo-seq-env
python trapo-seq.py --test
```

# Usage

First, assesing the quality of sequencing output is crucial. It would be appropriate to examine the length and quality distributions of the reads to check fragmentation.

[The general structure](#flowchart) of the pipeline (relation of modules and input/output files can be seen as flowchart)

# Input Data Analysis

## Read Length Histograms

If plasmids are extracted with sufficient purity, most of the reads are expected to be close to the length of the plasmid.
### Requirements
- **read_histogram** module is enough to handle basic settings, but there are other tools you can also use such as [NanoPlot](https://github.com/wdecoster/NanoPlot)

```
python trapo-seq.py read_histogram -f data/barcode01/barcode01_raw.fastq -o data/barcode01/read_len_hist
```
- **-f** fastq file path
- **-o** output file path

 The desired scenario is that there are mobile element insertions into the plasmids so the majority of reads will be longer than the original plasmid. If there are shorter reads than, this may be due to: 
 - fragmentations while doing library preparation, 
 - the enzyme cuts more than one restriction site
 - the inserted mobile genetic elements have restriction site. 
 
 For extra information please visit the [docs page](/docs#readme). So, it is crucial to pick restriction enzyme that does not cut mobile genetic elements and has only one restriction site on the plasmid. If the majority of the reads are not longer than the the original length of the plasmid, this may be explained by scenarios such as **mutations** on selective genes in plasmid or the mobile element **excisions**.

As an example output (Figure 1):

![example_output](/images/read_len-1.png)
**Figure 1:** Non weighted histogram graph of plasmid reads. The original size of this plasmid is 7kb, and there is a peak around 7 kb. Also, there is one peak around 8 kb. So, it can be said that insertion frequency is low, most of the plasmid reads do not have insertions. The peak around 5kb probably caused by the several unwanted enzyme cuts (insertions sometimes can cause new restriction sites).

## Filtering

It is essential to filter shorter reads because they are not helpful most of the time. We want long reads to capture the traposons as a whole, so it may be better to filter reads lower than 2kb length (average transposon length*2 can be used as prediction). But this threshold value can change according to sequencing data (length and quality distributions) and plasmid length. In an ideal world, it may be better to eliminate reads shorter than the plasmid length, however sometimes [different cases](/docs#undesired-cases) can happen.

### Requirements
- Basecalled and barcode demultiplexed .fastq files
- If you want to filter only length-based (you do not need an external program like [FiltLong](https://github.com/rrwick/Filtlong)) you can use **read_histogram** module:

```
python trapo-seq.py filter -f data/barcode01/barcode01_raw.fastq -o data/barcode01/read_len_hist -l 2000
```
- **-f** fastq file path
- **-o** output file path
- **-l** filtering threshold 


# Alignment

## Map

### Requirements
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)

You need to have minimap2 and samtools in the your path. So, try these commands in the terminal before going further:
```
minimap2 --version
samtools --version
```

If there is no error, you can continue. Please, be careful about **folder structure and naming**. (Suggested naming: barcode01/barcode01.fastq) or you can use the condition names which barcodes corressponding. **Folder names and fastq names must be matched**!
```
|-- data
  |-- barcode_01
    |-- barcode_01.fastq
  |-- barcode_02
    |-- barcode_02.fastq
```

```
python trapo-seq.py map -d data -p data/plasmid.fasta -g data/genome.fasta --force
```
- **-d** directory which barcode folders having fastq files 
- **-p** plasmid sequence data path (*.fasta)
- **-g** genome sequence data path (*.fasta)
- **--force** to overwrite results and temp files

In this script, the reads are mapped to the plasmid and also the genome. Using the both plasmid and genome is significant, because if both of them have alignments in a read that means there is a jumping gene from genome to the plasmid (If the plasmid and genome do not have homologous portions, and if the reads are not chimeric). Therefore, transposition ratio can be calculated.

If genomic DNA contamination is to be examined, a careful approach should be taken because transposed plasmids carry portions of genomic DNA. If a read do not have any plasmid portion, the read can be labeled as a DNA contamination.

During library prepation or in the post-analysis in sequencing protocol, different barcodes can get mixed or incorrect demultiplexing can happen (the amounts are expected to be low). This is also important to take into account before moving further especially if different strains or different plasmids are sequenced in same batch.

You can also look the docs section [for taking extra info about sam file and CIGAR strings](/docs#understanding-cigar-string). If you want to analyse your data with different tools than minimap2, the insertion finding algorithm may need some adjustments. I tested it with [bwa-mem](https://github.com/lh3/bwa), and it is working fine but I do not know others.


# Analysis
## Insertion Finder

This module will find insertions using bam files created by the map module. Just give the data folder as input and it will detect and analyze bam files.

```
python trapo-seq.py insert_finder_batch -d data
```
- **-d** directory which folders having bam files

Output of this step creates a tabular file (insertion_from_bam.tab) in the format of:

| read_id | query_start | query_end | insertion_length | ref_pos | quality | insert_type | is_reverse |
|---------|-------------|-----------|------------------|---------|---------|-------------|------------|
| read_1  |      100   |    1200       |    1180       | 4093      |   98  |  [IN, SC]   |   True     |
</br>

You can also look the docs section [to see different insert_types: IN, SC](/docs#extra-case-soft-clips)

## Manually curated Transposons

If there is no data about tranposons, you can still continue to this pipeline without entering **--is_fasta** argument (the program will used a dummy_is_file.fasta):

```
dummy_is_file.fasta:

>IS3_ISAba66
TGAACCGTACCGG...
```
And then you can investigate the insertions and corresponding genome locations to extract transposons and annotate through blasting with [ISFinder](https://www-is.biotoul.fr/index.php).

However, I recommend scanning the genome before the next step using tools like [ISFinder](https://www-is.biotoul.fr/index.php) or [ISescan](https://github.com/xiezhq/ISEScan) to identify transposons, and then proceeding with further analysis. And creating **IS_curated.fasta** file like this as input for **--is_fasta**: Sequence IDs must be in the format **(>ISfamily_ISName)** such as (IS5_ISPs1). 

```
IS_curated.fasta:

>IS5_ISPs1
ATGGCTGATGACAAD...
>IS3_ISPs5
GCTGATGACAAD...
```
Then, it is ready for the next module.

## Annotation of Insertions

In this step, the extracted insertions will be searched in the genome and manually curated IS FASTA file, and the results will be presented in a tabular format. Also, the output of this module can be used to label transposons. 

### Requirements
- [blastn and makeblastdb](https://www.ncbi.nlm.nih.gov/books/NBK569861/) in PATH

```
blastn -version
makeblastdb -version
```

```
python trapo-seq.py blast_annot_batch -d data -g data/genome.fasta --is_fasta data/IS_curated.fasta --threads 12 --no_temp
```
- **-d** directory which barcode folders having tab files (*insertion_from_bam.tab)
- **-g** genome sequence data path (*.fasta)
- **--is_fasta** manually curated and labeled ISes (*.fasta) Sequence IDs must be in the format (>ISfamily_ISName) such as (IS5_ISPs1). Please have a look [Manually curated Transposons](#manually-curated-transposons)
- **--threads** how many threads can be possible (default 2)
- **--no_temp** do not delete temp files, can be helpful for debugging 

The output of this module is a tabular file (*best_alignment.tab) in the format of:

| Query ID | Subject ID | Identity (%) | Score | E-value | Query Start | Query End | Subject Start | Subject End | Note | Explained | ref_pos | is_Reverse |
|----------|------------|--------------|-------|---------|-------------|-----------|---------------|-------------|------|-----------|---------|------------|
|   read_1    |   IS30      |    98   |   456   |   0   |     12     |   800    |    10    |   780   |  IS_DB  |  98     |   4000   |    True     |
|   read_2    |   genome_id      |    99   |   900   |   0   |     13     |   2100    |    15    |   2500   |  Genome  |  99     |   4250   |    False     |
|   read_3    |   no_blast_hit      |    N/A   |   N/A   |   N/A   |     0     |   0    |    0    |   0   | Contamination  |  0     |   5165   |    True     |
</br>

The most significant columns in the table are **Explained**, **Note**, and **ref_pos**. The **"Explained"** value indicates how well the insertion aligned(0-100), demonstrating its relevance and matching with the manually curated database (**IS_DB**) or (**Genome**). If it aligns with a part of the **genome** rather than the database ((**IS_DB**), it can potentially represent **a novel insertion**. The **"ref_pos"** value is the insertion point on the reference plasmid, marking the exact location where the insertion occurs. If the insertion is labeled as **Contamination**, the insertion is not coming from genome and can be an artifact related to sequencing, barcoding issues.

## Examining the novel inserts

One of the key advantages of this pipeline is its ability to identify **novel insertions** without relying exclusively on a database, allowing it to detect previously unreported insertions. For example, in the output table above, regions labeled as **'Genome'** in the **'Note'** column represent segments not found in the manually curated IS database. If a certain genomic region (The **'Subject Start' and 'Subject End'** columns indicate the location within the genome.) is observed in multiple reads and has **high 'Explained' (%)** values, you may want to extract this genomic section for further examination. You can search it against IS databases like [ISFinder](https://www-is.biotoul.fr/index.php) or use [BLAST on the NCBI website](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to see which genes are in. 

## Stats of Transposons

This module can be used to summarize output tabular files (*best_alignment.tab) from the previous module.

```
python trapo-seq.py is_stat -d data -o data/is_stats.rcp
```
- **-d** directory which barcode folders having tab files (best_alignment.tab)
- **-o** output file (*.rcp) (a human-readable and parsable format to summarize stats)

An example output file (*.rcp) of this module:

  - **>** -> name of condition
  - **?** -> number of transposons
  - **IS_Family, IS_Name, %, #** -> IS Family, IS Name, Percentage, Count 
  - **---** -> splitting conditions
```
is_stats.rcp:

>barcode01
?40
IS5, ISPa41, 50.00, 20
IS3, ISPen2, 25.00, 10
IS5, ISPs1, 25.00, 10
---
>barcode02
?200
IS3, ISPen2, 65.00, 130
IS5, ISPs1, 25.00, 50
IS5, ISPa41, 5.00, 10
IS5, ISPa26, 5.00, 10
```

The next step (**heatmap module**) needs this exact format (.rcp) as input. If you have any other transposon distribution data, you can convert to (.rcp) format to run **heatmap module**.

## Heatmap of Transposons

According to the output of previous module, heatmap module will create a heatmap of low and high frequency transposons. 

```
python trapo-seq.py heatmap -r data/is_stats.rcp -o data/heatmap.pdf -t data/heatmap.tsv
```
- **-r** stat file (.rcp)
- **-o** output pdf heatmap graph file (.pdf)
- **-t** output file (*.tsv) (to summarize heatmap data)

Example output:

![example_heatmap](/images/example_heatmap.png)
**Figure 2:** Distribution of transposons in different conditions (barcodes).
**Hint**: You can order the barcodes in .rcp file to change the order in the heatmap.

## Insertion Coordinates

The output (insertion_from_bam.tab) from [Annotation of Insertions](#annotation-of-insertions) can be used to extract insertion points (**ref_pos** in the .tab file). Also, [dr_logo](#dr-sequence-logos) module creates insertion locations plot.
## DR Finder

Transposition event can create **direct repeats (DR)**. Direct repeats in transposons are short, identical sequences of DNA that flank the insertion site of a transposon. These repeats are generated during the transposition process (**not generated in all transposition events!**), leaving overhangs that are subsequently filled in, creating duplicated sequences on either side of the transposon. Blasting reference sequence to the reads **can reveal overlapping regions** near the insertion point (**ref_pos**):  **The possible DRs.**

```
python trapo-seq.py dr_finder -d data -p data/plasmid.fasta
```
- **-d** directory which barcode folders having tab files (best_alignment.tab)
- **-p** trap plasmid sequence (.fasta)

An example output of DR Finder module (insertions.csv):

| Read ID	| Subject ID	| Query Start |	Query End |	Start Subject	| End Subject	| Insertion Point |	Repeat Length	| Overlap Sequence |
|----------|------------|-----------------|---------------|----------------|-------------|----------------|---------------|-------------|
| read1	| IS5_IS15 |    	8726	         | 9883	         | 4854	           | 4856        |	4854         |	3           |	CTA |
| read2	| no blast hit	|5416	            |6742           |	4839	         | 4850	       | 4839         |	12	         | TGCTTGGTTATG |
</br>

## DR Sequence Logos

Some transposons exhibit a high degree of specificity for insertion points, meaning they consistently **target particular sequences** or structural features in the host genome. This specificity can result in **conserved direct repeats** flanking the insertion as the transposase enzyme recognizes and targets **specific DNA motifs**. Conversely, some transposons display more flexibility in insertion sites, leading to less conserved or even absent direct repeats. This variability in insertion specificity and repeat conservation depends on the transposon's mechanism and the interaction between the transposase enzyme and the host genome. Such differences are significant in understanding the evolutionary impact and functional role of transposons in various genomic contexts

### Requirements
  - [MAFFT](http://mafft.cbrc.jp/alignment/software/) (MAFFT is used to align sequences, so it needs to be in the PATH)
```
mafft --version
```
```
python trapo-seq.py dr_logo -d data/insertions -p data/plasmid.fasta -o data/dr_logos
```
- **-d** directory which barcode folders having tab files (best_alignment.tab)
- **-o** output pdf heatmap graph file (.pdf)
- **-t** output file (*.tsv) (to summarize heatmap data)

Conservation of overlapping regions can reveal direct repeats and insertion point's specificity of transposons. To reveal the conservation of direct repeats; sequence logos with aligning the overlapping regions are made in this module for each transposon identified in the previous modules. The length of the logo depends on the alignment length*.

An example output for one of the transposon:

![example_sequence_logo](/images/example_sequence_logo.png)
**Figure 3:** As an example, TGCTTGGTTATG is a 12 bp insertion site for this transposon (**TGCTTGGTT** is highly conserved). The closest transposon in ISFinder has DR as TGCTGAATT (close but different), so new insertion points can be found with the help of [**trapo-seq**](https://github.com/recepcanaltinbag/trapo-seq).

*(for DR alignments of each transposon; most common length is chosen)


![example_insertion_points](/images/example_insertion_points.png)
**Figure 4:** Insertion coordinates. Main trend is 4000-5000 bp locations of this insertion. 


## Kde Plots 

You can visualize **the insertions from the BAM file** using the **kde_plot** module. This module can be run after the [Insertion Finder](#insertion-finder) step, as it requires the insertion output. It can be used to examine the distribution of insertions under different conditions. Instead of showing the exact insertion points, it displays **the general trend**; total number of bases in each read coming from the plasmid and genome.

```
python trapo-seq.py kde_mobile -b data/barcode01/03_sorted_mapped_to_plasmid.bam -o kde_plot
```
- **-b** bam file with prefix: 03*.bam
- **-o** output name (.pdf)

An example output:

![example_kde](/images/example_kde.png)
**Figure 5:** Kde Plot of insertions.

## Possible Excisions

Transposons can sometimes excise themselves after insertion, a process known as excision, in such cases the number of insertions may appear low. It may be useful to view insertions and deletions as weighted because transposons can create deletions while exiting.

```
python trapo-seq.py in_del_plot -i data/barcode01/insertion_from_bam.tab -b data/barcode01/03_sorted_mapped_to_plasmid.bam -o in_del_plot
```
- **-i** insertion file (.tab)
- **-b** bam file with prefix: 03*.bam
- **-o** output name (.pdf)

An example output:

![example_kde](/images/example_in_del_plot.png)
**Figure 6:** In-del Plot for insertions and deletions. Weighted insertions are higher when compared to deletions, so there is no significant excision. And, most of the insertions are between 3800-5000 bp region, for that case it is expected (This region contains the selective gene) but can be changed from experiment to experiment.   

## Mutation and Variant Analysis 

If  the insertion deficiency cannot be revealed in the in-del plot, you can perform mutation and variant analysis. For this, you can use an external program (SNP, variant analysis) for long reads. You can also use the (*_insertion_table.tsv, *_deletion_table.tsv) files created after the in_del_plot module.

## Map Distribution

The **map-dist** module can be used to see the distributions of reads based on their location in the genome and plasmid, so the mobile sections of the genome can be seen.

```
python trapo-seq.py map_dist -b data/barcode01/03_sorted_mapped_to_plasmid.bam -o map_dist
```
- **-b** bam file with prefix: 03*.bam
- **-o** output name (.pdf)

![example_map_dist](/images/example_map_dist.png)
**Figure 7:** Read distributions on plasmid and genome. For the plasmid; around 3000 there is a certain peak which is actually caused by restriction enzyme, and also peaks around 4000 which is caused by transpositions. For the genome; around ~2.4 Mbp and ~6.8 Mbp, sections are very mobile. This output can not be used for exact locations, but can give an idea about general trend and can be helpful for the undesired conditions.

# FAQ

## How Can I Cite trapo-seq?

Hopefully the paper will come soon :), for now, this github page. 

## License

This software in under [MIT License](https://github.com/recepcanaltinbag/trapo-seq/tree/main?tab=MIT-1-ov-file). Copyright (c) 2024 Recep Can Altınbağ

## Flowchart

![example_kde](/images/example_flowchart.png)

### Requirements based on each module 

  - minimap2
  - samtools
- insert_finder_batch
  - pysam
- blast_annot_batch
  - biopython
  - blastn, makeblastdb
- is_stat
  - pandas
  - biopython
- heatmap
  - seaborn
  - matplotlib
  - pandas
  - numpy
- dr_finder
  - biopython
  - pandas
- dr_logo
  - MAFFT
  - pandas
  - matplotlib
  - seaborn
  - biopython
  - numpy
  - logomaker
- in_del_plot
  - matplotlib
  - pysam
  - pandas
- kde_mobile
  - matplolib
  - pysam
  - seaborn
- map_dist
  - pysam
  - matplotlib
