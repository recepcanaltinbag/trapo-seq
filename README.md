![example_output](/images/logov4.png)

### Traposome Sequencing with Nanopore Data

The biggest advantage of the long read technology for traposome analysis is that it can read the plasmids as a whole, just one read can allow the analysis of mobile element that have entered the plasmid. This also helps to obtain high throughput traposome distribution. 

First, it would be appropriate to examine the length and quality distributions of the reads to check fragmentation.

## Input Data Analysis

## 1. Read Length Histograms

## Requirements
- Basecalled and barcode demultiplexed .fastq files
- [NanoPlot](https://github.com/wdecoster/NanoPlot)
- If you want to visualize without using NanoPlot you can use the Python script under extra folder (extra/00_read_histograms.py)

```
NanoPlot --fastq <fastq_file> -o <output_folder>
```

If plasmids are extracted with sufficient purity, most of the reads are expected to be close to the length of the plasmid. The desired scenario is that there are mobile element insertions into the plasmids so the majority of reads are longer than the original plasmid. If there are shorter reads than, this may be due to: fragmentations while doing library preparation, the enzyme cuts more than one restriction site or the inserted mobile genetic elements have restriction site. So, it is crucial to pick restriction enzyme that does not cut mobile genetic elements and has only one restriction site on the plasmid. If the majority of the reads are not longer than the the original length of the plasmid, this may be explained by scenarios such as mutations on selective genes in plasmid or the mobile element excisions.

As an example output (Figure 1):

There are fragmented reads around 1kb, but most of reads are longer than 5kb. Also, the original size of this trap plasmid is 7kb, but there are reads more than 7 kb such as around 8-10 kb. So, it can be said that there are two main insertions of different lengths, and this data is good for next steps.

![example_output](/images/1_LengthHistogramv2.png)

Figure 1: Non weighted histogram graph of plasmid reads.

## 2. Filtering


## Requirements
- Basecalled and barcode demultiplexed .fastq files
- [FiltLong](https://github.com/rrwick/Filtlong) or 
- If you want to filter only length-based (you do not need FiltLong) you can use the Python script under extra folder (extra/01_filtering_based_on_len.py)

```
filtlong --min_length 1000 merged.fastq | gzip > filtered_1kb.fastq

```

## 3. Alignment

## Requirements
- Filtered .fastq files
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)

In the directory where the demultiplexed barcode folders are located, you can use scripts/01_data_prep.sh. In this script, the alignment section has the code for mapping the fastq files inside the folders to plasmid and also genome.


## 4. Insertion Finder









