![example_output](/images/logov4.png)

### Traposome Sequencing with Long Read Technologies

trap-seq is a traposome sequencing pipeline for plasmid sequencing with long read technologies. The biggest advantage of the long read technology for traposome analysis is that it can read the plasmids as a whole, just one read can allow the complete sequence of the mobile elements  that have trapped by the plasmid. This also helps to obtain high throughput traposome data. 

First, assesing the quality of plasmid DNA library and the sequencing output is crucial. It would be appropriate to examine the length and quality distributions of the reads to check fragmentation.

## Input Data Analysis

## 1. Read Length Histograms

## Requirements
- Basecalled and barcode demultiplexed .fastq files
- [NanoPlot](https://github.com/wdecoster/NanoPlot) or If you want to visualize without using NanoPlot you can use the [Python script](/extra/00_read_histograms.py) under extra folder.

```
NanoPlot --fastq <fastq_file> -o <output_folder>
```

If plasmids are extracted with sufficient purity, most of the reads are expected to be close to the length of the plasmid. The desired scenario is that there are mobile element insertions into the plasmids so the majority of reads are longer than the original plasmid. If there are shorter reads than, this may be due to: fragmentations while doing library preparation, the enzyme cuts more than one restriction site or the inserted mobile genetic elements have restriction site. For extra information please visit the [docs page](/docs#readme). So, it is crucial to pick restriction enzyme that does not cut mobile genetic elements and has only one restriction site on the plasmid. If the majority of the reads are not longer than the the original length of the plasmid, this may be explained by scenarios such as mutations on selective genes in plasmid or the mobile element excisions.

As an example output (Figure 1):

There are fragmented reads around 1kb, but most of reads are longer than 5kb. Also, the original size of this trap plasmid is 7kb, but there are reads more than 7 kb such as around 8-10 kb. So, it can be said that there are two main insertions of different lengths, and this data is good for next steps.

![example_output](/images/1_LengthHistogramv2.png)

Figure 1: Non weighted histogram graph of plasmid reads.

## 2. Filtering

To eliminate the shorter reads, it is essential to filter. We want long reads to capture the traposons as a whole, so it may be better to filter reads lower than 2kb length. But this threshold value can change according to the histograms and your plasmid length. In an ideal world, it may be better to eliminate reads shorter than the plasmid length, however sometimes [different cases](/docs#undesired-cases) can happen.

## Requirements
- Basecalled and barcode demultiplexed .fastq files
- [FiltLong](https://github.com/rrwick/Filtlong) or If you want to filter only length-based (you do not need FiltLong) you can use the [Python script](/extra/01_filtering_based_on_len.py) under extra folder.

```
filtlong --min_length 1000 merged.fastq | gzip > filtered_1kb.fastq

```

## 3. Alignment

## Requirements
- Filtered .fastq files
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)

In the directory where the demultiplexed barcode folders are located, you can use [the data prep script](/scripts/01_data_prep.sh). In this script, the reads in the different folders are mapped on your plasmid and also the genome you used. Using the both plasmid and genome is significant, because if both of them have alignments in a read that means there is a jumping gene from genome to the plasmid (If the plasmid and genome do not have homologous portions, and if the read are not chimeric). Therefore, transposition ratio can be calculated.

If genomic DNA contamination is to be examined, a careful approach should be taken because transposed plasmids carry portions of genomic DNA. If a read do not have any plasmid portion, it can be labeled as a DNA contamination.

During barcoding, different barcodes can get mixed up with other barcodes or that incorrect classification occurs after demultiplexing. This is also important to take into account before moving further especially if barcodes are from different strains or different plasmids.

You can also look the docs section [to understand sam file and CIGAR strings](/docs#understanding-cigar-string). If you want to analyse your data with different tools than minimap2, the insertion finding algorithm can need some adjustments. I tested it with [bwa-mem](https://github.com/lh3/bwa), and it is working fine but I do not know others.


## 4. Insertion Finder









