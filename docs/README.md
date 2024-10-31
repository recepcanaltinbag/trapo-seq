### Overview of the Pipeline



### Enzyme Cut Histograms



#### Initial Case:

If the enzyme cuts at a single site, we would expect the reads to predominantly accumulate at a point corresponding to the plasmid length in the initial case (1): 

![example_output](/images/01_EnzymeCutsHistogram_1.svg)

#### Perfect Case:

Following transposition, multiple insertions into the plasmid will cause an increase in read lengths, which can be observed by examining the gel images or analyzing the histogram of the reads (2):

![example_output](/images/01_EnzymeCutsHistogram_2.svg)

#### Undesired Cases:

If the inserted transposed sequence contains a restriction site, this will cause the reads to become fragmented, which can be identified from gel images or read histograms (3). 

Additionally, if there are multiple restriction sites on the plasmid, a similar histogram pattern may appear (4). Therefore, it is crucial to choose the restriction enzyme carefully to prevent read fragmentation.

Additionally, DNA fragmentation can occur during plasmid extraction and library preparation processes. To minimize this, carefully following protocols is essential, as such fragmentation could negatively impact the analysis (5).

![example_output](/images/01_EnzymeCutsHistogram_3.svg)


### Understanding CIGAR string

For the read alignments to the reference, .sam format is mostly used and in this pipeline insertions are exracted with the help of analyzing .sam files.

Different alignment programs can use varying coordinate systems for start and end positions. In my analyses, I used `pysam`, which operates with a 0-based coordinate system. Therefore, I am displaying all coordinates in a 0-based format for consistency.

**Query Start = 0** or if there is a Soft Clip, **Query Start = S**

**Query End = Query Start + (Total size of M and I)**

**Alignment Length: Matches and Deletions (0,2)**

Insert type can be IN or SC. IN means that the alignment program can successfully classify insertion as I. However, for some cases, alignment programs may not be classify them and fragment the alignment into two queries. So, the insertion can be find in soft-clipped region (SC). I will discuss this in next section. 

The start positions can be based on either 1-based or 0-based indexing. In this example, it is 0-based. For 1-based indexing, the positions would be [3,7], while for 0-based indexing, they would be [2,7).

For example, minimap2 classifies a section as a insertion, but bwa-mem classifies 
same section as a soft clip, also there are other type of labels such as 
hard clip, exact match etc. So, it can be useful to analyze what kind of different CIGAR alphabet
is used in the .sam file.

![example_output](/images/CIGAR1.svg)

![example_output](/images/CIGAR2.svg)

![example_output](/images/CIGAR3.svg)


### Extra Case: Soft Clips

![example_output](/images/02_softClips_1.svg)

![example_output](/images/02_softClips_2.svg)




