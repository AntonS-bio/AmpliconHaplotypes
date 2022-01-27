# AmpliconHaplotypes
Analysis of haplotypes in amplicon data - only tested on MinIon/PacBio so far, but will extend to paried-end Illumina.
The Main.py does the following:
1. Takes all bams files and validates
2. For each bam file in the directory speficied in config/configname.py validates bam file
3. Determines the regions of reference genome (reference itself not required) to which the reads from bams map
4. Filters out spurrious regions (settings in configs) and merges overlapping regions (settings again in config) into one
5. For each region creates a numpy pickle file which contains a binary matrix with rows as reads and columns as ACTG- for each position
  5.a. Matrix has column "Sample" which specifies which samples the read (row) came from. Useful for colouring the points on plots later.

The matrices can be analysed in AmpliconReadClustering notebook. 
