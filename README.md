# AmpliconHaplotypes
Analysis of haplotypes in amplicon data - only tested on MinIon/PacBio so far, but will extend to paried-end Illumina.
The Main.py does the following:
1. Takes all bams files and validates
2. For each bam file in the directory speficied in config/configname.py validates bam file
3. Determines the regions of reference genome (reference itself not required) to which the reads from bams map
4. Filters out spurrious regions (settings in configs) and merges overlapping regions (settings again in config) into one
5. For each region creates a numpy pickle file which contains a binary matrix with rows as reads and columns containing 1,2,3,4,5 corresponding to ACTG- 
  a. this works with sklearn definition of hamming distances, but not with jaccard or other boolean distance measures.
6. Matrix has column "Sample" which specifies which samples the read (row) came from. Useful for colouring the points on plots later.

To run first create a config file in configs directory. Copy existing one to use as template and change as needed.  
Second, in Main.py change line that looks like "import configs.CONFIGFILENAME as Config" to be "import configs.CONFIGFILEYOUCREATED as Config"  
Finally, run Main.py
The results can be visualised in AmpliconReadClustering.ipynb notebook.  

Relies on:  
numpy  
pandas  
umap (McInnes, L, Healy, J, 2018) though this can be replaced with sklearn implementaiton of TSNE or other methods.  
Be careful not to mixup umap and umap-learn. The script requires the latter.  
matplotlib  
sklearn  
BioPython  
pysam  