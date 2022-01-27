class configData:
    def __init__(self):
        self.wd = "/mnt/storage5/anton/Amplicon/LeenIllumina/"  #directory in which any work will be done file will be
        self.bamDir = "/mnt/storage9/leen/amplicon-seq/2021_11_08_Ligation/"  #directory which will be scanned for bam files
        self.tempDir=self.wd+"/temp/" #directory for temp files
        self.minRegionDepth=10 #any aggragated region with fewer than this reads will be dropped from further analysis
        self.maxDistanceBetweenReads=20 #any regions where starts and ends are both fewer than this bases apart are considered one region
        self.minPosDepthAbsolute=5 #to accomodate high error rate in MinIon, the positions in nucleotide matrix supported by fewer than this reads are dropped
        self.minPosDepthRelative=0.01 #to accomodate high error rate in MinIon, the positions in nucleotide matrix supported by lower than this % of total reads dropped
        self.MaxReadsToOutput=1000 #use for downsampling of very large datasets        