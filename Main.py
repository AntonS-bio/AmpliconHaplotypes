from genericpath import exists
from matplotlib.pyplot import contour
import numpy as np
import sys
#sys.path.insert(1, '/mnt/storage5/anton/generalScripts/utilities/')
from Bio import SeqIO
import subprocess
import pysam as ps
from os import chdir, listdir
from os.path import isfile, splitext, join
import determineAlignmentRegions #this returns {chr: "regionstart_regionsend"} dictionary for all bam files bassed to it where regions are non-overlapping
#the point is to generate regions of alignment purely from bam files. Currently works for single end reads.
import writeFastas #this outputs fasta files and aligns them with mafft. It returns file names
import fastaToNumpy #this converts MSAs into numpy matrices, positions which satisfy both <1% of reads and <5 reads are dropped
import collectReads #loops through bams and collects reads which will be output to fasta, also applies downsampling
#import configs.Leen as Config
#import configs.LeenIllumina as Config
#import configs.Ashley as Config
import configs.LeenNewBarcodes as Config
#import configs.Holly as Config

config=Config.configData()


wd=config.wd
bamDir=config.bamDir
tempDir=config.tempDir
unvalidated_samples = [f for f in listdir(bamDir) if isfile(join(bamDir, f)) and splitext(f)[1]==".bam"]
#validate files
validsamples=[]
for sample in unvalidated_samples:
    try:
        bam = ps.AlignmentFile(bamDir+sample, "rb",check_sq=False)
        bam.close()
        if not exists(bamDir+sample+".bai"):
            subprocess.run(f'samtools index {bamDir}{sample}', shell=True, executable="/bin/bash", stdout=subprocess.PIPE)
        validsamples.append(sample)
    except:
        #sample bam is invalid
        print("Invalid bam: "+sample)
        #pass

#alignmentRegions={chr: [chrregionstart_chrregionend,chrregionstart_chrregionend]}}
#intervalDepths={chr: {interval: depth}}
#validsamples=validsamples[0:5]
alignmentRegions, intervalDepths=determineAlignmentRegions.run(validsamples,config)
print(alignmentRegions)
print(intervalDepths)


chdir(wd)
sampleCounter=0
allResults=[]

sampleChrRegionReads=collectReads.getReadsFromBams(validsamples, intervalDepths, config) #key=sample, value={key=chr, value={key=region, value=[read sequences]}}


alignedFiles={} #key=sample, value={key=chr, value={key=region, value=aligned filename}}
alignedFiles, MSArunTimes=writeFastas.run(sampleChrRegionReads, "allsamples", config)
print("Alignemnt runtimes: longer than most means likely problem with alignemnt")
print(MSArunTimes) #shows if some alignments took too long to run. Long run time indicates problems

#generate numpy matrix from the MSA fasta file
MSAfiles=set()
for sample in alignedFiles:
    for chr in alignedFiles[sample]:
        for region in alignedFiles[sample][chr]:
            MSAfiles.add(alignedFiles[sample][chr][region])



fastaToNumpy.run(MSAfiles, config)

