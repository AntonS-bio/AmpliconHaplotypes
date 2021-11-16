from matplotlib.pyplot import contour
import numpy as np
import sys
sys.path.insert(1, '/mnt/storage5/anton/generalScripts/utilities/')
from Bio import SeqIO
import subprocess
import pysam as ps
from os import chdir, listdir
from os.path import isfile, splitext, join
import determineAlignmentRegions #this returns {chr: "regionstart_regionsend"} dictionary for all bam files bassed to it where regions are non-overlapping
#the point is to generate regions of alignment purely from bam files. Currently works for single end reads.
import writeFastas #this outputs fasta files and aligns them with mafft. It returns file names
import fastaToNumpy #this converts MSAs into numpy matrices, positions which satisfy both <1% of reads and <5 reads are dropped
import configs.Leen as Config
config=Config.configData()


##DISUSED
##        subprocess.run("samtools view -h "+bamDir+sample+" | java -jar /mnt/storage5/anton/generalScripts/jvarkit/dist/sam2tsv.jar > sam2tsvMatrix.tsv", shell=True, stdout=subprocess.PIPE)


wd="/mnt/storage5/anton/Amplicon/Leen/"
bamDir="/mnt/storage9/leen/Minion/mapping/"
tempDir=wd+"/temp/"
unvalidated_samples = [f for f in listdir(bamDir) if isfile(join(bamDir, f)) and splitext(f)[1]==".bam"]
#validate files
validsamples=[]
for sample in unvalidated_samples:
    try:
        bam = ps.AlignmentFile(bamDir+sample, "rb",check_sq=False)
        bam.close()
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
#sys.exit()


chdir(wd)
sampleCounter=0
allResults=[]

sampleChrRegionReads={} #key=sample, value={key=chr, value={key=region, value=[read sequences]}}
#ChrRegionReads={} 
for sample in validsamples:
    sampleChrRegionReads[sample]={}
    for chr in intervalDepths:
        sampleChrRegionReads[sample][chr]={}
        for region in intervalDepths[chr]:
            sampleChrRegionReads[sample][chr][region]=[]

readStartEndToRegion={} #this is for speed as most reads will start and end at the same coordinates. key=(readstart,readend), value=(regionstart,regionend)
for sample in validsamples:
    bam = ps.AlignmentFile(bamDir+sample, "rb",check_sq=False)
    for read in bam.fetch():
        if bam.get_reference_name(read.reference_id) not in sampleChrRegionReads[sample]:
            continue #sometimes spurrious hits are generated on non-target chromosomes
        if (read.reference_start, read.reference_end) not in readStartEndToRegion:
            #determine which region the read falls into
            for region in sampleChrRegionReads[sample][bam.get_reference_name(read.reference_id)].keys():
                if abs(read.reference_start-region[0])<config.maxDistanceBetweenReads and abs(read.reference_end-region[1])<config.maxDistanceBetweenReads:
                    readStartEndToRegion[(read.reference_start, read.reference_end)]=region
                    break
            if (read.reference_start, read.reference_end) not in readStartEndToRegion: #this checks if read in in target regions, if not it is skipped
                continue
        sampleChrRegionReads[sample][bam.get_reference_name(read.reference_id)][readStartEndToRegion[(read.reference_start, read.reference_end)]].append(read.query_alignment_sequence.replace("N","-"))

#the reads from samples have been collected into sampleChrRegionReads[sample][chr][region]=[read seqeuences]
#print them to fasta and allign using mafft

alignedFiles={} #key=sample, value={key=chr, value={key=region, value=aligned filename}}
alignedFiles=writeFastas.run(sampleChrRegionReads, "allsamples", config)


#generate numpy matrix from the MSA fasta file
MSAfiles=set()
for sample in alignedFiles:
    for chr in alignedFiles[sample]:
        for region in alignedFiles[sample][chr]:
            MSAfiles.add(alignedFiles[sample][chr][region])



fastaToNumpy.run(MSAfiles, config)

sys.exit()


#         if readCounter>10:
#             allResults=allResults+analyseBamTsvMatrix.run(matrix,key,sample)

#             np.savetxt(outDir+key+".tsv", matrix, delimiter="\t", fmt='%d')

# with open(outDir+"Leen_PF_frequencies.tsv", "w") as output:
#     for i in range(len(allResults)):
#         output.write('\t'.join(str(f) for f in allResults[i])+"\n")