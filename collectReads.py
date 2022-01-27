import pysam as ps
import random


def getReadsFromBams(validsamples, intervalDepths, config):

    sampleChrRegionReads={}
    for sample in validsamples:
        sampleChrRegionReads[sample]={}
        for chr in intervalDepths:
            sampleChrRegionReads[sample][chr]={}
            for region in intervalDepths[chr]:
                sampleChrRegionReads[sample][chr][region]=[]

    readStartEndToRegion={} #this is for speed as most reads will start and end at the same coordinates. key=(readstart,readend), value=(regionstart,regionend)
    processedReads=set()#some reads mapped to multple locations, only primary will be used
    for sample in validsamples:
        pairedReadIDs={}
        pairedReadStartEnd={}
        pairedReadSequences={}
        bam = ps.AlignmentFile(config.bamDir+sample, "rb",check_sq=False)
        for read in bam.fetch():
            readStart=read.reference_start
            readEnd=read.reference_end
            if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped and not  read.is_secondary: # is_secondary is false when read is unmapped
                if read.query_name in pairedReadIDs:
                    # this is the second of read pairs
                    # the two reads can overlap, sometimes perfectly
                    # make final read a composite of both                
                    fullRead="-"*(max(readStart,readEnd,pairedReadStartEnd[read.query_name][1])-min(readStart,readEnd,pairedReadStartEnd[read.query_name][0]))
                    if min(readStart,readEnd)<pairedReadStartEnd[read.query_name][0]:
                        #the current read comes first in the pair
                        fullRead=read.query_alignment_sequence.replace("N","-")+fullRead[len(read.query_alignment_sequence.replace("N","-")):]
                        fullRead=fullRead[ 0: len(fullRead)-pairedReadStartEnd[read.query_name][0] ] + pairedReadSequences[read.query_name]
                        pairedReadSequences[read.query_name]=fullRead
                    else:
                        fullRead=pairedReadSequences[read.query_name]+fullRead[len(pairedReadSequences[read.query_name]):]
                        fullRead=fullRead[0:len(fullRead)-len(read.query_alignment_sequence.replace("N","-"))] + read.query_alignment_sequence.replace("N","-")
                        pairedReadSequences[read.query_name]=fullRead
                    pairedReadIDs[read.query_name]=2
                    #the second read of the pair has been processed, so set the start/end to be values for the whole pair to be consistent with single end reads
                    readStart=min(readStart,readEnd,pairedReadStartEnd[read.query_name][0])
                    readEnd=max(readStart,readEnd,pairedReadStartEnd[read.query_name][1])
                elif read.is_paired and not read.is_unmapped  and not read.is_secondary:
                    pairedReadIDs[read.query_name]=1
                    pairedReadStartEnd[read.query_name]=( min(readStart,readEnd), max(readStart,readEnd) )
                    pairedReadSequences[read.query_name]=read.query_alignment_sequence.replace("N","-")

            if not read.is_paired or (read.query_name in pairedReadIDs and pairedReadIDs[read.query_name]==2):
                if bam.get_reference_name(read.reference_id) not in sampleChrRegionReads[sample]:
                    continue #sometimes spurrious hits are generated on non-target chromosomes
                if (read.reference_id, readStart, readEnd) not in readStartEndToRegion:
                    #determine which region the read falls into
                    for region in sampleChrRegionReads[sample][bam.get_reference_name(read.reference_id)].keys():
                        if abs(readStart-region[0])<config.maxDistanceBetweenReads and abs(readEnd-region[1])<config.maxDistanceBetweenReads:
                            readStartEndToRegion[(read.reference_id,readStart, readEnd)]=region
                            break
                    if (read.reference_id, readStart, readEnd) not in readStartEndToRegion: #this checks if read in in target regions, if not it is skipped
                        continue
                if (read.query_name in pairedReadIDs and pairedReadIDs[read.query_name]==2):
                    sampleChrRegionReads[sample][bam.get_reference_name(read.reference_id)][readStartEndToRegion[(read.reference_id,readStart, readEnd)]].append(pairedReadSequences[read.query_name].replace("N","-"))
                else:
                    sampleChrRegionReads[sample][bam.get_reference_name(read.reference_id)][readStartEndToRegion[(read.reference_id,readStart, readEnd)]].append(read.query_alignment_sequence.replace("N","-"))


    #check if for some sampels reads need to be downsampled
    for sample in sampleChrRegionReads:
        for chr in sampleChrRegionReads[sample]:
            for region in sampleChrRegionReads[sample][chr]:
                if len(sampleChrRegionReads[sample][chr][region])>config.MaxReadsToOutput:
                    subset=random.sample(range(0, len(sampleChrRegionReads[sample][chr][region]) ), config.MaxReadsToOutput)
                    sampleChrRegionReads[sample][chr][region]=[sampleChrRegionReads[sample][chr][region][i] for i in subset]


    return sampleChrRegionReads
    #the reads from samples have been collected into sampleChrRegionReads[sample][chr][region]=[read seqeuences]
    #print them to fasta and allign using mafft

# from collections import Counter
# import configs.LeenNewBarcodes as Config
# config=Config.configData()
# intervalDepths={'Pf3D7_03_v3': {(221325, 221660): 26698, (221326, 221353): 1529, (221639, 221663): 14, (221326, 221400): 24}, 'Pf3D7_13_v3': {(1465038, 1465398): 12553}, 'Pf3D7_07_v3': {(403511, 404010): 998, (403512, 403553): 48}}
# intervalDepths={'Pf3D7_13_v3': {(1465038, 1465398): 12553}}
# validsamples=["sample2.bam"]
# getReadsFromBams(validsamples,intervalDepths,config)