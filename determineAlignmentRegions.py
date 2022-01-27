from numpy import append
import pysam as ps
from os import chdir, listdir
from os.path import isfile, splitext, join
from matplotlib import pyplot as plt
from collections import Counter

#these are 1-based array indices, python uses 0-based
#1-QNNAME
#2-FLAG
#3-Reference sequence name 
#4-1-based leftmost mapping (POS)
#5-MARQ, 255=quality not available



def selectNonOverlappingIntervals(intervalsCount, config):
    aggregatedIntervalDepths={}#for each removed interval, it's depth is inhertied by the interval that knocked it out
    #this works for single read sequencing, not sure about paired-end
    multiReadIntervals=[]
    for interval in intervalsCount.keys():
        if intervalsCount[interval]>1:
            #dismiss any interval with fewer than 2 reads mapping to it.                
            multiReadIntervals.append(interval)
            aggregatedIntervalDepths[interval]=intervalsCount[interval]
    #do an elimination loop for all intervals@
    count=0
    i=0
    while i<len(multiReadIntervals) and count<10000: #the count limit is just-in-case circuit breaker
        #actual selection of non-overlapping intervals.
        base_interval_start=multiReadIntervals[i][0]
        base_interval_end=multiReadIntervals[i][1]
        k=i+1
        while k<len(multiReadIntervals):
            evaluated_interval_start=multiReadIntervals[k][0]
            evaluated_interval_end=multiReadIntervals[k][1]
            #check if intervals overlap, there are four cases, but use simple difference between starts and ends.
            if abs(base_interval_start-evaluated_interval_start)<config.maxDistanceBetweenReads and abs(base_interval_end-evaluated_interval_end)<config.maxDistanceBetweenReads:
                k=k-1 #because element is removed from iterated list, the counter should not change for this loop
                #intervals overlap, keep the most frequent one
                if abs(base_interval_start-base_interval_end) < (evaluated_interval_start-evaluated_interval_end):
                    multiReadIntervals.remove( (base_interval_start, base_interval_start) )
                    aggregatedIntervalDepths[(evaluated_interval_start,evaluated_interval_end)]+=aggregatedIntervalDepths[ (base_interval_start, base_interval_end) ]
                    i=i-1 #because base interval is removed, the next interval is same index as current one.  
                    break
                    #the "base_interval_values" can be reset, but I don't think this is necessary because the alorithm is not based on exact matches
                else:
                    multiReadIntervals.remove( (evaluated_interval_start,evaluated_interval_end) )
                    aggregatedIntervalDepths[ (base_interval_start, base_interval_end) ]+=aggregatedIntervalDepths[(evaluated_interval_start,evaluated_interval_end)]
            k+=1
            count+=1
        i+=1
    keys=set(aggregatedIntervalDepths.keys())
    for interval in keys:
        if interval not in multiReadIntervals or aggregatedIntervalDepths[interval]<config.minRegionDepth: #fewer than N reads map to the region
            aggregatedIntervalDepths.pop(interval)
    return multiReadIntervals, aggregatedIntervalDepths


chrRegions={} #key=chr, value=[start, end]
pairedReadIDs={} #key=read id, tuple=min(read1, read2), max(read1,read2)
def run(samples, config):
    for sample in samples:
        bam = ps.AlignmentFile(config.bamDir+sample, "rb",check_sq=False)
        for read in bam.fetch():
            if not read.is_unmapped:
                if bam.get_reference_name(read.reference_id) not in chrRegions:
                    chrRegions[bam.get_reference_name(read.reference_id)]=[]
                if read.is_paired  and not read.is_unmapped and not  read.is_secondary:# and read.is_proper_pair:
                    ###Illumina data
                    #print(read.query_name)
                    if read.query_name not in pairedReadIDs:
                        pairedReadIDs[read.query_name]=( min(read.reference_start,read.reference_end), max(read.reference_start,read.reference_end) )
                    else:
                        #this is the second read of a pair
                        #using min of three numbers allows to disregard the mapping orentation.
                        startOfPair=min(read.reference_start,read.reference_end, pairedReadIDs[read.query_name][0])
                        endOfPair=max(read.reference_start,read.reference_end, pairedReadIDs[read.query_name][1])
                        chrRegions[bam.get_reference_name(read.reference_id)].append( (startOfPair,endOfPair) ) 
                else:
                    #single end sequencing or 3rd gen sequencing
                    chrRegions[bam.get_reference_name(read.reference_id)].append( (read.reference_start,read.reference_end) ) 

    joinedIntervals={} #key=chr, value=[start,end] of collapsed overlapping intervals    
    joinedDepths={} #key=chr, value={interval: depth} of collapsed overlapping intervals   
    for chr in chrRegions.keys():
        joinedIntervals[chr]=[] #value=[start,end] of collapsed overlapping intervals
        #iterate over the regions to collapse them into overlapping regions. 
        intervalsCount=Counter(chrRegions[chr])
        joinedIntervals[chr],joinedDepths[chr]=selectNonOverlappingIntervals(intervalsCount, config)
    #remove chromosomes without intervals - this is likely due to very low coverage of some regions
    for chr in joinedDepths:
        if len(joinedDepths[chr])==0:
            joinedIntervals.pop(chr)
    return joinedIntervals, joinedDepths

