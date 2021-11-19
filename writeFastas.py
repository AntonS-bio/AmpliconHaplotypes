import subprocess
from os import mkdir
from os.path import exists
import shutil
import time

def run(sampleChrRegionReads, fileContent, config):

    tempDir=config.tempDir

    if exists(tempDir):
        shutil.rmtree(tempDir)
    mkdir(tempDir)
    alignedFiles={} #key=sample, value={key=chr, value={key=region, value=aligned filename}}
    filesToAlign=set()
    for sample in sampleChrRegionReads:
        alignedFiles[sample]={}
        for chr in sampleChrRegionReads[sample]:
            alignedFiles[sample][chr]={}
            for region in sampleChrRegionReads[sample][chr]:
                if fileContent=="singleSample":
                    filePrefix=f'{sample}{chr}_{region[0]}_{region[1]}'
                    writetype="w"
                else:
                    #each fasta contans all reads from all samples for a given region
                    filePrefix=f'{chr}_{region[0]}_{region[1]}'
                    writetype="a"
                alignedFiles[sample][chr][region]=f'mafft{filePrefix}'
                filesToAlign.add(filePrefix)
                with open(tempDir+f'{filePrefix}.fasta',writetype) as mafftinput:
                    i=1
                    for read in sampleChrRegionReads[sample][chr][region]:
                        mafftinput.write(f'>{sample} read_{str(i)}\n')
                        mafftinput.write(read+"\n")
                        i+=1
    
    runTimes={} #keeps track of run time, mostly mafft, per region. Long run indicates likely problems
    for filePrefix in filesToAlign:
        
        startTime = time.time()
        print("processing "+filePrefix)        
        subprocess.run(f'mafft --quiet {tempDir}/{filePrefix}.fasta > {tempDir}/mafft{filePrefix}.fasta', shell=True, executable="/bin/bash", stdout=subprocess.PIPE)
        runTimes[filePrefix]=startTime-time.time()
    return alignedFiles, runTimes