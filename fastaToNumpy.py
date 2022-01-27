import numpy as np
import pandas as pd
from Bio import SeqIO

def run(MSAfiles, config):
    tempDir=config.tempDir
    for filePrefix in MSAfiles:
        print(filePrefix)
        #check the length and depth of MSA to pre-allocate numpy matrix
        MSAreads=0
        MSAreadlength=0
        for record in SeqIO.parse(f'{tempDir}/{filePrefix}.fasta', "fasta"):
            MSAreads+=1
            if MSAreadlength==0:
                MSAreadlength=len(str(record.seq))

        #generate numpy matrix
        symbols=5 # ACTG and -
        matrix=np.zeros( ( MSAreads, (MSAreadlength) ), dtype="b")
        samplesVector=[] #this will be a column in output matrix
        readCounter=0

        for record in SeqIO.parse(f'{tempDir}/{filePrefix}.fasta', "fasta"):
            samplesVector.append(record.id.split(" ")[0].replace(">",""))            
            for i in range( len(record.seq) ):
                if record.seq[i]=="a":
                    matrix[ readCounter , i ]=1
                elif record.seq[i]=="t":
                    #matrix[ readCounter , i+MSAreadlength ]=1
                    matrix[ readCounter , i ]=2
                elif record.seq[i]=="c":
                    #matrix[ readCounter , i+MSAreadlength*2 ]=1
                    matrix[ readCounter , i ]=3
                elif record.seq[i]=="g":
                    #matrix[ readCounter , i+MSAreadlength*3 ]=1
                    matrix[ readCounter , i ]=4
                elif record.seq[i]=="-":
                    #matrix[ readCounter , i+MSAreadlength*4 ]=1
                    matrix[ readCounter , i ]=5
                else:
                    print("Problem")
            readCounter+=1

        #drop the columns that have fewer than 1% and 5 bases ##This is no longer required because instead of column per possible base, each columns is 1-5 for bases
        #this work with sklearn definition of hamming distances
        Positions=["Pos_"+str(f) for f in range(0,MSAreadlength) ]
        #As=["A"+str(f) for f in range(0,MSAreadlength) ]
        #Ts=["T"+str(f) for f in range(0,MSAreadlength) ]
        #Gs=["G"+str(f) for f in range(0,MSAreadlength) ]
        #Cs=["C"+str(f) for f in range(0,MSAreadlength) ]
        #Blanks=["Blank"+str(f) for f in range(0,MSAreadlength) ]
        #data=pd.DataFrame(matrix, columns=[As+Ts+Cs+Gs+Blanks]) #converting numpy to 
        data=pd.DataFrame(matrix, columns=Positions) #converting numpy to 
        #data=data[data.columns[ (data.sum()>config.minPosDepthAbsolute) | (data.sum()>int(MSAreads*config.minPosDepthRelative) ) ]]
        data["Sample"]=samplesVector
        data.to_pickle(f'{tempDir}/{filePrefix}.pkl')
        #print(data.shape)


