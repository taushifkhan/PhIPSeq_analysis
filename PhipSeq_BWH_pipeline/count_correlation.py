#/gpfs/software/genomics/PhipSeq/env python3.6
import pandas as pd
import glob
import json
import numpy as np
import gzip
import scipy.stats as ST


def count_stat(sample,tr1,tr2):
    """
    count files
    tr1 : technical repeat 1
    tr2 : technical repeat 2
    """
    #import ipdb; ipdb.set_trace()
    try:
        tr1_file = np.loadtxt(gzip.GzipFile(tr1,'r'),skiprows=1,usecols=1,dtype=int,delimiter=",")
        tr2_file = np.loadtxt(gzip.GzipFile(tr2,'r'),skiprows=1,usecols=1,dtype=int,delimiter=",")
        corr = ST.pearsonr(tr1_file,tr2_file)
        return [sample,np.sum(tr1_file),np.sum(tr2_file),corr[0],corr[1],np.median(tr1_file),np.median(tr2_file)\
,np.max(tr1_file),np.max(tr2_file),\
np.mean(tr1_file),np.mean(tr2_file)]
    except:
        print ("[ERROR] file not found!! ")
        return []


#write to file
def count_statistics(countDir,sampleIds,countFile="count_stat.txt",errorFile="count_error.txt"):
    stat_fl = open(countFile,'w')
    stat_fl.write("SampleID,TR1_sum,TR2_sum,pearson_coff,pearson_pval,TR1_median,TR2_median,TR1_max,TR2_max,TR1_mean,TR2_mean\n")
    error_fl = open('count_error.txt','w')

    for k in sampleIds:
        _samples= [countDir+k+"_TR1.count.gz",countDir+k+"_TR2.count.gz"]
        print (k, _samples)
        res = count_stat(k,_samples[0],_samples[1])
        if len(res):
            res = [str(i) for i in res]
            stat_fl.write("{}\n".format(",".join(res)))
        else:
            error_fl.write("{}\n".format(k))
    stat_fl.close()
    error_fl.close()



if __name__ =="__main__":
    import sys
    sm = pd.read_csv(sys.argv[1])
    print ("Number of samples:",sm.sid.unique().shape)
    count_statistics('alignmentCount/',sm.sid.unique())
