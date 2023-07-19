import glob
import pandas as pd
import json

data_path = {
        'pathROOT':'/gpfs/projects/SDR_PPM1-1220-150017_nmarr/analysis/SidraOccHealth_vir3/run41_NS500_HKYTMBGX5/fastq_2018-XX-XX/',
        'rundir':['QBBrun32'],        
        'metadata':'sample_metadata.csv',
                }

data_meta = pd.read_csv(data_path['metadata'])
print ("Number of samples to be analyzed : {}".format(data_meta.shape[0]))

def searchPool(sampleId,rootpath=data_path['pathROOT'],pdirs=data_path['rundir']):
    sample_pool = []
    dir_host = []
    #try:
    for k in pdirs:
       _path =rootpath+k+'/fastq/{}*.fastq'.format(sampleId)
       for i in glob.glob(_path):
           sample_pool.append(i)
           dir_host.append(k)
       if len(sample_pool):
           _sd = list(set(dir_host))
           if len(_sd)==1:
               return sample_pool,1,list(set(dir_host))[0]
           else:
               return sample_pool, len(_sd),','.join(_sd)
    print ("error")
    return sample_pool,0,dir_host


def sample_type(sampleName):
    if sampleName[:3] == 'SIR':
        return 'QBB'
    elif sampleName[0] == "B":
        return "Bead"
    elif sampleName[:3] == "OCC":
        return "OCC"
    else:
        return sampleName
    

sample_fastqs = {}
sample_dir = {}
file_error = []
fERR = open('errorFile.txt','w')
allDir_info = []
dir_summ = []


for k in data_meta.runId:
    _tr = data_meta[data_meta.runId==k].name_as.values[0]
    sampleId = data_meta[data_meta.runId==k].sid.values[0]
    print (k,_tr,sampleId)
    _sampletype = sample_type(sampleId)
    _d = data_meta[data_meta.runId==k].name_as.values[0]
    _pool,len_dir,_dir  = searchPool(k)
    if len(_pool):
        sample_fastqs[_d] = _pool
        sample_dir[_d] = _dir
        allDir_info.append([sampleId,_tr,k,_pool,_dir,_sampletype])
        dir_summ.append([sampleId,_tr,k,len(_pool),len_dir,_dir,_sampletype])
        if len(_pool) > 5:
            print (_d)
            import ipdb; ipdb.set_trace();
    else:
        file_error.extend(k)
        fERR.write("%s\n"%k)
fERR.close()

k = pd.DataFrame(allDir_info,columns=['sid','name_as','runId','fastqFiles','run_dirName','sampleType'])
k2 = pd.DataFrame(dir_summ,columns=['sid','name_as','runId','#fastqs','#dirs','run_dirName','sampleType'])
print (k.head())
k.to_csv('../sample_metaData_fastq.csv',index=False)
k2.to_csv('./sample_metaData_summary.csv',index=False)
print ("Sample File error: {}\nIds in file ".format(len(file_error)))
unqSample = list(data_meta.sid.unique())
print ("NUmber of Unique samples: {} Analyzing: {}".format(len(unqSample),len(list(sample_fastqs.keys()))))  
