#/gpfs/software/genomics/PhipSeq/env python3.6
import pandas as pd
import glob
import json
import gzip
# declare paths and scripts

data_path = {
        'metadata':'SLED_sample_map_filtered.csv',
        'phipDir':'phipDir/',
        'summDir':'resStatistics/',
        'countDir':'alignmentCount/',
                }

def merge(dirarray,outFile,index="id",compressed=0):
    """
    __doc_string__ : merge the files
    @param
    dirarray : list of files
    outfile : merged out file
    index : except for species nad oranism all indexed as "id"

    """
    k = pd.DataFrame()
    try:
        listFile = tuple(dirarray)[0]
        concatfile = tuple(outFile)[0]
        for i in listFile:
            if k.shape[1]==0:
                if compressed:
                    k = pd.read_csv(gzip.open(i)).set_index(index)
                else:
                    k = pd.read_csv(i).set_index(index)
            else:
                if compressed:
                    _ = pd.read_csv(gzip.open(i)).set_index(index)
                else:
                    _ = pd.read_csv(i).set_index(index)

                k = k.join(_)
        print (k.shape,)
        k.to_csv(concatfile)
        print (listFile,concatfile)
        return 1
    except:
        import ipdb; ipdb.set_trace();

data_meta = pd.read_csv(data_path['metadata'])
unqSample = data_meta.sid.unique()
print ("Number of samples to be analyzed : {}\nUnique samples: {}".format(data_meta.shape[0],len(unqSample)))

_sf = {}
for _name in unqSample:
    _sf[_name] = data_path["phipDir"]+_name+".zhit.csv"

### pipeline starts here

rule all:
    input:
        data_path['summDir']+"count_peptide.csv",
        data_path['summDir']+"pval_peptide.csv",
        data_path['summDir']+"hit_peptide.csv",

rule compress_hit:
    input:lambda wildcards: _sf[wildcards.name]
    output:data_path['phipDir']+"{name}.zhit.csv.gz"
    shell:"gzip -c {input} > {output}"


rule hit_merge:
    input:
        cfile=expand("{dataDir}{name}.count.gz",dataDir=data_path['countDir'],name=data_meta.name_as.values),
        pfile=expand("{dataDir}{name}.zpval.csv",dataDir=data_path['phipDir'],name=data_meta.name_as.values),
        hfile=expand("{dataDir}{name}.zhit.csv",dataDir=data_path['phipDir'],name=list(unqSample)),
    output:
        cfile=data_path['summDir']+"count_peptide.csv",
        pfile=data_path['summDir']+"pval_peptide.csv",
        hfile=data_path['summDir']+"hit_peptide.csv",
    run:
        merge({input.cfile},{output.cfile},index="id",compressed=1)
        merge({input.hfile},{output.hfile},index="id")
        merge({input.pfile},{output.pfile},index="id")
