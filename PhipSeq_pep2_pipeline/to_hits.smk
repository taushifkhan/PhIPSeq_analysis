#/gpfs/software/genomics/PhipSeq/env python3.6
import pandas as pd
import glob
import json
import gzip
# declare paths and scripts

data_path = {
        'metadata':'SLED_sample_map_filtered.csv',  # metadata csv file see readme
        'mergeDir':'fastQmerged/',	# path to store merged fastq files
        'countDir':'alignmentCount/',	# path to store alignmed bam and count files
        'phipDir':'phipDir/',		# path to store enrichment results
        'summDir':'resStatistics/',	
        'reference_pep':'refSeq/peptidome', # directory has reference sequence unzip the directories in refSeq and provide the path
        'input':'avg_sled_inputlib.csv.gz', # generate the input count for after the alignment and use this for enrichment
                }

bwh_path={
		'pval':'scripts/calc_zipval.py',
        'zhit':'scripts/call_hits.py',
        'zscore':'scripts/calc_scores.py'
}

data_meta = pd.read_csv(data_path['metadata'])
unqSample = data_meta.sid.unique()
print ("Number of samples to be analyzed : {}\nUnique samples: {}".format(data_meta.shape[0],len(unqSample)))

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


#write to file
def count_statistics(countDir,sampleIds,countFile="count_stat.txt",errorFile="count_error.txt"):
    stat_fl = open(countFile,'w')
    stat_fl.write("SampleID,TR1_sum,TR2_sum,pearson_coff,pearson_pval,TR1_median,TR2_median,TR1_max,TR2_max,TR1_mean,TR2_mean\n")
    error_fl = open('count_error.txt','w')

    for k in sampleIds:
        count_samples[k]=[coutDir+k+"_TR1.count.gz",data_path['countDir']+k+"_TR2.count.gz"]
        print (k, count_samples[k])
        res = count_stat(k,count_samples[k][0],count_samples[k][1])
        if len(res):
            res = [str(i) for i in res]
            stat_fl.write("{}\n".format(",".join(res)))
        else:
            error_fl.write("{}\n".format(k))

    print ("Number of error samples: ",len(error_sample))
    stat_fl.close()
    error_fl.close()


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

### pipeline starts here

rule all:
    input:
        expand("{dataset}{sid}.zhit.csv",dataset=data_path['phipDir'],sid=unqSample)

_sf = {}
for _name in data_meta.name_as.values:
    _sf[_name] = data_path["mergeDir"]+_name+'.fastq.gz'
    #if _name[-2:] == "R1":
    #    _sf[_name] = data_path["mergeDir"]+_name[:-3]+'-TR1.fastq.gz'
    #else:
    #    _sf[_name] = data_path["mergeDir"]+_name[:-3]+'-TR2.fastq.gz'

rule bow_align:
    input:
        fq=lambda wildcards: _sf[wildcards.name],
        ref=data_path['reference_pep'],
    output:
        sam=data_path['countDir']+"{name}.sam"
    shell:
        """
        module load Bowtie2/2.3.4.3
        bowtie2 --very-sensitive-local -x {input.ref}pep2 {input.fq} -S {output.sam}
        """

rule compress_sam:
    input:data_path['countDir']+"{name}.sam",
    output:data_path['countDir']+"{name}.sam.gz"
    shell:
        """
        gzip {input}
        du -sh {output}
        """
rule bam_sort:
    input:
        sam=data_path['countDir']+"{name}.sam.gz"
    output:
        tmp_bam=temp(data_path['countDir']+"{name}.bam"),
        tmp_sorted_bam=temp(data_path['countDir']+"{name}.sorted.bam")
    shell:
        """
        module load SAMtools/1.9
        samtools view -Sb {input.sam} > {output.tmp_bam}
        samtools sort {output.tmp_bam}> {output.tmp_sorted_bam}
        """
rule gen_alignment:
    """
    ** for mac the last sed rule
    samtools idxstats {output.sorted_bam} | cut -f 1,3 |sed -e '/^\*/d' -e '1i\'$'\n'id,file''|tr "\\t" ","> {output.countfl}
    ** for linux
    samtools idxstats {output.sorted_bam} | cut -f 1,3 | sed -e '/^\*\t/d' -e "1 i id\file" | tr "\\t" "," > {output.countfl}
    	samtools idxstats {output.sorted_bam} | cut -f 1,3 | sed -e '/^\*\t/d' -e "1 i id\{input.fq}" | tr "\\t" "," > {output.countfl}
    """
    input:
        bam=data_path['countDir']+"{name}.sorted.bam"
    output:
        tmp_count=temp(data_path['countDir']+"{name}.count"),
    shell:
        """
        module load SAMtools/1.9
        samtools idxstats {input.bam} | cut -f 1,3 |sed -e '/^\*/d' -e '1i\\'$'\\n'id,{wildcards.name}''|tr "\\t" ","> {output.tmp_count}
        """

rule gen_compressedcount:
    input:data_path['countDir']+"{name}.count"
    output:data_path['countDir']+"{name}.count.gz"
    shell:"gzip {input}"
#expand("{outdir}{name}.count.gz",outdir=data_path['mergeDir'],name=list(sample_fastqs.keys()))
count_gz={}
for _sample in data_meta.sid.unique():
    tr = data_meta[data_meta.sid==_sample].name_as.values
    for k in tr:
        count_gz[k]=data_path['countDir']+k+".count.gz"

rule count_stat:
    output:"count_stat.txt"
    run:
       count_statistics(data_path['countDir'],data_meta.sid.unique(),countFile={output})


rule bwh_pvalue:
    input:
        lambda wildcards: count_gz[wildcards.name]
        #fq=data_path['countDir']+"{name}.count.gz"
    output:
        pval=data_path['phipDir']+"{name}.zpval.csv",
    log:"log/pval/{name}/"
    shell:
        """
        module load python/3.4.9
        python {bwh_path[pval]} {input} {data_path[input]} {log} > {output.pval}
        """

rule compress_pval:
    input:data_path['phipDir']+"{name}.zpval.csv"
    output:data_path['phipDir']+"{name}.zpval.csv.gz"
    shell:"gzip -c {input} > {output}"

tr_samples = {}
for k in unqSample:
    tr_samples[k]=[data_path['phipDir']+k+"-TR1.zpval.csv.gz",data_path['phipDir']+k+"-TR2.zpval.csv.gz"]

rule bwh_zhit:
    input:
        lambda wildcards: tr_samples[wildcards.name]
    output:
        csv = data_path['phipDir']+"{name}.zhit.csv",
    log:"log/zhit/{name}/"
    params:
        threshold=2.3
    shell:
        """
        module load python/3.4.9
        python {bwh_path[zhit]} {input} {log} {params.threshold} > {output.csv}
        """
