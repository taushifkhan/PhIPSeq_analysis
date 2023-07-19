#/gpfs/software/genomics/PhipSeq/env python3.6
import pandas as pd
import glob
import json
import gzip
import count_correlation as cCorr
# declare paths and scripts

data_path = {
        'metadata':'sample_metaData_fastq.csv',
        'mergeDir':'fastQmerged/',
        'countDir':'alignmentCount/',
        'phipDir':'phipDir/',
        'summDir':'resStatistics/',
        'reference_vir3':'/gpfs/projects/tmedicine/tkhan/bwh/reference_seq/',
        'input':'/gpfs/projects/SDR_PPM1-1220-150017_nmarr/analysis/SidraOccHealth_vir3/run41_NS500_HKYTMBGX5/fastq_2018-XX-XX/input/input.csv.gz',
        'cleanVir3':'/gpfs/projects/tmedicine/tkhan/bwh/metadata/VIR3_clean.csv.gz',
        'nhits_sample':'/gpfs/projects/tmedicine/tkhan/virScan_qbb/nhits/nhits_samples_qbb.csv.gz',
        'nhits_beads':'/gpfs/projects/tmedicine/tkhan/virScan_qbb/nhits/nhits_beads.csv.gz',
        'errorFile':'errorIds_sample.txt',
        'fastqPath':'/gpfs/projects/tmedicine/tkhan/virScan_qbb/fastq_path_QBB.json',
                }

bwh_path={
		'pval':'/gpfs/projects/tmedicine/tkhan/bwh/scripts/calc_zipval.py',
        'zhit':'/gpfs/projects/tmedicine/tkhan/bwh/scripts/call_hits.py',
        'zscore':'/gpfs/projects/tmedicine/tkhan/bwh/scripts/calc_scores.py'
}

data_meta = pd.read_csv(data_path['metadata'])
print ("Number of samples to be analyzed : {}".format(data_meta.shape[0]))

unqSample = data_meta.sid.unique()

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
        data_path['summDir']+"count_peptide.csv",
        data_path['summDir']+"pval_peptide.csv",
        data_path['summDir']+"hit_peptide.csv",
        data_path['summDir']+"spp_score.csv",
        data_path['summDir']+"org_score.csv"


_sf = {}
for _name in data_meta.name_as.values:
    _t = data_meta[data_meta.name_as==_name].fastqFiles.values[0][2:-2].strip().split(",")
    _sf[_name] = [i.strip().strip("'") for i in _t]

print (_sf)
rule merge_fastq:
    input:
        lambda wildcards: _sf[wildcards.name]
    output:
        m1=data_path['mergeDir']+"{name}.fastq",
    shell:
        """
        cat {input} > {output.m1}
        du -sh {output.m1}
        """
#import ipdb; ipdb.set_trace();

rule compressed:
    input:data_path['mergeDir']+"{name}.fastq"
    output:data_path['mergeDir']+"{name}.fastq.gz"
    shell:"gzip {input}"

rule bow_align:
    input:
        fq=data_path['mergeDir']+"{name}.fastq.gz",
        ref=data_path['reference_vir3']
    output:
        sam=data_path['countDir']+"{name}.sam"
    shell:
        """
        module load Bowtie2/2.3.4.3
        bowtie2 --very-sensitive-local -x {input.ref}vir3 {input.fq} -S {output.sam}
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
        module load phipseq-env/24Jan19
        python {bwh_path[pval]} {input} {data_path[input]} {log} > {output.pval}
        """

rule compress_pval:
    input:data_path['phipDir']+"{name}.zpval.csv"
    output:data_path['phipDir']+"{name}.zpval.csv.gz"
    shell:"gzip -c {input} > {output}"

tr_samples = {}
for k in unqSample:
    tr_samples[k]=[data_path['phipDir']+k+"_TR1.zpval.csv.gz",data_path['phipDir']+k+"_TR2.zpval.csv.gz"]

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
        module load phipseq-env/24Jan19
        python {bwh_path[zhit]} {input} {log} {params.threshold} > {output.csv}
        """


rule compress_hit:
    input:data_path['phipDir']+"{name}.zhit.csv"
    output:data_path['phipDir']+"{name}.zhit.csv.gz"
    shell:"gzip -c {input} > {output}"

rule bwh_zscore:
    input:
        hitfile=data_path['phipDir']+"{name}.zhit.csv.gz"
    output:
        sp=data_path['phipDir']+"{name}.spp.zscore.csv",
        org=data_path['phipDir']+"{name}.org.zscore.csv"
    shell:
        """
        module load phipseq-env/24Jan19
        python {bwh_path[zscore]} {input.hitfile} {data_path[cleanVir3]} {data_path[nhits_beads]} \
        {data_path[nhits_sample]} Species 7 > {output.sp}

        python {bwh_path[zscore]} {input.hitfile} {data_path[cleanVir3]} {data_path[nhits_beads]} \
        {data_path[nhits_sample]} Organism 7 > {output.org}
        """

# merging all and compressing

#rule compress_pval:
#    input:data_path['phipDir']+"{name}.zpval.csv"
#    output:data_path['phipDir']+"{name}.zpval.csv.gz"
#    shell:"gzip {input}"

#rule compress_hit:
#    input:data_path['phipDir']+"{name}.zhit.csv"
#    output:data_path['phipDir']+"{name}.zhit.csv.gz"
#    shell:"gzip {input}"



rule hit_merge:
    input:
        cfile=expand("{dataDir}{name}.count.gz",dataDir=data_path['countDir'],name=data_meta.name_as.values),
        pfile=expand("{dataDir}{name}.zpval.csv",dataDir=data_path['phipDir'],name=data_meta.name_as.values),
        hfile=expand("{dataDir}{name}.zhit.csv",dataDir=data_path['phipDir'],name=list(unqSample)),
        sfile=expand("{dataDir}{name}.spp.zscore.csv",dataDir=data_path['phipDir'],name=list(unqSample)),
        ofile=expand("{dataDir}{name}.org.zscore.csv",dataDir=data_path['phipDir'],name=list(unqSample))
    output:
        cfile=data_path['summDir']+"count_peptide.csv",
        pfile=data_path['summDir']+"pval_peptide.csv",
        hfile=data_path['summDir']+"hit_peptide.csv",
        sfile=data_path['summDir']+"spp_score.csv",
        ofile=data_path['summDir']+"org_score.csv"
    run:
        merge({input.cfile},{output.cfile},index="id",compressed=1)
        merge({input.hfile},{output.hfile},index="id")
        merge({input.pfile},{output.pfile},index="id")
        merge({input.sfile},{output.sfile},index="Species")
mpo        merge({input.ofile},{output.ofile},index="Organism")
