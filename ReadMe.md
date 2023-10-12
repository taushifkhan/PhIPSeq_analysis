### Manual for PhIP-Seq Bioinformatics Analysis

###
@author : taushif khan
@email  : tkhan@sidra.org | taushif.khan@jax.org
@Date   : 21 Aug 2019
###

### About 
Given a sequencing output and metadata for sampleIds and fastq files, the pipeline
propose to provide stepwise analysis and save results in systematic way. The

### How to run

copy the directory PhipSeq_BWH_pipeline/ to location of your choice.
from terminal. go to the pipeline directory.
```
cd PhipSeq_BWH_pipeline/
module load env-python/3.6.6
```
1. Prepare sample metadata
the program will take csv file as input. with following fields as columns (at-least),

```
<sid>,<name_as>,<runId>

sid: sample Identity <unique name of a sample in study>
name_as : clone (or technical repeats) of the sample all the fastq files will be merged on this name.
runId : clone name with primer identifier for identification file to be merged.
```
save the file as csv and proceed to next step.

2. Add data to sample metadata [go to : /DataFile_prePro]
  ```
  cd DataFile_prePro/
  ```
  a. setup Correct Path: open fasqPaths.py from terminal and edit
  assuming the fastqs are in differnt directory of a root directory
  data_path{
     'pathROOT': 'full path of root directory',
     'rundir': ['dir1','dir2',...],
     'metadata': 'a csv file with required columns'}

  with all set up run,
  ```
  python fastqPaths.py
  ```
  The program will search presence of files in give location using the 'runId' and
  add a column "fastqFiles". Additional information on Directory name and sample type my be added.
  Output will be written in ../sample_metaData_fastq.csv
  if a sampleId do not find any fastq files those will be written in errorFile.txt

3. Run snakemake

after successful execution of previous step you will have metadat file in parent
snakemake directory : sample_metaData_fastq.csv

This file has all fastq files for different sample ids. Now its time to run the pipeline.

! Before running make sure following by reading data_path in the *.smk file*:
A. data_path
1. metadata: <we have generated in last step and should be in the same directory>
2. outputDirectories : mergeDir, countDir , phipDir and summDir are output Directories, where results at different stages are stored. Can changed to user required name.
3. referenceDirectory : refence_vir3
4. cleanVir3 : metadata with annotation of protein and species
5. nhits_samples: important see section : "Generate nhits for VIRAL EXPOSURE"
6. nhits_beads : (5)
B. bwh_path : path for python codes from BWH.

4. Launch Snakemake:
make a dry run to make sure files are present
```
snakemake -np -s build_phipSeq.smk
```

Better to initiate a screen and run snake make with cluster specification
```
screen
snakemake -p -s build_phipSeq.smk -j --cluster 'bsub -n 20 -J BWH -P takhBWH -e bwh_test_vir3.err -o bwh_test_vir3.out'

<get out of screen>
shift + a + d
```

we use 20 nodes with error and output file of the run on specific location.

### Generate nhits for VIRAL EXPOSURE

```
cd pythonModules/
python get_nhits.py
cp *.gz ../nhits/
```

to evaluate the virus exposure, following peptides are dropped based on background and samples,
1. peptides that are present in atleast 3 mockIp.negative controls
2. peptides that are not present in atleast 2 samples.

from science paper [2015]: eliminated peptides that come at least 3 of 22 immunuprecipitations with beads alone . We also filtered out any peptide that were not enriched in at least two of the sample"


### Requirements:

The workflow is written in snakemake. Snakemake uses python for providing input \
output flow to number of different programs that might be used in the whole pipeline.
For PhIP-seq analysis primarily we will use,
a. Aligner : Bowtie2
b. SamsTools
c. python

1. To run on HPC:
    '''
    module load env-python/3.6.6
    '''
    this will load all the required packages and environment to run this workflow.

2. To set-up a local run :contact me

### Preparing metadata

This is most important for this pipeline to run.

### Steps:
The pipeline runs with different sets of instructions executed in top to bottom.
These sets are called "rules". Each rule describes a step in PhIP-seq analysis with
exclusive set of input and output.

The workflow checks availability of defined input in rule "all". As per the
presence of different files, and executes following steps:

Job counts:
	count	jobs
	1	all
	192	bam_sort
	192	bow_align
	192	bwh_pvalue
	96	bwh_zhit
	96	bwh_zscore
	96	compress_hit
	192	compress_pval
	192	compress_sam
	192	compressed
	192	gen_alignment
	192	gen_compressedcount
	1	hit_merge
	192	merge_fastq
	2018
