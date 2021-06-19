# unmappedReadsAnalysis
- Andrew Davidson
- aedavids@ucsc.edu
- 5/1/21
- BME 237 Applied RNA Bioinformatics Class project

cheat sheet how to run batch jobs
```
setsid sh -c ' set -x;  runSTAR.sh ' >  runSTAR.sh.out.`~/extraCellularRNA/bin/dateStamp.sh` 2>&1 &

ps -e -o pid,ppid,pgid,command,user |head -n 1; ps -e -o pid,ppid,pgid,command,user |grep aedavids

pstree aedavids

kill -TERM -[pgid]
```
Claim is we have low mapping rates

## explore kras.ipsc data set
in vitro so should be easier to identify techincal issues. i.e. unlikely bacterial contamination

1. bin/mineSalmonLogs.sh  > data/salmonMappingRate.tsv
   - produces tsv meta data file 
   - sample name, ref, right fastq, left fastq, mapping rate
   - samples where run using different references.
   - make sure we are comparing apples to apples
2. sample with the lowest mapping rate
   ```
   36.8215
       kras.ipsc/data/exo.data/gen1c.day.5.exo.input.data/kras.3/quantFiles		
       /public/groups/kimlab/gen.v31.k.31.v.0.14.index
       
   37.5851
       kras.ipsc/data/exo.data/gen1c.day.5.exo.input.data/ctrl.1/quantFiles
	   /public/groups/kimlab/gen.v31.k.31.v.0.14.index
   ```

3. samples with highest mapping rate

    ```
    95.1761 
        kras.ipsc/data/bulk.data/day.7/kras.2/gencode.salmon.out
        /public/groups/kimlab/indexes/gen.32.ucsc.rmsk.index
        
    95.1761 
        kras.ipsc/data/bulk.data/day.7/kras.2/gencode.te.locus.salmon.out
        /public/groups/kimlab/indexes/gen.32.ucsc.rmsk.index
    ```
    
4. reference has big impact on mapping rate
   ```
   50.7606	
       kras.ipsc/data/exo.data/gen1c.day.5.exo.input.data/kras.3/gencode.salmon.out
       /public/groups/kimlab/indexes/gencode.32.v.1.index
       
   75.0142
       kras.ipsc/data/exo.data/gen1c.day.5.exo.input.data/kras.3/gencode.te.locus.salmon.out
       /public/groups/kimlab/indexes/gen.32.ucsc.rmsk.index
       
   36.8215
       kras.ipsc/data/exo.data/gen1c.day.5.exo.input.data/kras.3/quantFiles
	   /public/groups/kimlab/gen.v31.k.31.v.0.14.index
   ```
   
## split fastq files into mapped and unmapped reads
1. bin/salmonUnmapped.sh salmonIndexDir right.fastq left.fastq outdir
   - will create a quant.sf file from mapped reads
   ```
   d   = mapping type was determined as mapping to a decoy sequence
   u   = The entire pair was unmapped. No mappings were found for either the left or right read.
   m1  = Left orphan (mappings were found for the left (i.e. first) read, but not the right).
   m2  = Right orphan (mappinds were found for the right read, but not the left).
   m12 = Left and right orphans. Both the left and right read mapped, but never to the same transcript.
   ```
    
   - $outdir /aux_info/unmapped_names.txt 
   ```
   $head aux_info/unmapped_names.txt 
   NB551675:38:HTMGLAFXY:1:11101:22494:8399 u
   NB551675:38:HTMGLAFXY:1:11101:10728:8400 d
   NB551675:38:HTMGLAFXY:1:11101:21165:8401 d
   NB551675:38:HTMGLAFXY:1:11101:6081:8410 d
   NB551675:38:HTMGLAFXY:1:11101:19777:8412 u
   NB551675:38:HTMGLAFXY:1:11101:4640:8415 d
   NB551675:38:HTMGLAFXY:1:11101:11482:8418 d
   NB551675:38:HTMGLAFXY:1:11101:21024:8419 d
   NB551675:38:HTMGLAFXY:1:11101:25790:8423 u
   NB551675:38:HTMGLAFXY:1:11101:8557:8423 u
   ```
   - is d really m12? did we find a bunch of novel transcripts?
   ```
   ctrl.1.unmapped]$ cut -d " " -f 2 aux_info/unmapped_names.txt | sort | uniq -c
   519916 d
   39097 m1
   34534 m2
   747447 u
   ```

2. ???? selectUnmappedReadsFromFastq.sh
3. ??? selectMappedReadsFromFastq.sh
   
   

## bin/ table of contents


## Install 3rd party tools

salmon: The easiest way is to download the tarball from https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz. conda does not work you will get a very old version even if you use "-c bioconda"

download samtools source and build
```
wget https://sourceforge.net/projects/samtools/files/samtools/1.12/samtools-1.12.tar.bz2/download
tar -xvf download 
cd samtools-1.12/
./configure 
less README 
make
```

download seqtk. Makes it easy to extract unmapped reads by name/id from fastq
```
cd bin
git clone git@github.com:lh3/seqtk.git seqtk.git
cd seqtk.git
make
cd ..
ln -s seqtk.git seqtk
```

STAR: download linux and compile
```
$ wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz
$ tar -zxf 2.7.9a.tar.gz 
$ ln -s STAR-2.7.9a/bin/Linux_x86_64_static/STAR .
```

## creating STAR index

down load genome and annotations
```
$ cd /private/groups/kimlab/genomes.annotations
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
$ wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.p13.genome.fa.gz
$ gunzip gencode.v37.annotation.gtf.gz
$ gunzip GRCh38.p13.genome.fa.gz
```

build index
```
output=/private/groups/kimlab/indexes
data=/private/groups/kimlab/genomes.annotations
idxName=gne.38.p.13.v37.annotation.STAR.2.7.9a.index

# run the following using nohub
STAR --runMode genomeGenerate runThreadN 10 \
    --genomeDir ${output}/${indexName}/ \
    --genomeFastaFiles ${data}/GRCh38.p13.genome.fa \
    --sjdbGTFfile ${data}/gencode.v37.annotation.gtf 
```


## Starting env
some shell scripts may use tools like fastqc that are installed in the conda environment extraCellularRNA




# bin/ Table of Contents
- createTmpFile.sh
  - defines utility functions for creating a tmp file that will be deleted
  automatically if script is interupted
  
- findCandidateUnmappedDirs.sh
  - utlity script. useful for finding output directories you want to delete so you can
  - re-run part of the pipeline
  
- createStarIndex.sh
  - single use, creates reference need for STAR

## step 1) Salmon meta data mining prep pipeline
- samonUnmapped.sh
  - runs salmon quant with --writeUnmappedNames

- runSalmonAll.sh
  - call salmonUnmapped.sh on a list of with ref index gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx
  - examples of file list
  ```
  /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/output_forward_paired.fq.gz
  /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.2/output_forward_paired.fq.gz
  ```
  
  ```
  setsid sh -c ' set -x; \
      runSalmonAll.sh \
      /scratch/aedavids/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx \
      ../data/aale.laud.exo/aale.laud.exo.foward.paired.fq.gz.list.txt ' \
      2>&1 > runSalmonAll.sh.`../../extraCellularRNA/bin/dateStamp.sh`.out &
  ```
  
- makeUnmappedFastq.sh
  - parses the ${salmonOut}/aux_info/unmapped_names.txt file and create a separate
  pair of fastq files for each unamapped category, d, u, m1, m2, m12, plus a fastq
  file containing all the unmapped reads and one containing all the mapped reads
  
  - calls selectUnmappedReadsFromFastq.sh


  
## step 2) mineSalmonLogs.sh pipe line
- mineSalmonLogs.sh 
  - outputs tsv with columns
  - sampleName mappingRate salmonOut readCount index mates1 mates2
  
- masterCountReads.sh
  - for each of the Salmon unmapped reads type, d, u, m1, m2, m12, unmapped
  counts the number of reads
  
- countReads.sh
  - counts the number of reads in a fastq file

mineSalmonLogs.sh, masterCountReads.sh -> countReads.sh

## down stream analysis
- runSTARonUnmappedSalmon.sh
  - runs STAR with arguments 
    - --outReadsUnmapped Fastx # create fastq file with unmapped reads
    - --outSAMtype BAM SortedByCoordinate
  - will also create a bam index file
  
- runSTARUnmappedFastqc.sh
  - run fastqc on bam files created by STAR
  
- runSalmonUnmappedFastqc.sh
  - run fastqc on the unmapped fastq files created by makeUnmappedFastq.sh
  
- masterFastqc.sh
  - runs runSTARonUnmappedSalmon.sh and  runSTARUnmappedFastqc.sh

- runHTSeqCount.sh
  - finds all the Aligned.sortedByCoord.out.bam files produced by STAR and run htseq-count
  
- runMultiqc.sh
  - use to select directories multiqc should search
  

### unmappedReadsAnalysis/jupyterNotebooks

- mappingRatesPlots.ipynb
  - horizontal box plots of mapping rates
  - y axis is sample
    - data point u_b_5_k (unmapped bulk day 5 kras replicants)
    - ytick label is sample + salmon idx
  - x axis
    - percent

-  uniqueUnMappedCounts.ipynb
   - must be run on private servers
   - searches for htseq-count output files for unmapped reads for a given salmon idx
   - output uniqueUnMappedCounts.tsv
     - 3 columns, htseqCountFile, sampleId, numUnmappedUniqueGenes
     
   - uniqueMappedCounts.ipynb
     - counts the number of unique genes 
       -input file created by extraCellularRNA/R/notebooks/kras.ipsc.DESeq.normalize.v2.Rmd

-  geneCountPlots.ipynb

### extraCellularRNA/R/notebooks
    - R notebooks in a separate git repo
    - kras.ipsc.DESeq.normalize.v2.Rmd
      - uses DESeq to create normalized gene counts used in uniqueMappedCounts.ipynb
    
