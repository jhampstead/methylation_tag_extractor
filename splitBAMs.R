#!/usr/bin/Rscript
require(data.table)
require(foreach)
require(iterators)
require(doParallel)

## Get input BAM as command line argument
args = commandArgs(trailingOnly=TRUE)
bam_path = args[1]

## Ensure the BAM is indexed; if the BAM is unindexed, index it
if( !file.exists(paste0('/data/',bam_path,'.bai')) ) {
  
  system(paste0('/usr/local/bin/samtools index /data/',bam_path))
  
}

## Subset input BAM into chromosomes
## Index the split BAMs
chromosomes = c(paste0('chr',1:22))

bam_name = tail(unlist(strsplit(bam_path,split='/')),1)
bam_subdirs = paste(head(unlist(strsplit(bam_path,split='/')),-1),collapse='/')
split_dir = paste0('/data/',bam_subdirs,'/split_bams/')
dir.create(split_dir,showWarnings = FALSE)

numCores <- detectCores() - 1
registerDoParallel(numCores)
foreach( chr = chromosomes ) %dopar% {
  
  system(paste0('/usr/local/bin/samtools view -b /data/',bam_path,' ',chr,' > ',split_dir,gsub('.bam','',bam_name),'.',chr,'.bam'))
  system(paste0('/usr/local/bin/samtools index ',split_dir,gsub('.bam','',bam_name),'.',chr,'.bam'))
  
}