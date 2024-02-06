#!/usr/bin/Rscript
require(data.table)
require(foreach)
require(iterators)
require(doParallel)

## Get input BAM as command line argument
args = commandArgs(trailingOnly=TRUE)
bam_path = args[1]

## Read strings that match pattern C[+m5]C[+m7]C[+m10] (have an tagged 5mC call)
modifiedReads <- system(paste0('/usr/local/bin/samtools mpileup -M /data/',bam_path,' | grep "\\[+m[0-9]\\]"'),
                         intern=TRUE)
modifiedReads <- strsplit(modifiedReads,split='\t')

numCores <- detectCores() - 1
registerDoParallel(numCores)
modifiedReadsDT <- foreach( i = 1:length(modifiedReads), .combine = 'rbind' ) %dopar% {
  
  return(data.table(t(modifiedReads[[i]])))
  
}

## Exclude reads with indels that shift bases relative to the reference sequence
excludeIndels <- function( s ) {
  
  if( grepl('\\+[0-9]*|\\-[0-9]*',s) ) {
    
    val <- grep('\\+[0-9]{1,20}$|\\-[0-9]{1,20}$',unlist(strsplit(s,split='[ATCGNatcgn\\*\\]]',perl=TRUE)),value=TRUE)
    bases <- unlist(strsplit(s,split='(?<=[ATCGNatcgn\\*\\]](?!\\[))',perl=TRUE))
    
    for( v in val ) {
      
      r = abs(as.numeric(v)) - 1
      pos = rev(grep(v,bases,fixed=TRUE))
      
      for( p in pos ) {
        
        vec = seq(p,p+r,1)
        bases <- bases[-vec]
        
      }
      
    }
    
  } else {
    
    bases <- s
    
  }
  
  s_out <- paste0(bases,collapse='')
  return(s_out)
  
}

modifiedReadsDT <- modifiedReadsDT[,V5:=excludeIndels(V5),by=V5]

## Get per-read 5mC probabilities
getMethylationProbabilities <- function( s ) {
  
  tmp <- as.numeric(gsub('\\+m','',as.character(do.call(cbind, strsplit(s, "\\[+|\\]", perl = TRUE)))))
  out <- paste((tmp[!is.na(tmp)])/256,collapse=',')
  
  
}

modifiedReadsDT <- modifiedReadsDT[,methylation_probabilities:=getMethylationProbabilities(V5),by=V5]

## Get the modified base for each read
getModifiedBase <- function( s ) {
  
  in_split <- unlist(strsplit(s,split=''))
  matches <- grep('\\[',unlist(strsplit(s,split='')),perl=TRUE) - 1
  out <- paste(unique(toupper(in_split[matches])),collapse=';')
  return(out)

}

modifiedReadsDT <- modifiedReadsDT[,modified_base:=getModifiedBase(V5),by=V5]

## Remove reads without a C or G modified base; these are either: 
## - Reads where modifications occurred on bases removed during indel filtering, or
## - Reads where both a C and G base are reported modified (unreliable/misaligned reads)
modifiedReadsDT <- modifiedReadsDT[modified_base %in% c('C','G')]

## Create CpG dinucleotide IDs
modifiedReadsDT <- modifiedReadsDT[,strand:=ifelse(modified_base=='C','Forward','Reverse'),by=modified_base]
modifiedReadsDT[,continguous_positions:=c(0, abs(diff(as.numeric(modifiedReadsDT[,V2]))) == 1)]
modifiedReadsDT[,cpg_id:=ifelse(continguous_positions==1 & modified_base=='G',paste0('cg_',gsub('chr','',V1),'_',shift(V2,1L,type="lag"),'_',V2),
                                ifelse(continguous_positions==0 & modified_base=='G',paste0('cg_',gsub('chr','',V1),'_',as.numeric(V2)-1,'_',V2),paste0('cg_',gsub('chr','',V1),'_',V2,'_',as.numeric(V2)+1)))]
modifiedReadsDT[,start:=ifelse(continguous_positions==1 & modified_base=='G',shift(V2,1L,type="lag"),
                               ifelse(continguous_positions==0 & modified_base=='G',as.numeric(V2)-1,V2))]
modifiedReadsDT[,end:=ifelse(continguous_positions==1 & modified_base=='G',V2,
                             ifelse(continguous_positions==0 & modified_base=='G',V2,as.numeric(V2)+1))]

modifiedReadsDTOut <- unique(modifiedReadsDT[,c('cpg_id','V1','start','end','strand','methylation_probabilities')])
colnames(modifiedReadsDTOut) <- c('cpg_id','chr','start','end','strand','methylation_probabilities')

## Write outfile
bam_name = tail(unlist(strsplit(bam_path,split='/')),1)
bam_subdirs = paste(head(grep('split_bams',unlist(strsplit(bam_path,split='/')),invert=TRUE,value=TRUE),-1),collapse='/')
methprob_dir = paste0('/data/',bam_subdirs,'/read_methylation_probabilities/')
dir.create(methprob_dir,showWarnings = FALSE)
outfile_path <- paste0(methprob_dir,gsub('.bam','',bam_name),'.ReadMethylationProbabilities.txt')

write.table(modifiedReadsDTOut,outfile_path,
            col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
