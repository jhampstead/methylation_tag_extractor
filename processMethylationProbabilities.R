#!/usr/bin/Rscript
require(data.table)
require(foreach)
require(iterators)
require(doParallel)

## Get input probabilities file as command line argument
args = commandArgs(trailingOnly=TRUE)
prob_path = args[1]
include_bisulfite = args[2]

methylseq_filtered <- fread('/data/HG002/HG002_methylseq_filtered.txt')
methylation_read_probabilities <- fread(paste0('/data/',prob_path))

combinedTable <- merge(methylseq_filtered,methylation_read_probabilities,by=c('chr','end'),all.x=FALSE,all.y=FALSE)

if( include_bisulfite==TRUE ) {
  
  getBisulfiteCalls <- function( c, s, e ) {
    
    return(paste(methylseq_filtered[chr==c & start==s & end==e,prop_methylated_reads]/100,collapse=','))
    
  }
  
  combinedTable <- combinedTable[,bisulfite_prop:=getBisulfiteCalls(chr,start.x,end),by=c('chr','start.x','end')]
  
  cpg_ids <- unique(combinedTable[,cpg_id])
  
  numCores = detectCores()
  registerDoParallel(numCores)
  
  probabilitiesTable <- foreach( cpg = cpg_ids, .combine = 'rbind' ) %dopar% {
    
    calls <- ifelse(length(unique(combinedTable[cpg_id==cpg,strand]))==2,'Both strands',
                    ifelse(unique(combinedTable[cpg_id==cpg,strand])=='Forward','Forward only','Reverse only'))
    
    return(data.table('cpg_id' = cpg,
                      'calls' = calls,
                      'lrs_probability' = as.numeric(unlist(strsplit(combinedTable[cpg_id==cpg,methylation_probabilities],split=','))),
                      'methylation_status_bisulfite' = unique(combinedTable[cpg_id==cpg,methylation_status]),
                      'bisulfite_prop' = as.numeric(unlist(strsplit(unique(combinedTable[cpg_id==cpg,bisulfite_prop]),split=',')))))
    
  }
  
} else {
  
  cpg_ids <- unique(combinedTable[,cpg_id])
  
  numCores = detectCores()
  registerDoParallel(numCores)
  
  probabilitiesTable <- foreach( cpg = cpg_ids, .combine = 'rbind' ) %dopar% {
    
    calls <- ifelse(length(unique(combinedTable[cpg_id==cpg,strand]))==2,'Both strands',
                    ifelse(unique(combinedTable[cpg_id==cpg,strand])=='Forward','Forward only','Reverse only'))
    
    return(data.table('cpg_id' = cpg,
                      'calls' = calls,
                      'lrs_probability' = as.numeric(unlist(strsplit(combinedTable[cpg_id==cpg,methylation_probabilities],split=',')))))

  }
  
}

## Write outfile
filepath = unlist(strsplit(prob_path,split='/'))
filename = paste0(gsub('.txt','',filepath[length(filepath)]),'_Beta.txt')
dir.create(paste0('/data/',filepath[1],'/mixture_model_probabilities/'),showWarnings = FALSE)

write.table(probabilitiesTable,paste0('/data/',filepath[1],'/mixture_model_probabilities/',filename),
            col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')