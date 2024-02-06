#!/usr/bin/Rscript
require(data.table)
require(foreach)
require(iterators)
require(doParallel)
require(stats)

## Get input probabilities file as command line argument
args = commandArgs(trailingOnly=TRUE)
sample_name = args[1]

## Get beta distribution mixture model parameters from infile
beta_model <- fread('/data/HG002/betaDistributionMixture.model',header=TRUE)
meth_shape1 <- as.numeric(beta_model[mixture=='Methylated' & shape==1,n])
meth_shape2 <- as.numeric(beta_model[mixture=='Methylated' & shape==2,n])
unmeth_shape1 <- as.numeric(beta_model[mixture=='Unmethylated' & shape==1,n])
unmeth_shape2 <- as.numeric(beta_model[mixture=='Unmethylated' & shape==2,n])
cat(paste0('Beta mixture parameters:\n\tMethylated distribution: shape1 = ',meth_shape1,
           '; shape2 = ',meth_shape2,
           '\n\tUnmethylated distribution: shape1 = ',unmeth_shape1,
           '; shape2 = ',unmeth_shape2,'\n'))

## Add beta distribution log likelihood ratios and predictions for LLR < -2 (Unmethylated), -2 < LLR < 2 (Uncertain), LLR > 2 (Methylated)
container_data_path = '/ifs/data/research/projects/juliet/projects/methylation_lrs/'
probabilityFiles <- list.files('/data',
                               pattern=paste0(sample_name,'.chr[0-9]*.ReadMethylationProbabilities.txt$'),
                               full.names=TRUE,
                               recursive=TRUE)
probabilityFiles <- gsub(container_data_path,'',probabilityFiles)
cat(paste0('File list:\n\t',probabilityFiles))

numCores = detectCores() - 1
registerDoParallel(numCores)

for( file in probabilityFiles ) {
  
  infile <- fread(file,header=TRUE)
  cpg_ids <- unique(infile[,cpg_id])
  
  llrTable <- foreach( cpg = cpg_ids, .combine = 'rbind', .inorder = TRUE ) %dopar% {
    
    chr <- infile[cpg_id==cpg,chr]
    start <- infile[cpg_id==cpg,start]
    end <- infile[cpg_id==cpg,end]
    strands <- unique(infile[cpg_id==cpg,strand])
    
    strandTable <- foreach( s = strands, .combine = 'rbind' ) %do% {
      
      max_val = 0.9970938
      min_val = 0.00290625
      delta = 0.001
      lrs_probabilities <- as.numeric(unlist(strsplit(infile[cpg_id==cpg & strand==s,methylation_probabilities],split=',')))
      lrs_probabilities <- ifelse(lrs_probabilities==0,min_val - delta,
                                  ifelse(lrs_probabilities==1,max_val + delta,lrs_probabilities))
      
      meth_pred <- pbeta(lrs_probabilities, meth_shape1, meth_shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
      unmeth_pred <- pbeta(lrs_probabilities, unmeth_shape1, unmeth_shape2, ncp = 0, lower.tail = FALSE, log.p = FALSE)
      llr <- paste(log(meth_pred/unmeth_pred),collapse=',')
      pred <- paste(ifelse(llr %in% c(Inf,-Inf),'Uncertain',
                           ifelse(llr > 2,'Methylated',
                                  ifelse(llr < -2,'Unmethylated','Uncertain'))),collapse=',')
      
      return(data.table('cpg_id' = cpg,
                        'chr' = chr,
                        'start' = start,
                        'end' = end,
                        'strand' = s,
                        'methylation_probabilities' = paste(lrs_probabilities,collapse=','),
                        'llrs' = llr,
                        'predictions' = pred))
      
    }
    
    return(strandTable)

  }
  
  ## Write outfile
  outfilePath <- unique(gsub('.chr[0-9]*.','.',probabilityFiles))
  
  if( !file.exists(outfilePath) ) {
    
    write.table(llrTable,outfilePath,
                col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
    
  } else {
    
    write.table(llrTable,outfilePath,
                col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t',append=TRUE)
    
  }
  
}
