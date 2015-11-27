######################################################
#Pipeline step #1: Normalisation and scaling
#This file has R scripts to transform single cell qPCR data to scaled mRNA counts for use in the remaining portions of the pipeline
######################################################

####################   All variables you need to change for running the pipeline on your data is here ################

#These are the names of the files you would like to normalise; place the corresponding files in the "data" folder:
#The maximum Ct by default is 28.  If this is not the case, please change max_val_Ct accordingly
files=c("CLP_WT_Raw", "GMP_WT_Raw", "HSC_WT_Raw", "LMPP_WT_Raw", "PreM_WT_Raw")

#These are the names of the items in the final list; they correspond to the file names above
popNames=c("CLP", "GMP", "HSC", "LMPP", "PreM")

#These are the housekeeping genes to normalise to:
housekeepers=c("Polr2a", "Ubc")

#genes of interest: the subset of genes that will be normalised and utilised for downstream analysis
genesOfInterest=c("Erg", "Eto2", "Fli1", "Gata1", "Gata2", "Gfi1", "Gfi1b", "hHex", "Ldb1", "Lmo2", "Lyl1", "Meis1", "Mitf", "Nfe2", "PU.1", "Runx1", "SCL", "Tel")
######################################################################################################################


#######  Advanced parameters ######
max_val_Ct=28 #This is the cap on the Ct value.  
max_trusted_Ct=15 #This is the maximum trusted Ct score after normalisation to housekeeping genes
max_mRNA=200 #This is the largest value in the hypergeo lookup table
###################################

#function to normalise data (please load this function before proceeding)
#expression is a table of raw qPCR Ct values, indexed by [cell, gene]
#house is the list of housekeeping genes
#filterGene is a gene for which you would like to remove any cells that have zero qPCR signal from downstream analysis
#max_val is the listed qPCR value when no signal was detected by qPCR
#reset_max is useful if you would like to replace max_val in the normalised dataset
#The output is a normalised dataset
normalize <- function(expression, house, max_val=28, reset_max=28){
  expression[which(expression>max_val, arr.ind=2)]=max_val
  #expression=expression[which(expression[,filterGene]!=max_val),]
  expression2 = (expression - max_val) - (rowMeans(expression[house]) - max_val)
  expression2[which(expression==max_val, arr.ind=2)]=reset_max
  return(expression2);
}


#read values from file and normalise, output is in: results_tf_normal
results_tf_normal=list();
for(file in files){
  print(file)
  geneExp=read.table(paste("data/", file, sep=""), header=TRUE);
  geneExp=normalize(geneExp, housekeepers, max_val=max_val_Ct) # genes normalised to
  results_tf_normal[[file]] <- geneExp;
}
results_tf_normal_all=do.call(rbind, results_tf_normal)
#find TF-specific scaling factor that will ensure that all mRNAcounts will be under 200 in subsequent procedures
removePolIIgonners=results_tf_normal_all[which(apply(results_tf_normal_all, 1, function(i){ length(which(i[housekeepers]==28))==0})),]
tf_scalefactors=apply(removePolIIgonners[, genesOfInterest]
                      , 2, function(i){max_mRNA/(2^(max_trusted_Ct-min(i)))}) 

#This the most crucial bit: 

mRNAcounts=lapply(c(1:length(results_tf_normal)), function(pop){
  removePolIIgonners=results_tf_normal[[pop]][which(apply(results_tf_normal[[pop]], 1, function(i){ length(which(i[housekeepers]==28))==0})),]
  results_tf_only=removePolIIgonners[,genesOfInterest ] 
  apply(results_tf_only, 1, function(i){
    temp=(2^(max_trusted_Ct-i)*tf_scalefactors) 
    sapply( temp, function(j){round(j)})
  })  
})

names(mRNAcounts)=popNames


