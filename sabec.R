#############################
# SABEC: Clustering algorithm
#############################

####################   All variables you need to change for running the pipeline on your data is here ################
data=do.call(cbind, mRNAcounts) #Your mRNA counts
numPops=5 # number of clusters you expect (K)
outname="clusterTest"
######################################################################################################################

#######  Advanced parameters ######
likelihoods=hypergeo # the table of P(mRNA | Kon, Koff, Kt), precalculated in Mathematica
tempChange=10 #Temperature parameter that influences the speed of convergence
numIter=100 #max number of iters (should converge much quicker than that)
###################################

#Note: this should take about 40 minutes to run.

numGenes=length(rownames(data))
paramOpts=length(likelihoods[1,])
#Sort cells into hypothesized populations randomly
class=sample(c(1:numPops), length(data[1,]), replace=TRUE); 
tableNum=1
#now do the iterative portion of the algorithm
for(iterations in c(1:numIter)){
  #outMats sets up the datastructure that will contain outMat for each population 
  outMats=lapply(c(1:numPops), function(pop){
    matrix(rep(0, paramOpts*numGenes), nrow=paramOpts, ncol=numGenes)
  })
  
  #This fills up outMats
  num=0
  for(count in c(1:length(likelihoods[,1]))){    
    countVectors=lapply(c(1:numPops), function(pop){
      apply(data[,which(class==pop)], 1, function(cells){
        length(which(cells==num));
      })})
    likelihoodVector=likelihoods[count,]
    outMats= lapply(c(1:numPops), function(pop){
      outMats[[pop]]+(as.numeric(likelihoodVector) %*% t(as.numeric(countVectors[[pop]])))
    })
    num=num+1;
  }
  
  #record previous cluster assignments
  write.table(class, paste(outname, "_cluster_", tableNum, sep=""))
  
  #params selects the maximum likelihood parameter estimates for each cluster
  params=lapply(c(1:numPops), function(pop){
    apply(as.matrix(outMats[[pop]]), 2, function(i){ which.max(i) })
  })
  
  #cellLikelihoods assigns the likelihood of a cell to come from each population
  cellLikelihoods=sapply(c(1:numPops), function(pop){
    apply(data, 2, function(i){
      sum(as.numeric(sapply(c(1:numGenes), function(j){
        likelihoods[(1+i[j]), params[[pop]][j]]
      })))
    })
  })
  preclass=class
  if(length(unique(preclass))!=length(unique(class))){
    break;
  }
  #assign new classes probabilistically
  class=apply(cellLikelihoods, 1, function(p){
    p_temp=(-sum(p)+p)
    p_temp=p_temp/sum(p_temp)
    p_temp=p_temp^(tempChange*tableNum)
    sample(which(p_temp>0), 1, prob=p_temp[which(p_temp>0)])
  })
  
  #check to see whether the algorithm has converged (fewer than 5% classificiation changes)
  percent_equal=length(which(preclass==class))/length(class)
  if(percent_equal>0.95){
    break;
  }
  tableNum=tableNum+1;
}


#####Save output data
for(pop in c(1:numPops)){
  #print("start saving table")
  filename=paste("tupleIndex_", outname,"_", pop, sep="_")
  #save table
  write.table(apply(outMats[[pop]], 2, function(col){ which.max(col)}), filename)
  
  filename=paste("likelihoods_", outname,"_", pop, sep="_")
  #save table
  write.table(apply(outMats[[pop]], 2, function(col){ max(col)}), filename)
  
  
}


