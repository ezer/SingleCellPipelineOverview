
##################################
# Overview of Single Cell Analysis Pipeline
##################################

##############################################################
#Step 0: load pre-calculated tables of P( mRNA | Kon, Koff, Kt)
##############################################################
#load("hypergeo.RData")

hyperSubs=lapply(c(1:10), function(i){
  load(paste("hyper_", i, ".RData", sep=""))
  hyper_temp
})
hypergeo=do.call(cbind, hyperSubs)
load("allTuples.RData")
#This was generated using Mathematica.  See hypergeoParamSaving workbook for example.
#hypergeo: a matrix that represents the lookup table of Probability(mRNA count | Kon, Koff, Kt).  The rows represent the mRNA count and the columns represent the parameter set
#allTuples: a matrix that contains the parameter values corresponsing to the columns in hypergeo.  The columns are Kon, Koff, and Kt.


##############################################################
#Step 1: normalise raw single cell qPCR data and turn Ct into a "transcript count" metric
##############################################################

#See normaliseSingleCellData.R
#Beforehand, make tab-delimited tables for the raw Ct values, 
#where each column represents a different gene and each row represents a different gene

#The output is a list called mRNAcounts, where each element in the list is a table of transcript counts for each input file

#############################################################
#Step 2: SABEC: Simulated Annealing for Bursty Expression Clustering
#############################################################

#For an example of running the clustering algorithm (once) see: sabec.R
#However, it is important to run the code many times to get accurate results.  It is recommended to do this on the cluster
#See: FILES for examples of how to run this on a Condor Cluster.  

#The output files are as follows:
#A) outname_cluster_id files: the current classification for iteration "id" of the algorithm; the highest number is the one that corresponds to the "best" clustering
#B) tupleIndex__outname__id: this lists the index of allTuples corresponding to the most likely set of kinetic parameters for each gene in cluster "id"
#C) likelihoods__outname__id: the log-likelihoods of each gene having these kinetic parameters in cluster "id"

#Sample output files are found in the folder "sampleOutput"

#############################################################
#Step 3: EPiK: Estimating Parameters in Kinetics
#############################################################

#For an example of running EPiK, see: epik.R








#data2 transforms data into a new format: 
#Rows represent genes, columns represent mRNA count, and the values represent the number of cells that have a particular mRNA count of a particular gene
  data2=sapply(c(0:200), function(i){ #The reason this is 0 to 200 is that this is the range for which the parameters have been pre-calculated in hypergeo
    apply(pop, 1, function(j){
      length(which(j==i))
    })
  })
  

#This the most crucial bit: 
#outmat is a matrix that represents the likelihoods of each parameter set, given the data, for each gene
#Rows represent genes, columns represent a particular parameter set
#by parameter set, I mean the combination of (Kon, Koff, and Kt)
  outmat= data2 %*% as.matrix(hypergeo)
  
#Here I identify the parameter set that has the maximum likelihood for each gene
  ml=apply(outmat, 1, function(i){
    which.max(i)
  })
 
#here I lookup the specific value of the most likely parameter set, which is contained in allTuples
#tup: each row represents the gene, the columns represent the most likely values of (Kon, Koff, and Kt), respectively
  tup=allTuples[ml,]
  rownames(tup)=rownames(outmat)





#In order to see how to do this on a Camgrid cluster, see:
#camgrid_0_to_20_EstimatedParamsLikeData_6and7clusters.R
#job_0_to_20_EstimatedParamsDataLike_6and7clusters.txt


########################
# merging hundreds of iterations of the clustering algorithm
########################
i=100
file=(Sys.glob(paste("Estim*LikeData*randomSet_", i, "_iter_", 1, "_d*", sep="")))
arr=read.table(file[1], header=TRUE)
names=rownames(arr)
print(names)
f=c("6_or_7", "6_or_7", "8_or_9", "8_or_9", "10_or_11", "10_or_11")
clust=c(6,7,8,9,10)
#sLikeData_6_or_7clusters_0_to_20_Table_sim_testPops_7_randomSet_9_iter_6_d_evenbiggerH
output_realData=lapply(c(1:5), function(id){
  run=f[id]
  clust_num=clust[id]
  sapply(c(1:50), function(j){
    file=(Sys.glob(paste("Estim*LikeData*", run, "*Pops_", clust_num, "_randomSet_", i, "_iter_", j, "_d*", sep="")))
    print(file[1])
    read.table(file[1], header=TRUE)
  })
})

########################
# EPiK: estimating kinetic parameter changes
########################



