#############################
# EPiK: Predicts changes in kinetic parameters
#############################

####################   All variables you need to change for running the pipeline on your data is here ################
data=do.call(cbind, mRNAcounts) #Your mRNA counts
gene="Gfi1b"
pop1_ids=c(1:123) #indices of one of the subpopulations you want to compare
pop2_ids=c(124:247) #indices of the second subpopulation you want to compare
######################################################################################################################

#######  Advanced parameters ######
likelihoods=hypergeo 
max_mRNA=200
load("min_thresholds_marg.RData") #min_thresholds_marg is the threshold for assigning parameter changes via the MP method
load("min_thresholds_sub.RData") #min_thresholds_sub is the threshold for assigning parameter changes via the subset method
###################################

#output#
#bicOut= BIC predictions
#mpOut= MP predictions: values corresponding to Kon, Koff, Kt. 
#mpOutBest= the set of predicted changes by the mp method, as determined by the thresholds set in min_thresholds_marg
#subOut= subset predictions: values corresponding to Kon, Koff, Kt
#subOutBest= the set of predicted changes by subset method, as determined by the thresholds set in min_thresholds_sub

###Likelihood calculation at core of method
pops=list(pop1_ids, pop2_ids)
data2=sapply(c(0:max_mRNA), function(i){
  sapply(pops, function(j){
    length(which(data[gene,j]==i))
  })
})

outmat= data2 %*% as.matrix(likelihoods)

###Setup indices into nice tables for quick calculations with matrix multiplications
p=paste(allTuples[,1], allTuples[,2])
kon_koff=unique(p)
kon_koff=sapply(kon_koff, function(i){which(p==i)})

#get all (koff, kt) pairs and all indices in allTuples
p=paste(allTuples[,2], allTuples[,3])
koff_kt=unique(p)
koff_kt=sapply(koff_kt, function(i){which(p==i)})

#get all (kon, kt) pairs and all indices in allTuples
p=paste(allTuples[,1], allTuples[,3])
kon_kt=unique(p)
kon_kt=sapply(kon_kt, function(i){which(p==i)})

#get all (kon, kt) pairs and all indices in allTuples
p=allTuples[,1]
kon_same=unique(p)
kon_same=sapply(kon_same, function(i){which(p==i)})

p=allTuples[,2]
koff_same=unique(p)
koff_same=sapply(koff_same, function(i){which(p==i)})

#get all (koff, kt) pairs and all indices in allTuples
p=allTuples[,3]
kt_same=unique(p)
kt_same=sapply(kt_same, function(i){which(p==i)})



############################BIC portion##################################################

######First, lets get likelihood, assuming every parameter is the same

ml_none=max(outmat[1,]+outmat[2,])

######Now, lets get likelihood assuming one parameter has changed

out_temp1=apply(koff_kt, c(1,2), function(i){outmat[1,i]})
out_temp2=apply(koff_kt, c(1,2), function(i){outmat[2,i]})
ml_kon=max(apply(out_temp1, 2, function(i){max(i)})+apply(out_temp2, 2, function(i){max(i)}))

out_temp1=apply(kon_kt, c(1,2), function(i){outmat[1,i]})
out_temp2=apply(kon_kt, c(1,2), function(i){outmat[2,i]})
ml_koff=max(apply(out_temp1, 2, function(i){max(i)})+apply(out_temp2, 2, function(i){max(i)}))

out_temp1=apply(kon_koff, c(1,2), function(i){outmat[1,i]})
out_temp2=apply(kon_koff, c(1,2), function(i){outmat[2,i]})
ml_kt=max(apply(out_temp1, 2, function(i){max(i)})+apply(out_temp2, 2, function(i){max(i)}))

#####Now, let's calculate the likelihoods for two TFs changing

out_temp1=apply(kt_same, c(1,2), function(i){outmat[1,i]})
out_temp2=apply(kt_same, c(1,2), function(i){outmat[2,i]})
ml_kon_koff=max(apply(out_temp1, 2, function(i){max(i)})+apply(out_temp2, 2, function(i){max(i)}))

out_temp1=apply(koff_same, c(1,2), function(i){outmat[1,i]})
out_temp2=apply(koff_same, c(1,2), function(i){outmat[2,i]})
ml_kon_kt=max(apply(out_temp1, 2, function(i){max(i)})+apply(out_temp2, 2, function(i){max(i)}))

out_temp1=apply(kon_same, c(1,2), function(i){outmat[1,i]})
out_temp2=apply(kon_same, c(1,2), function(i){outmat[2,i]})
ml_koff_kt=max(apply(out_temp1, 2, function(i){max(i)})+apply(out_temp2, 2, function(i){max(i)}))

#####Now let's calculate the likelihoods for all three TFs changing

ml_all=max(outmat[1,])+max(outmat[2,])

##########Turn Likelihoods into BICs

mls=c(ml_none, ml_kon, ml_koff, ml_kt, ml_kon_koff, ml_kon_kt, ml_koff_kt, ml_all)
names(mls)=c("none", "kon", "koff", "kt", "kon, koff", "kon, kt", "koff, kt", "all")
bicFactor=log(2*length(data[,1])) -log(2*pi)
bics=bicFactor*c(0, 1, 1, 1, 2, 2, 2, 3)+(-2*mls)
bicOut=names(mls)[which.min(bics)]



################################### Marginal Probability-based strategy #############################
pop1=outmat[1,]
pop2=outmat[2,]

#probability of Kt being the same

out_temp1=apply(kt_same, c(1,2), function(i){pop1[i]})
out_temp2=apply(kt_same, c(1,2), function(i){pop2[i]})

val=max(out_temp1)-10
sum_temp1=log(sum(exp(out_temp1-val)))+val
out_temp1=out_temp1-sum_temp1
val=max(out_temp1)-10
temp1=sapply(c(1:length(out_temp1[1,])), function(k){
  x=(sum(sapply(c(1:length(out_temp1[,1])), function(c){
    exp(out_temp1[c,k]-val)  
  })))
  log(x)+val
})

val=max(out_temp2)-10
sum_temp2=log(sum(exp(out_temp2-val)))+val
out_temp2=out_temp2-sum_temp2
val=max(out_temp2)-10
temp2=sapply(c(1:length(out_temp2[1,])), function(k){
log(sum(sapply(c(1:length(out_temp2[,1])), function(c){
    exp(out_temp2[c,k]-val)  
  })))+val
})

#now calculate probability that both populations have that kt
#aka: log-scale multiplication (which becomes addition)
temp=temp1+temp2

#now calculate the sum of these possibilities, in non-log scale
#in non-log scale these -inf become 0, so no worries
p_kt=log(sum(sapply(temp, function(i){ exp(i-max(temp)+10)})))+max(temp)-10


#Do the same for Koff
out_temp1=apply(koff_same, c(1,2), function(i){pop1[i]})
out_temp2=apply(koff_same, c(1,2), function(i){pop2[i]})
val=max(out_temp1)-10
sum_temp1=log(sum(exp(out_temp1-val)))+val
out_temp1=out_temp1-sum_temp1
val=max(out_temp1)-10
temp1=sapply(c(1:length(out_temp1[1,])), function(k){
  x=(sum(sapply(c(1:length(out_temp1[,1])), function(c){
    exp(out_temp1[c,k]-val)  
  })))
  log(x)+val
})

val=max(out_temp2)-10
sum_temp2=log(sum(exp(out_temp2-val)))+val
out_temp2=out_temp2-sum_temp2
val=max(out_temp2)-10
temp2=sapply(c(1:length(out_temp2[1,])), function(k){
  log(sum(sapply(c(1:length(out_temp2[,1])), function(c){
    exp(out_temp2[c,k]-val)  
  })))+val
})

#now calculate probability that both populations have that koff
#aka: log-scale multiplication (which becomes addition)
temp=temp1+temp2
#now calculate the sum of these possibilities, in non-log scale
#in non-log scale these -inf become 0, so no worries
p_koff=log(sum(sapply(temp, function(i){ exp(i-max(temp)+10)})))+max(temp)-10

#Same put for Kon
out_temp1=apply(kon_same, c(1,2), function(i){pop1[i]})
out_temp2=apply(kon_same, c(1,2), function(i){pop2[i]})
val=max(out_temp1)-10
sum_temp1=log(sum(exp(out_temp1-val)))+val
out_temp1=out_temp1-sum_temp1
val=max(out_temp1)-10
temp1=sapply(c(1:length(out_temp1[1,])), function(k){
  x=(sum(sapply(c(1:length(out_temp1[,1])), function(c){
    exp(out_temp1[c,k]-val)  
  })))
  log(x)+val
})

val=max(out_temp2)-10
sum_temp2=log(sum(exp(out_temp2-val)))+val
out_temp2=out_temp2-sum_temp2
val=max(out_temp2)-10
temp2=sapply(c(1:length(out_temp2[1,])), function(k){
  log(sum(sapply(c(1:length(out_temp2[,1])), function(c){
    exp(out_temp2[c,k]-val)  
  })))+val
})


#now calculate probability that both populations have that kon
#aka: log-scale multiplication (which becomes addition)
temp=temp1+temp2

#now calculate the sum of these possibilities, in non-log scale
#in non-log scale these -inf become 0, so no worries
p_kon=log(sum(sapply(temp, function(i){ exp(i-max(temp)+10)})))+max(temp)-10

#output
mpOut=c(p_kon, p_koff, p_kt)
names(mpOut)=c("kon", "koff", "kt")
mpOutBest=names(mpOut)[which(mpOut<min_thresholds_marg)]




########################## Subset method ########################

calc_subsetPrediction<-function(mat1, mat2, subsamples){
  output_subset_realData=sapply(c(1:length(mat1[1,])), function(gene){
    subs=sapply(c(1:100), function(i){
      data=cbind(mat1[sample(length(mat1[,1]), subsamples),gene],
                 mat2[sample(length(mat2[,1]), subsamples),gene])
      data=apply(data, c(1,2), function(i){if(i>max_mRNA){max_mRNA}else{i}})
      
      data2=sapply(c(0:max_mRNA), function(i){
        
        sapply(c(1:2), function(j){
          
          length(which(data[,j]==i))
          
        })
        
      })
      
      
      
      outmat= data2 %*% as.matrix(likelihoods)
      
      t=apply(outmat, 1, function(k){which.max(k)})
      
      allTuples[t,]
      
    })
    
    sapply(c(1:3), function(param){
      
      a=sapply(c(1:100), function(i){subs[param,i][[1]][1]})
      
      b=sapply(c(1:100), function(i){subs[param,i][[1]][2]})
      
      ks.test(a, b)$statistic  
      
    })
    
  })
  
  return(output_subset_realData)
  
}

subsamples=31
mat1=data[,pop1_ids]
mat2=data[,pop2_ids]
allSubSets=calc_subsetPrediction(t(mat1), t(mat2), subsamples)
subOut=allSubSets[,which(rownames(data)==gene)]
names(subOut)=c("kon", "koff", "kt")
subOutBest=names(subOut)[which(subOut>min_thresholds_sub)]
