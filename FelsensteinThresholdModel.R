### Felsenstein threshold model


# libraries ---------------------------------------------------------------
library(phytools)
library(castor)
library(geiger)
library(ggplot2)
library(aplot)

# set up ------------------------------------------------------------------

#tree
tree <- pbtree(n = 5) #pure birth tree with 10 tips
tree$tip.label <- 1:Ntip(tree)
plot(tree)

# liabilities
l1 <- fastBM(tree, 
                       sig2 = 1,
                       internal = T) #simulate l1 liabilities on the tree
ancL1 <- l1[(Ntip(tree) + 1 ):(Ntip(tree) + Nnode(tree))] #subset liabilities at the node
taxL1 <- l1[1:Ntip(tree)] #subset liabilities at the tips

l2 <- l1 * 1.2 +  runif(1, min = -0.5, max = 0.5) #generate set of libilities that is correlated with l1
ancL2 <- l2[(Ntip(tree) + 1 ):(Ntip(tree) + Nnode(tree))] #subset liabilities at the node
taxL2 <- l2[1:Ntip(tree)] #subset liabilities at the tips

#threshold
t1 <- runif(1, min = -0.5, max = 0.5) #threshold for trait 1
t2 <- runif(1, min = -0.5, max = 0.5) #threshold for trait 2

#observed phenotypes:
phenotypeL1 <- list()
phenotypeL2 <- list()
for(i in 1:Ntip(tree)){
  p1 <- c()
  p2 <- c()
  for(j in 1:1000){ #sampling 1000 individuals in each population
    liability1 <- rnorm(1, mean = taxL1[i], sd = 1)
    if(liability1 >= t1){
      p1[j] <- 1  
    }else{
      p1[j] <- 0
    }
    liability2 <- rnorm(1, mean = taxL2[i], sd = 1)
    if(liability2 >= t1){
      p2[j] <- 1  
    }else{
      p2[j] <- 0
    }
  }
  phenotypeL1[[i]] <- p1
  phenotypeL2[[i]] <- p2
}

#example genotype/phenotype
ibm <- c("#648FFF","#785EF0","#DC267F", "#FE6100", "#FFB000")
p1 <- ggplot()+
  geom_density(data = as.data.frame(rnorm(1000, mean = l1[1], sd = 1)), aes(x = rnorm(1000, mean = l1[1], sd = 1)), fill = ibm[1])+
  geom_vline(xintercept = t1)
p2 <- ggplot()+
  geom_histogram(data = as.data.frame(phenotypeL1[[1]]), aes(x=phenotypeL1[[1]], y = after_stat(count / sum(count))), fill = ibm[2])+
  scale_x_continuous(breaks = c(0,1))
p3 <- p1/p2
p3  


# inferring ancestral liabilities -----------------------------------------

ng <- 1000 #number of MCMC generations
burnIn <- 0.1*ng
pf <- 10 #print freq
acceptableDiff <- 0.1

propFunc <- function(vec){
  return(sum(vec == 0) / length(vec)) 
} #proportion of individuals in pop with phenotype 0
acceptFunc <- function(vec1, vec2){
  if(all(abs(vec1 - vec2) < acceptableDiff)){
    return(TRUE)
  }else{
    return(FALSE)
  }
  
} #acceptance/rejection of proposed new liability
phenotypeFunc <- function(testL){
  p1 <- c()
  for(j in 1:1000){ 
    liability1 <- rnorm(1, mean = testL, sd = 1)
    if(liability1 >= t1){
      p1[j] <- 1  
    }else{
      p1[j] <- 0
    }
  }
  return(p1)
} #generates phenotypes based on liabilities


prop0vec1 <- lapply(phenotypeL1, propFunc) #propotion 0 in phenotype L1 (simulated above; our observed data)

startLiabilities <- fastBM(tree,
                           internal = T) #starting liabilities

resList <- list()
for(i in 1:ng){
  j <- 1 #counter var
  repeat{
    testPhenotypes <- lapply(lapply(startLiabilities, phenotypeFunc), propFunc) #test phenotypes based on the test liabilities
    if(acceptFunc(unlist(testPhenotypes[1:Ntip(tree)]), unlist(prop0vec1))){
      #choose node to update:
      p <- sample(1 : (Nnode(tree) + Ntip(tree)), 1) #sample one of the tree's nodes
      pAnc <- tree$edge[which(tree$edge[,2] == p), 1] #get the singular ancestor of that node
      pDesc <- tree$edge[which(tree$edge[,1] == p), 2] #get the descendants of that node
      if(length(pDesc) > 1 &&  length(pAnc) == 1){ #if p is an internal node that is not the root
        pAncVal <- startLiabilities[which(names(startLiabilities) == pAnc)] #ancestor value
        pDescVal1 <- startLiabilities[which(names(startLiabilities) == pDesc[1])] #descendant 1 value
        pDescVal2 <- startLiabilities[which(names(startLiabilities) == pDesc[2])] #descendant 2 value
        mean <- ((1/castor::get_pairwise_distances(tree, p, pAnc)) * pAncVal + 
                   (1/castor::get_pairwise_distances(tree, p, pDesc[1])) * pDescVal1 + 
                   (1/castor::get_pairwise_distances(tree, p, pDesc[1])) * pDescVal2
        ) / ((1/castor::get_pairwise_distances(tree, p, pAnc)) + 
               (1/castor::get_pairwise_distances(tree, p, pDesc[1])) + 
               (1/castor::get_pairwise_distances(tree, p, pDesc[1])))
        variance <- 1/ ((1/castor::get_pairwise_distances(tree, p, pAnc)) + 
                          (1/castor::get_pairwise_distances(tree, p, pDesc[1])) + 
                          (1/castor::get_pairwise_distances(tree, p, pDesc[1])))
        startLiabilities[which(names(startLiabilities) == p)] <- rnorm(1, mean, sqrt(variance))
      }else if(length(pAnc) == 1 && length(pDesc) == 0){
        #tip; mean = pAnc val; var = sqrt branch length
        pAncVal <- startLiabilities[which(names(startLiabilities) == pAnc)]
        startLiabilities[which(names(startLiabilities) == p)] <- rnorm(1, pAncVal, sqrt(castor::get_pairwise_distances(tree, p, pAnc)))
      }else if(length(pAnc) == 0 && length(pDesc) > 0){
        pDescVal1 <- startLiabilities[which(names(startLiabilities) == pDesc[1])]
        pDescVal2 <- startLiabilities[which(names(startLiabilities) == pDesc[2])]
        mean <- ( (1/castor::get_pairwise_distances(tree, p, pDesc[1])) * pDescVal1 + 
                   (1/castor::get_pairwise_distances(tree, p, pDesc[1])) * pDescVal2
        ) / (  (1/castor::get_pairwise_distances(tree, p, pDesc[1])) + 
               (1/castor::get_pairwise_distances(tree, p, pDesc[1])))
        variance <- 1/ ( (1/castor::get_pairwise_distances(tree, p, pDesc[1])) + 
                          (1/castor::get_pairwise_distances(tree, p, pDesc[1])))
        startLiabilities[which(names(startLiabilities) == p)] <- rnorm(1, mean, sqrt(variance))
      }
      resList[[i]] <- startLiabilities
      break
    }else{
      
      if(i == 1 || j > 50000){ #mercy kill j>50000
        startLiabilities <- fastBM(tree,
                                   internal = T)
      }else{
        p <- sample(1 : Ntip(tree), 1)
        pAnc <- tree$edge[which(tree$edge[,2] == p), 1]
        pAncVal <- startLiabilities[which(names(startLiabilities) == pAnc)]
        testLiab <- startLiabilities
        testLiab[which(names(startLiabilities) == p)] <- rnorm(1, pAncVal, sqrt(castor::get_pairwise_distances(tree, p, pAnc)))
        
        testPhenotypesInt <- lapply(lapply(testLiab, phenotypeFunc), propFunc) #test phenotypes based on the test liabilities
    
        if(acceptFunc(unlist(testPhenotypesInt[1:Ntip(tree)]), unlist(prop0vec1))){
          startLiabilities <- testLiab
        }
      }
      
      
      if(j %% 100 == 0){
        print(paste("NG: ", i, "on rep ", j, " with no suitable liability found"))
      }
    }
    j <- j+1
  }
  
  if(i %% pf == 0){
    print(paste("NG: ", i))
  }
}

dat <- c()
for(i in burnIn:ng){
  dat<- c(dat, resList[[i]][8])
}
hist(dat) 
l1[8]

datDF <- data.frame("mean_liability" = dat)
ggplot()+
  geom_density(data = datDF, aes(x=mean_liability), fill = ibm[1])+
  geom_vline(xintercept = l1[8])

