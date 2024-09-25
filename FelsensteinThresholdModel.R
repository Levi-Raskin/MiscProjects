### Felsenstein threshold model


# libraries ---------------------------------------------------------------
library(phytools)
library(castor)
library(geiger)
library(ggplot2)
library(aplot)

# set up ------------------------------------------------------------------

#tree
tree <- phytools::pbtree(n = 10)
plot(tree)

# liabilities
l1 <- phytools::fastBM(tree, 
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
  for(j in 1:100){ #sampling 100 individuals in each population
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

ng <- 1000
burnIn <- 0.1*ng
pf <- 10
sf <- 10
acceptWiggles <- 0.5
searchSpread <- 0.05

propFunc <- function(vec){
  return(sum(vec == 0) / length(vec)) 
}
acceptFunc <- function(vec1, vec2){
  if(all(abs(vec1 - vec2) < acceptWiggles)){
    return(TRUE)
  }else{
    return(FALSE)
  }
  
}
phenotypeFunc <- function(testL){
    p1 <- c()
    for(j in 1:100){ 
      liability1 <- rnorm(1, mean = testL, sd = 1)
      if(liability1 >= t1){
        p1[j] <- 1  
      }else{
        p1[j] <- 0
      }
    }
    return(p1)
}
liabProposal <- function(vec, liabsToChange){
  vec2 <- vec
  vec2[liabsToChange] <- vec2[liabsToChange] + rnorm(length(liabsToChange), 0, searchSpread)
  return(vec2)
}

prop0vec1 <- lapply(phenotypeL1, propFunc)

t1Start <- rnorm(1, 0, 1)
liabStart <- phytools::fastBM(tree, 
                             sig2 = 1,
                             internal = T)
threshStart <- rnorm(1, 0, 1)

threshRes <- c()
liabRes <- list()

liabsToChange <- 1:20

for(i in 1:ng){
  liabsToChange <- 1:20
  j <- 1
  repeat{
    j <- j+1
    if(j > 10000 && i < 100){ 
      testLiab <- phytools::fastBM(tree, 
                                  sig2 = 1,
                                  internal = T)
      testPhenotype <- lapply(testLiab[1:Ntip(tree)], phenotypeFunc)
      print("random search")
    }else{
      testLiab <- liabProposal(liabStart,liabsToChange)
      testPhenotype <- lapply(testLiab[1:Ntip(tree)], phenotypeFunc)
      if(searchSpread < 2){
        searchSpread <- searchSpread*1.0001
      }
    }
   
    if(acceptFunc(unlist(lapply(testPhenotype, propFunc)), unlist(prop0vec1)) == TRUE){
      liabStart <- testLiab
      liabRes[[i]] <- testLiab
    
      if(i %% pf == 0){
        print(paste("NG: ", i))
      }
      if(acceptWiggles > 0.05){
        acceptWiggles <- acceptWiggles*0.99
      }
      searchSpread <- 0.05
      break
    }else{
      liabsToChange <- c()
      for(k in 1:length(testPhenotype)){
        if(abs(lapply(testPhenotype, propFunc)[[k]] - prop0vec1[[k]]) >= acceptWiggles){
          liabsToChange <- c(liabsToChange, k)
        }
      }
    }
  }
  
}

dat <- c()
for(i in burnIn:ng){
  dat[i] <- liabRes[[i]][3]  
}
hist(dat)
l1[3]

