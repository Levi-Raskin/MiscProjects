library(phytools)

#read in your data in Nexus format here
data <- read.nexus.data("C:/INPUT/FILEPATH/HERE.nex")

#set the output file's path
file <- file("C:/OUTPUT/FILEPATH/HERE.txt")

numTaxa <- length(data)
numChar <- length(data[[1]])

string <- paste(numTaxa, " ", numChar, collapse = "")

for(i in 1:length(data)){
  string <- c(string, paste(names(data)[i], " ", paste(data[[i]], collapse = "")))
}

writeLines(
  string,
  file,
  sep = "\n")

close(file)