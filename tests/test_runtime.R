library(IMaGES)
library(graph)
library(igraph)
library(sfsmisc)
library(lavaan)
library(Rgraphviz)
library(stats)
library(utils)

autism_driver <- function() {
  #get file locations
  
  sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source);
  
  #get filenames 
  filenames <- list.files("test/steve", pattern="autism*", full.names=TRUE)
  
  matrices = list()
  
  #import data
  for (i in 1:length(filenames)) {
    matrices[[i]] <- as.matrix(read.table(filenames[[i]], header=TRUE))
  }
  
  
  
  #run IMaGES on data
  results = new("IMaGES", matrices = matrices, penalty=3)
  
  
}

start.time <- Sys.time()

data(IMData)

#res2 <- IMaGES(matrices=data.list)
autism_driver()

end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste("Time: ", time.taken))


