## ----setup, include = FALSE, eval=TRUE-----------------------------------
knitr::opts_chunk$set(
	eval = TRUE,
	collapse = TRUE,
	comment = "#>",
	include = FALSE
)

## ----echo=TRUE, warning=FALSE, include=TRUE------------------------------
require(IMaGES)

data(IMData)

#run IMaGES
im.fits <- IMaGES(matrices=data.list, penalty=3, num.markovs=5)


## ----eval=FALSE, warning=FALSE, include=TRUE-----------------------------
#  require(IMaGES)
#  
#  ## Load predefined data
#  data(IMData)
#  
#  #run IMaGES
#  im.fits <- IMaGES(matrices=data.list, penalty=3, num.markovs=5)
#  
#  #plot global graph and all individual graphs with own SEM data
#  plotAll(im.fits)

## ----eval=FALSE, warning=FALSE, include=TRUE-----------------------------
#  require(IMaGES)
#  ## Load predefined data
#  data(IMData)
#  
#  #run IMaGES
#  im.fits <- IMaGES(matrices=data.list, penalty=3, num.markovs=5)
#  
#  #plot global graph alongside Markov Equivalence Class
#  plotMarkovs(im.fits)

## ----eval=FALSE, warning=FALSE, include=TRUE-----------------------------
#  require(IMaGES)
#  
#  ## Load predefined data
#  data(IMData)
#  
#  #run IMaGES
#  im.fits <- IMaGES(matrices=data.list, penalty=3, num.markovs=5)
#  
#  #plot individual graph
#  plotIMGraph(im.fits$results$.single.graph[[1]])
#  

