# IMaGES

This is the repository for the R implementation of the IMaGES algorithm. This project was overseen by Professor Stephen Jose Hanson of the Rutgers University Brain Imaging Center. The repository started as a fork of [pcalg](https://cran.r-project.org/web/packages/pcalg/index.html) and is now a standalone product. The additional code and changes were written/made by Noah Frazier-Logue.

This algoritm elaborates on the GES algorithm by using a global score across the supplied datasets and operating over the datasets concurrently to determine the representative graph(s) with the best goodness of fit.

**NOTE**: This software is in beta! If you come across any issues while using this package or have any suggestions for improvement, submit a pull request.

### Installation

To install from this repository, simply run these commands in an R shell:

```
> library(devtools)
> install_github("noahfl/IMaGES")
```

TODO: Add stuff about CRAN when that becomes relevant.

### Usage

```R
#matrices should be a list of >= 1 datasets with an optional header
matrices <- list(matrix1, matrix2,...)

im.results <- new("IMaGES", matrices=matrices, penalty=3, num.markovs=5)

#plot individual graph, in this case the global graph
plotIMGraph(im.results$results$.global)

#plot Markov Equivalence Class (size specified by num.markovs)
plotMarkovs(im.results)

#plot global graph with SEM data, and all individual datasets' SEM data
#imposed on the global graph
plotAll(im.results)
```
