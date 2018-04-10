library(IMaGES)
library(graph)
library(igraph)
library(sfsmisc)
library(lavaan)
library(Rgraphviz)


## Load predefined data


get_gmg <- function() {
  set.seed(50)
  p <- 8
  n <- 2000
  ## true DAG:
  vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
  gGtrue <- randomDAG(p, prob = 0.3, V = vars)
  #gmG  <- list(x = rmvDAG(n, gGtrue, back.compatible=TRUE), g = gGtrue)
  gmG8 <- list(x = rmvDAG(n, gGtrue),                       g = gGtrue)
  return(gmG8)
}

randomDAG <- function (n, prob, lB = 0.1, uB = 1, V = as.character(1:n))
{
  stopifnot(n >= 2, is.numeric(prob), length(prob) == 1,
            0 <= prob, prob <= 1,
            is.numeric(lB), is.numeric(uB), lB <= uB)
  edL <- vector("list", n)
  nmbEdges <- 0L
  for (i in seq_len(n - 2)) {
    listSize <- rbinom(1, n - i, prob)
    nmbEdges <- nmbEdges + listSize
    edgeList <- sample(seq(i + 1, n), size = listSize)
    weightList <- runif(length(edgeList), min = lB, max = uB)
    edL[[i]] <- list(edges = edgeList, weights = weightList)
  }
  ## i=n-1 separately
  ## (because of sample(7,1/) is actually sample(1:7,1) and not 7)
  listSize <- rbinom(1, 1, prob)
  if (listSize > 0) {
    nmbEdges <- nmbEdges + 1
    edgeList <- n
    weightList <- runif(1, min = lB, max = uB)
  } else {
    edgeList <- integer(0)
    weightList <- numeric(0)
  }
  edL[[n-1]] <- list(edges = edgeList, weights = weightList)
  if (nmbEdges > 0) {
    edL[[n]] <- list(edges = integer(0), weights = numeric(0))
    names(edL) <- V
    new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
  }
  else
    new("graphNEL", nodes = V, edgemode = "directed")
}

rmvDAG <-
  function(n, dag,
           errDist = c("normal", "cauchy", "t4", "mix", "mixt3", "mixN100"),
           mix = 0.1, errMat = NULL, back.compatible = FALSE,
           use.node.names = !back.compatible)
  {
    ## Purpose: Generate data according to a given DAG (with weights) and
    ## given node distribution (rows: number of samples; cols: node values in
    ## topological order)
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## - n     : Number of samples
    ## - dag   : Graph object containing the DAG and weights
    ## - errDist: "normal" or "mix" for pure standard normal node distribution
    ##           or mixing with standard cauchy
    ## - mix   : Percentage of points sampled from standard cauchy
    ## ----------------------------------------------------------------------
    ## Author: Markus Kalisch, Date: 26 Jan 2006;  Martin Maechler
    
    ## check input &  initialize variables
    stopifnot(is(dag, "graph"),
              (p <- length(nodes(dag))) >= 2)
    
    ##  as(.,"matrix") now {for some versions of 'graph' pkg} is 0/1
    ## weightMatrix <- t(as(dag,"matrix"))
    weightMatrix <- if(back.compatible) wgtMatrix.0(dag) else wgtMatrix(dag)
    
    ## check if top. sorted
    nonZeros <- which(weightMatrix != 0, arr.ind = TRUE)
    if (nrow(nonZeros) > 0) {
      if (any(nonZeros[,1] - nonZeros[,2] < 0) || any(diag(weightMatrix) != 0))
        stop("Input DAG must be topologically ordered!")
    }
    
    errDist <- match.arg(errDist)
    if(grepl("^mix", errDist))
      eMat <- function(outs) { # (n,p)
        X <- c(rnorm(n*p - length(outs)), outs)
        matrix(sample(X), nrow = n)
      }
    if(is.null(errMat)) {
      ## generate errors e_i
      errMat <-
        switch(errDist,
               "normal" = matrix(rnorm  (n*p),  nrow = n),
               "cauchy" = matrix(rcauchy(n*p),  nrow = n),
               "t4" =     matrix(rt(n*p, df = 4), nrow = n),
               "mix"    = eMat(rcauchy(round(mix*n*p))),
               "mixt3"  = eMat(     rt(round(mix*n*p), df = 3)),
               "mixN100"= eMat(  rnorm(round(mix*n*p), sd = 10)))
    }
    else { ## check & use 'errMat' argument:
      stopifnot(!is.null(dim.eM <- dim(errMat)),
                dim.eM == c(n,p), is.numeric(errMat))
    }
    if(use.node.names)
      colnames(errMat) <- nodes(dag) # == colnames(weightMatrix)
    
    ## compute X matrix X_i
    if (sum(weightMatrix) > 0) {
      X <- errMat
      for (j in 2:p) { ## uses X[*, 1:(j-1)] -- "recursively" !
        ij <- 1:(j-1)
        X[,j] <- X[,j] + X[, ij, drop = FALSE] %*% weightMatrix[j, ij]
      }
      X
    }
    else
      errMat
  }

wgtMatrix.0 <- function(g, transpose = TRUE)
{
  ## Purpose: work around "graph" package's  as(g, "matrix") bug
  ## ----------------------------------------------------------------------
  ## ACHTUNG: mat_[i,j]==1 iff j->i,
  ## whereas with as(g,"matrix") mat_[i,j]==1 iff i->j
  ## ----------------------------------------------------------------------
  ## Arguments: g: an object inheriting from (S4) class "graph"
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, based on Seth Falcon's code;  Date: 12 May 2006
  
  ## MM: another buglet for the case of  "no edges":
  if(numEdges(g) == 0) {
    p <- length(nd <- nodes(g))
    return( matrix(0, p,p, dimnames = list(nd, nd)) )
  }
  ## Usual case, when there are edges:
  if(!("weight" %in% names(edgeDataDefaults(g))))
    edgeDataDefaults(g, "weight") <- 1L
  w <- unlist(edgeData(g, attr = "weight"))
  ## we need the *transposed* matrix typically:
  tm <- if(transpose) t(as(g, "matrix")) else as(g, "matrix")
  ## now is a 0/1 - matrix (instead of 0/wgts) with the 'graph' bug
  if(any(w != 1)) ## fix it
    tm[tm != 0] <- w
  ## tm_[i,j]==1 iff i->j
  tm
}



wgtMatrix <- function(g, transpose = TRUE) {
  res <- as(g, "matrix") # from 'graph' package, now reliable (we hope)
  if (transpose) ## default!
    t(res) else res
}

#' Auxiliary function reading an edge list (as used in the constructors
#' of DAGs) out of an adjacency matrix or a graphNEL object
#' @param from adjacency matrix, graphNEL object, or object inherited
#'  from ParDAG
#' @return list of in-edges; length of list = number of vertices,
#' entries for i-th vertex = indices sources of in-edges
inEdgeList <- function(from)
{
  if (is.matrix(from)) {
    p <- nrow(from)
    stopifnot(p == ncol(from))
    lapply(1:p, function(i) which(from[, i] != 0))
  } else if (class(from) == "graphNEL") {
    nodeNames <- graph::nodes(from)
    edgeList <- lapply(graph::inEdges(from), function(v) match(v, nodeNames))
    names(edgeList) <- NULL
    edgeList
  } else if (length(grep(".*ParDAG", class(from)) == 1)) {
    from$.in.edges
  }else {
    stop(sprintf("Input of class '%s' is not supported.", class(from)))
  }
}


make_data <- function(prob) {
  
  set.seed(40)
  p <- 8
  n <- 5000
  ## true DAG:
  vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
  gGtrue <- randomDAG(p, prob = prob, V = vars)
  dataset <- list(x = rmvDAG(n, gGtrue), g = gGtrue)
  #dataset<-dataset$x
}


#create the DAGs generated by gmG8 plus added noise
create_im_dags <- function(num_sets) {
  gmG8 <- get_gmg()
  #initial seed for generation of dataset
  start_seed <- 50
  data_list <- list()
  p <- 8

  n <- 2000
    ## true DAG:
  vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
  gGtrue <- gmG8$g
  
  for (i in 1:num_sets) {
   #set.seed(start_seed)
    #gGtrue <- randomDAG(p, prob = 0.3, V = vars)
    #inject noise into DAGs using rnorm
    set8 <- list(x = gmG8$x + matrix(rnorm(p*n,mean=0,sd=1),n,p), g = gGtrue)
    #set8 <- list(x = rmvDAG(n, gGtrue)+ matrix(rnorm(p*n,mean=0,sd=runif(1,0,0.5)),n,p), g = gGtrue)
    #set8 <- list(x = rmvDAG(n, gGtrue), g = gGtrue)
    data_list[[i]] <- set8
    #increment start seed
    #start_seed = start_seed + 20
  }
  return(data_list)
}

create_all <- function(num_sets) {
  gmG8 <- get_gmg()
  #initial seed for generation of dataset
  start_seed <- 50
  set.seed(start_seed)
  data_list <- list()
  
  p <- 8
  n <- 2000
  ## true DAG:
  vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
  gGtrue <- gmG8$g
  
  for (i in 1:num_sets) {
    set.list = list()
    for (k in 1:i) {
      #gGtrue <- randomDAG(p, prob = 0.3, V = vars)
      #inject noise into DAGs using rnorm
      set8 <- list(x = gmG8$x + matrix(rnorm(p*n,mean=0,sd=1),n,p), g = gGtrue)
      #set8 <- list(x = rmvDAG(n, gGtrue)+ matrix(rnorm(p*n,mean=0,sd=runif(1,0,0.5)),n,p), g = gGtrue)
      #set8 <- list(x = rmvDAG(n, gGtrue), g = gGtrue)
      set.list[[k]] <- set8
      #increment start seed
      #start_seed = start_seed + 20
    }
    data_list[[i]] <- set.list
  }
  return(data_list)
}

create_all_alt <- function(num_sets) {
  gmG8 <- get_gmg()
  #initial seed for generation of dataset
  start_seed <- 50
  set.seed(start_seed)
  data_list <- list()
  
  p <- 8
  n <- 2000
  ## true DAG:
  vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
  gGtrue <- gmG8$g
  
  set.list <- list()
  for (i in 1:num_sets) {
    #set.list = list()
    #for (k in 1:i) {
      #gGtrue <- randomDAG(p, prob = 0.3, V = vars)
      #inject noise into DAGs using rnorm
    set8 <- list(x = gmG8$x + matrix(rnorm(p*n,mean=0,sd=2),n,p), g = gGtrue)
      #set8 <- list(x = rmvDAG(n, gGtrue)+ matrix(rnorm(p*n,mean=0,sd=runif(1,0,0.5)),n,p), g = gGtrue)
      #set8 <- list(x = rmvDAG(n, gGtrue), g = gGtrue)
    set.list[[i]] <- set8
    data_list[[i]] <- set.list
      #increment start seed
      #start_seed = start_seed + 20
  }
  return(data_list)
}

#create and return list of score objects to be passed into IMaGES
create_scores <- function(datasets) {
  scores <- list()
  for (i in 1:length(datasets)) {
    scores[[i]] <- new("GaussL0penObsScore", datasets[[i]]$x)
  }
  return(scores)
}

#function used for evaluating original, GES-style runs on individual datasets
run_originals <- function(datasets) {
  fits = list()
  for (i in 1:length(datasets)) {
    images <- new("IMaGES", scores = list(datasets[[i]]), penalty=3)
    fits[[i]] <- images
  }
  return(fits)
}


## Plot the estimated essential graph and the true DAG
plot_graph <- function(fit) {
  if (require(Rgraphviz)) {
    par(mfrow=c(1,2))
    plot(fit, main = "Estimated CPDAG")
    plot(gmG8$g, main = "True DAG")
  } else { ## alternative:
    str(fit, max=2)
    as(as(fit$essgraph,"graphNEL"),"Matrix")
  }
}


#calculates precision, recall, and f-measure of graphs by finding:
#    true positives (edges in both graphs)
#    false positives (edges in the generated graph but not in the true DAG)
#    false negatives (edges that should have been in the generated graph but weren't)
#and using these for the calculations
find_error <- function(graph) {
  
  gmG8 <- get_gmg()
  
  true <- inEdgeList(gmG8$g)
  
  graph <- inEdgeList(graph)
  
  positives <- 0
  true_positives <- 0
  false_positives <- 0
  false_negatives <- 0
  print(graph)
  for (i in 1:length(graph)) {
      for (j in 1:length(graph[[i]])) {

        if (((length(graph[[i]]) > 0) && (length(true[[i]]) > 0)) && (graph[[i]][[j]] %in% true[[i]])) {
          positives <- positives + 1
        }
        else if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && !(graph[[i]][[j]] %in% true[[i]])) {
          false_positives <- false_positives + 1
        }
      }

    if (length(true[[i]]) > 0) {
      for (k in 1:length(true[[i]])) {
        if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && (true[[i]][[k]] %in% graph[[i]])) {
          true_positives <- true_positives + 1
        }
        else if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && !(true[[i]][[k]] %in% graph[[i]])) {
          false_negatives <- false_negatives + 1
        }
      }
    }
  }
  
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  
  f_measure <- 2 * ((precision * recall) / (precision + recall))
  return(f_measure)
}

find_both <- function(graph) {
  
  gmG8 <- get_gmg()
  
  true <- inEdgeList(gmG8$g)
  
  graph <- inEdgeList(graph)
  
  positives <- 0
  true_positives <- 0
  false_positives <- 0
  false_negatives <- 0
  print(graph)
  for (i in 1:length(graph)) {
    for (j in 1:length(graph[[i]])) {
      
      if (((length(graph[[i]]) > 0) && (length(true[[i]]) > 0)) && (graph[[i]][[j]] %in% true[[i]])) {
        positives <- positives + 1
      }
      else if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && !(graph[[i]][[j]] %in% true[[i]])) {
        false_positives <- false_positives + 1
      }
    }
    
    if (length(true[[i]]) > 0) {
      for (k in 1:length(true[[i]])) {
        if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && (true[[i]][[k]] %in% graph[[i]])) {
          true_positives <- true_positives + 1
        }
        else if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && !(true[[i]][[k]] %in% graph[[i]])) {
          false_negatives <- false_negatives + 1
        }
      }
    }
  }
  
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  
  print("made it here")
  
  res <- list(.precision=precision, .recall=recall)
  
  #f_measure <- 2 * ((precision * recall) / (precision + recall))
  return(res)
}


#######
#create accuracy measure
#precision, recall, f-measure
#0 or 1 for each edge/direction

#function used for plotting the error across multiple iterations of IMaGES and showing
#how it decreases due to the "IMaGES effect"
plot_error <- function(results) {
  
  inv_measures <- list()
  
  for (i in 1:length(results)) {
    #calculate f-measures for each graph in the set
    f_list <- list()
    # for (k in 1:length(results[[i]])) {
    #   #use 1 minus error for the sake of presentation
    #   print(results[[i]][[k]])
    #   f_list[[k]] <- 1 - find_error(results[[i]][[k]]$results$.global$.in.edges)
    # 
    # }
    #print(f_list)
    #use mean error (although should all be the same) for graph
    #inv_measures[[i]] <- mean(unlist(f_list))
    print(results[[i]]$results$.global$.graph)
    inv_measures[[i]] <- 1 - find_error(results[[i]]$results$.global$.graph)
    print(inv_measures[[i]])
  }

  #cannot be in list form for plotting
  plot_measures <- unlist(inv_measures)
  print(plot_measures)
  
  png(filename="f-measure.png", width=800, height=400)
  
  
  #dev.new(width=10, height=5)
  #plot(plot_measures, type="o", col="blue", main="Error", ylim=c(0,0.5))
  noise = 1.8
  plot(plot_measures, main=paste("IMaGES Error for noise value", noise), type="o", col='blue', ylim=c(0,1), xlab='', ylab='')
  at <- seq(from=0, to=length(results), by=length(results)/20)
  title(xlab="Number of datasets")
  title(ylab="Error")
  axis(side = 1, at = at)
  
  dev.off()
}

plot_both <- function(results) {
  
  prec_measures <- list()
  rec_measures <- list()
  
  for (i in 1:length(results)) {
    #calculate f-measures for each graph in the set
    #f_list <- list()
    # for (k in 1:length(results[[i]])) {
    #   #use 1 minus error for the sake of presentation
    #   print(results[[i]][[k]])
    #   f_list[[k]] <- 1 - find_error(results[[i]][[k]]$results$.global$.in.edges)
    # 
    # }
    #print(f_list)
    #use mean error (although should all be the same) for graph
    #inv_measures[[i]] <- mean(unlist(f_list))
    print(results[[i]]$results$.global$.graph)
    both <- find_both(results[[i]]$results$.global$.graph)
    prec_measures[[i]] <- 1 - both$.precision
    rec_measures[[i]] <- 1 - both$.recall
    #print(inv_measures[[i]])
  }
  
  #cannot be in list form for plotting
  precision_measures <- unlist(prec_measures)
  print(precision_measures)
  
  png(filename="precision.png", width=800, height=400)
  
  
  #dev.new(width=10, height=5)
  #plot(plot_measures, type="o", col="blue", main="Error", ylim=c(0,0.5))
  #noise = 0.01
  plot(precision_measures, main=paste("IMaGES Precision Error"), type="o", col='blue', ylim=c(0,1), xlab='', ylab='')
  at <- seq(from=0, to=length(results), by=length(results)/20)
  title(xlab="Number of datasets")
  title(ylab="Error")
  axis(side = 1, at = at)
  
  dev.off()
  
  recall_measures <- unlist(rec_measures)
  print(recall_measures)
  
  png(filename="recall.png", width=800, height=400)
  
  
  plot(recall_measures, main=paste("IMaGES Recall Error"), type="o", col='blue', ylim=c(0,1), xlab='', ylab='')
  at <- seq(from=0, to=length(results), by=length(results)/20)
  title(xlab="Number of datasets")
  title(ylab="Error")
  axis(side = 1, at = at)
  
  dev.off()
}


#driver for individual GES-like runs
driver <- function() {
  #change to how many graphs you want
  num_sets <- 3
  
  gmG8 <- get_gmg()
  
  #generate DAGS
  # dags <- create_im_dags(num_sets)
  # #create score objects
  # scores <- create_scores(dags)
  # #find GES-like fits using IMaGES
  # orig_fits <- run_originals(scores)
  # 
  # #print(orig_fits[[1]][[1]][[2]])
  # 
  # #plot in.edges for each graph
  # for (i in 1:length(orig_fits)) {
  #   plotIMGraph(orig_fits[[i]]$results$.global)
  # }
  
  #now do the same thing for IMaGES
  
  #create DAGS
  im_run_dags <- create_im_dags(num_sets)
  
  #create score objects
  im_run_scores <- create_scores(im_run_dags)
  #run IMaGES
  im_fits <- new("IMaGES", scores = im_run_scores, penalty=2)
  
  #plot results
  par(mfrow=c(1,2))
  plotIMGraph(im_fits$results$.global)
  plotAll(im_fits)
  
  plotMarkovs(im_fits)
}

driver_prob <- function() {
  #change to how many graphs you want
  num_sets <- 3
  
  #now do the same thing for IMaGES
  
  #create DAGS
  #im_run_dags <- create_im_dags(num_sets)
  
  dataset1 <- make_data(0.3)
  dataset2 <- make_data(0.3)
  dataset3 <- make_data(0.3)
  
  #create score objects
  im_run_scores <- create_scores(list(dataset1,dataset2,dataset3))
  #im_run_scores <- create_scores(list(dataset1))
  #run IMaGES
  im_fits <- new("IMaGES", scores = im_run_scores, penalty=3)
  
  
  plotAll(im_fits)
  plotIMGraph(im_fits$results$.global)
  
  #plot results
  # for (i in 1:length(im_fits$results)) {
  #   plot_graph(im_fits$results[[i]][[2]])
  # }
}

#driver for calculation of errors across runs of increasing size
plot_driver <- function() {
  #change to number of sets to iterate up to
  num_sets <- 25
  
  #generate gmG8 data
  #gmG8 <- get_gmg()

  #stores fits for each set size
  result_sets <- list()
  
  for (k in 1:num_sets) {
    #create DAGs
    im_run_dags <- create_im_dags(k)
    #create score objects
    im_run_scores <- create_scores(im_run_dags)
    #run IMaGES
    im_fits <- new("IMaGES", scores = im_run_scores, penalty=3)
    #append results to result_sets
    result_sets[[k]] <- im_fits
    
    print(k)
    #plotIMGraph(im_fits$results$.global)
    
    
  }
  #calculates errors for each of the result sets
  plot_error(result_sets)
  plot_both(result_sets)
  
  # #plots individual sets (might creash computer as it's a lot of plots)
  # for (k in 1:length(result_sets)) {
  #   for (i in 1:length(result_sets[[k]])) {
  #     plot_graph(result_sets[[k]]$results[[i]][[2]])
  #   }
  # }
  
  
}

test_driver <- function() {
  #change to number of sets to iterate up to
  num_sets <- 25
  
  #generate gmG8 data
  #gmG8 <- get_gmg()
  
  #stores fits for each set size
  result_sets <- list()
  
  #for (i in 1:num_sets) {
    #create DAGs
  im_run_dags <- create_all_alt(num_sets)
  #print(im_run_dags)
    #create score objects
  for (k in 1:length(im_run_dags)) {
    im_run_scores <- create_scores(im_run_dags[[k]])
    #run IMaGES
    im_fits <- new("IMaGES", scores = im_run_scores, penalty=20)
    #append results to result_sets
    result_sets[[k]] <- im_fits
    print(k)
    #plotIMGraph(result_sets[[k]]$results$.global)
    
  }
    
    
  #calculates errors for each of the result sets
  plot_error(result_sets)
  plot_both(result_sets)
  
  # #plots individual sets (might creash computer as it's a lot of plots)
  # for (k in 1:length(result_sets)) {
  #   for (i in 1:length(result_sets[[k]])) {
  #     plot_graph(result_sets[[k]]$results[[i]][[2]])
  #   }
  # }
  
  
}

#driver for running IMaGES on autism data. works but the plot still isn't showing up properly
#it might be due to the fact that the labels aren't included?
autism_driver <- function() {
  #get file locations
  
  library(graph)
  library(igraph)
  library(sfsmisc)
  library(lavaan)
  library(Rgraphviz)
  
  sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source);
  
  #get filenames 
  #filenames <- list.files("test/14-19", pattern="SB*", full.names=TRUE)
  filenames <- list.files("test/", pattern="VA*", full.names=TRUE)
  
  matrices = list()
  
  #import data
  for (i in 1:length(filenames)) {
  #for (i in 1:2) {
    #this might be causing the plotting issue but i'm not yet sure of a workaround
    matrices[[i]] <- as.matrix(read.table(filenames[[i]], header=TRUE))
    #test1 <- as.matrix(read.table(filenames[[i]]))#, skip=1))
    #test2 <- as.matrix(read.table(filenames[[i]]), skip=1)
    
  }
  

  
  #run IMaGES on data
  results = new("IMaGES", matrices = matrices, penalty=3)
  
  
  plotIMGraph(results$results$.global)
  #plotIMGraph(results$results$.alt)
  plotMarkovs(results)
  plotAll(results)
  
  par(mfrow=c(2,5))
  for (i in 1:length(matrices)) {
    #corrplot(cor(matrices[[i]]), type = "upper", 
    #         tl.col = "black", tl.srt = 45)#, is.corr=FALSE)
    cor(matrices[[i]])
  }
  
  #plot resulting DAGs
  # for (i in 1:length(results)) {
  #   par(mfrow=c(1,2))
  #   plot(results$results[[i]][[2]], main = "Estimated CPDAG")
  # }

}

powerball_driver <- function() {
  #get file locations
  sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source);
  
  #get filenames 
  filenames <- list.files("test/powerball", pattern="pb*", full.names=TRUE)
  matrices = list()
  
  #import data
  for (i in 1:length(filenames)) {
    #for (i in 1:2) {
    #this might be causing the plotting issue but i'm not yet sure of a workaround
    matrices[[i]] <- as.matrix(read.table(filenames[[i]], header=FALSE))
    #test1 <- as.matrix(read.table(filenames[[i]]))#, skip=1))
    #test2 <- as.matrix(read.table(filenames[[i]]), skip=1)
    
  }
  
  #run IMaGES on data
  results = new("IMaGES", matrices = matrices, penalty=3, num.markovs=6)
  
  plotIMGraph(results$results$.global)
  
  plotMarkovs(results)
  
  #plot resulting DAGs
  # for (i in 1:length(results)) {
  #   par(mfrow=c(1,2))
  #   plot(results$results[[i]][[2]], main = "Estimated CPDAG")
  # }
  
}

test_dataset <- function() {
  #get file locations
  sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source);
  
  #get filenames c
  filenames <- list.files("test/d9", pattern="dataset*", full.names=TRUE)
  matrices = list()
  
  #import data
  for (i in 1:length(filenames)) {
    #for (i in 1:2) {
    #this might be causing the plotting issue but i'm not yet sure of a workaround
    matrices[[i]] <- as.matrix(read.table(filenames[[i]], header=FALSE))
    #test1 <- as.matrix(read.table(filenames[[i]]))#, skip=1))
    #test2 <- as.matrix(read.table(filenames[[i]]), skip=1)
    
  }
  
  #run IMaGES on data
  results = new("IMaGES", matrices = matrices, penalty=3, num.markovs=5)
  
  
  
  plotIMGraph(results$results$.global)
  plotAll(results)
  plotMarkovs(results)
  
  #plot resulting DAGs
  # for (i in 1:length(results)) {
  #   par(mfrow=c(1,2))
  #   plot(results$results[[i]][[2]], main = "Estimated CPDAG")
  # }
  
}

convert <- function(from) {
  edgeList <- lapply(from$.in.edges, function(v) from$.nodes[v])
  names(edgeList) <- from$.nodes
  result <- new("graphNEL",
                nodes = from$.nodes,
                edgeL = edgeList,
                edgemode = "directed")
  return(reverseEdgeDirections(result))

}

