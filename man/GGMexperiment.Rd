% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LITree-class.R
\name{GGMexperiment}
\alias{GGMexperiment}
\title{R6 Class for running experiment of Gaussian Graphical Model}
\description{
The class aims at comparing multiple Gaussian Graphical Model Inference procedure when a ground
truth is available
}
\section{Usage}{

\preformatted{experiment = GGMexperiment$new(X.list, adjmat, methods=c("glasso","em.latent.trees"))
}
}

\section{Arguments}{

\code{- X.list} A list of data.frame, which will be used as input data for the inference

\code{- adjmat} The ground truth adjacency matrix used for the evaluation (roc curve)

\code{- nb.missing.var} Number of missing variable (0 by default)

\code{- methods} A vector of characters list the methods to be tested (see GGMfit for the list of possible methods)

\code{- fit.number} Number of evaluation point (20 by default)

\code{- nb.sample}  Number of data frame in X.list (set automatically while initializing the object)

\code{- K.score.array}  Array of prediction of edges. The array is 3 dimensional (2 first dimension for storing a results, third dimension for compiling all results)

\code{- prediction} A dataframe with 3 columns (prediction, label and method) used for methods evaluation
}

\section{Methods}{

\code{$new(X.list=NULL,adjmat=NULL,nb.missing.var=0,methods="glasso",fit.number=20)} Initialize the experiment

\code{$run(bagging=FALSE)} Running the experiment (with or without using bagging)

\code{$roc.plot()} Plot the roc curves (all methods on the same plot)
}

\examples{
\dontrun{
star.graph <- graphModel$new(type = "starerdos",size=30, p.or.m = 0.05)
star.model <- GGMmodel$new(graph=star.graph)
plot(star.model)
star.model$randomSample(n=50)
testingGlasso<-GGMexperiment$new(X.list = list(star.model$getX()), adjmat = star.model$getAdjmat())
testingGlasso$run()
print(glasso.auc<-testingGlasso$auc())
testingGlasso$roc.plot()
}

}
