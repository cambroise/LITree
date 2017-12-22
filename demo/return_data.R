X<-load("../data/sp100return.Rda")
testingEM.latent.tree <- GGM.fit$new(X,method="em.latent.trees",nb.missing.var=1)
testingEM.latent.tree$run()


