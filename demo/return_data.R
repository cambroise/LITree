X<-read.table("../data/sp100_returns.txt")
testingGlasso<-GGM.fit$new(X)
testingEM.latent.tree <- GGM.fit$new(X,method="em.latent.trees",nb.missing.var=1)
testingEM.latent.tree$run()
