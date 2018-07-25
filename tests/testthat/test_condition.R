test_that("ConditionBadCichlid", {
  utils::data("cichlid")
  parameters.to.try <- c(sigma.sq=6, mu=mean(cichlid$data)) #terrible values
  expect_warning(calculate.likelihood.result <- CalculateLikelihood(x=parameters.to.try , data=cichlid$data, phy=cichlid$phy, flow=cichlid$flow, allow.extrapolation=FALSE, do.kappa.check=TRUE, measurement.error=0))
  #V.modified <- GetVModified(x=x, phy=cichlid$phy, flow=cichlid$flow, actual.params=free.parameters[which(free.parameters)], measurement.error=NULL)
#  matrix.condition <- kappa(V.modified, exact=TRUE) high
})

test_that("ConditionBadNicotiana", {
  utils::data("nicotiana")
  parameters.to.try <- c(sigma.sq=6, mu=mean(nicotiana$data)) #terrible values
  expect_warning(calculate.likelihood.result <- CalculateLikelihood(x=parameters.to.try, data=nicotiana$data, phy=nicotiana$phy, flow=nicotiana$flow, allow.extrapolation=FALSE, do.kappa.check=TRUE, measurement.error=0))
  #V.modified <- GetVModified(x=x, phy=cichlid$phy, flow=cichlid$flow, actual.params=free.parameters[which(free.parameters)], measurement.error=NULL)
#  matrix.condition <- kappa(V.modified, exact=TRUE) high
})

test_that("BasicRun",{
  utils::data("cichlid")
  result <- BMhyb(cichlid$data, cichlid$phy, cichlid$flow, n.points=100, get.se=TRUE, plot.se=FALSE, measurement.error=0,n.random.start.points=100)
  expect_equal(class(result), "data.frame")
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 22)
})

test_that("Issue 13 is solved", {
  create_paper_network <- function(gamma, t1, t2, t3){
      phy <- ape::read.tree(text = paste0("((R:", t3, ",Y:", t3, "):", t1 + t2, ",X:", t1 + t2 + t3, ");"))
      network <- list(phy = phy,
                      flow = data.frame(donor = "X",
                                        recipient = "R",
                                        gamma = gamma,
                                        time.from.root.donor = t1,
                                        time.from.root.recipient = t1 + t2))
      network$flow$donor <- as.character(network$flow$donor)
      network$flow$recipient <- as.character(network$flow$recipient)
      return(network)
  }

  gamma <- 0.5
  t1 <- 0.3; t2 <- 0.4; t3 <- 0.3; # unit height
  network <- create_paper_network(gamma, t1, t2, t3)

  #PlotNetwork(network$phy, network$flow)
  #axis(1, at = c(0, t1, t1+t2, t1+t2+t3), labels = c("0", "t1", "t1+t2", "t1+t2+t3"))

  sigma2 = 1
  x <- c(sigma.sq = sigma2, mu = 0, SE = 0)
  vcv_BMhyb <- GetVModified(x, network$phy, network$flow, measurement.error=0)
  expect_equal(vcv_BMhyb, structure(c(0.65, 0.7, 0.15, 0.7, 1, 0, 0.15, 0, 1), .Dim = c(3L,
+ 3L), .Dimnames = list(c("R", "Y", "X"), c("R", "Y", "X")))) #from Issue 13
})

test_that("Issue 14 is solved", {
  gamma <- 0.5
  t1 <- 0.3; t2 <- 0.4; t3 <- 0.3;
  phy <- ape::read.tree(text = paste0("((R:", t3, ",Y:", t3, "):", t1 + t2, ",X:", t1 + t2 + t3, ");"))
  ## Network
  don_recp <- expand.grid(c("X"), c("Y", "R"))
  network <- list(phy = phy,
                  flow = data.frame(donor = don_recp[,1],
                                    recipient = don_recp[,2],
                                    gamma = rep(gamma, 2),
                                    time.from.root.donor = rep(t1, 2),
                                    time.from.root.recipient = rep(t1, 2)))
  network$flow$donor <- as.character(network$flow$donor)
  network$flow$recipient <- as.character(network$flow$recipient)
  ## Plot

  sigma2 = 1
  x <- c(sigma.sq = sigma2, mu = 0, SE = 0)
  expect_equal(GetVModified(x, network$phy, network$flow, measurement.error=0), structure(c(0.85, 0.55, 0.15, 0.55, 0.85, 0.15, 0.15, 0.15, 1
), .Dim = c(3L, 3L), .Dimnames = list(c("R", "Y", "X"), c("R",
"Y", "X"))))
})

test_that("VH matters", {
  gamma <- 0.5
  t1 <- 0.3; t2 <- 0.4; t3 <- 0.3;
  phy <- ape::read.tree(text = paste0("((R:", t3, ",Y:", t3, "):", t1 + t2, ",X:", t1 + t2 + t3, ");"))
  ## Network
  don_recp <- expand.grid(c("X"), c("Y", "R"))
  network <- list(phy = phy,
                  flow = data.frame(donor = don_recp[,1],
                                    recipient = don_recp[,2],
                                    gamma = rep(gamma, 2),
                                    time.from.root.donor = rep(t1, 2),
                                    time.from.root.recipient = rep(t1, 2)))
  network$flow$donor <- as.character(network$flow$donor)
  network$flow$recipient <- as.character(network$flow$recipient)
  ## Plot

  sigma2 = 1
  x <- c(sigma.sq = sigma2, mu = 0, vh=0, SE=0)


  V0 <- GetVModified(x, network$phy, network$flow, measurement.error=0)
  vh.add = 3
  y <- c(sigma.sq = sigma2, mu = 0, vh=vh.add,SE = 0)
  V1 <- GetVModified(y, network$phy, network$flow, measurement.error=0)
  expect_equal(V1[1,1]-vh.add, V0[1,1])
  expect_equal(V1[3,3], V0[3,3])
})

#
# test_that("MergingTrees", {
#   # idea is a network, decomposed into trees -- does it come back
#   phy1 <- ape::read.tree(text="(((A:2,B:2):1,C:3):5,D:8);")
#   phy2 <- ape::read.tree(text="(((B:2,C:2):1,A:3):5,D:8);")
#   phy.graph <- MergeTreesIntoPhyGraph(c(phy1, phy2))
#   free.parameters<-rep(TRUE, 5)
#   names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
#   free.parameters[which(names(free.parameters)=="bt")]<-FALSE
#   free.parameters[which(names(free.parameters)=="vh")]<-FALSE
#   x <- c(1, 0, 1)
#   results <- GetVandMFromIgraph(x, phy.graph, free.parameters)
#   expect_equal(nrow(results$V.modified), 4)
# })
#
# test_that("NonNetworkWorks", {
#   # idea is a network, decomposed into trees -- does it come back
#   phy1 <- ape::read.tree(text="(((A:2,B:2):1,C:3):5,D:8);")
#   phy2 <- phy1
#   phy.graph <- MergeTreesIntoPhyGraph(c(phy1, phy2))
#   free.parameters<-rep(TRUE, 5)
#   names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
#   free.parameters[which(names(free.parameters)=="bt")]<-FALSE
#   free.parameters[which(names(free.parameters)=="vh")]<-FALSE
#   x <- c(1, 0, 1)
#   results <- GetVandMFromIgraph(x, phy.graph, free.parameters)
#   expect_lte(max(abs(results$V.modified - ape::vcv(phy1))),1e-8)
# })
