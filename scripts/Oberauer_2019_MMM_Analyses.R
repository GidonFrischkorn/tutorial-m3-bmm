#### Read frequency data, fit with MMM


graphics.off()
rm(list=ls(all=TRUE))

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

require(R2jags)
source("plotPostKO.R")
source("plot.credint.R")
source("HDIofMCMC.R")
source("lineplot.ci3.R")

experiment <- 2
model <- 2
plotoption <- 2 # 0: Plots window; 1: R Graphics device; 2: PDF
traceP <- T  # trace predicted probabilities
traceS <- T  # trace subject parameters

modelname <- paste0("PairsRSS.MMM", model, ".txt")
data <- read.table(paste0("PairsRSS", experiment, "_freq.dat"), header=T)
aggdat <- aggregate(fcorrect ~ setsize+setsize+rsizeList+rsizeNPL, data=data, FUN=mean)
nconds <- dim(aggdat)[1]
nsubj <- length(unique(data$id))
N <- unique(data$N)
Setsizes <- sort(unique(aggdat$setsize))
RSSlist <- sort(unique(aggdat$rsizeList))
RSSnpl <- sort(unique(aggdat$rsizeNPL))

Fdata <- array(NA, dim=c(nconds, nsubj, 3))
SSidx <- NA   # consecutive indices for set sizes (not actual set sizes!)
RSSidx <- matrix(NA, nconds, 2)  # consecutive indices for RSS list / RSS new

for (cond in 1:nconds) {
  subdat <- subset(data, setsize==aggdat[cond, "setsize"] & rsizeList==aggdat[cond, "rsizeList"] & rsizeNPL==aggdat[cond, "rsizeNPL"])
  Fdata[cond,,1] <- subdat[,"fcorrect"]
  Fdata[cond,,2] <- subdat[,"fother"]
  Fdata[cond,,3] <- subdat[,"fnpl"]
  SSidx[cond] <- which(Setsizes == aggdat[cond, "setsize"])
  RSSidx[cond,1] <- which(RSSlist == aggdat[cond, "rsizeList"])
  RSSidx[cond,2] <- which(RSSnpl == aggdat[cond, "rsizeNPL"])
}

          
muparameters <- c("muC", "muA")
sgparameters <- c("sgC", "sgA")
if (model == 2) muparameters <- c(muparameters, "muB")
if (model == 2) sgparameters <- c(sgparameters, "sgB")
subjectparameters <- c("C", "A")
parameters <- c(muparameters, sgparameters, "loglik")
if (model == 2) parameters <- c(parameters, "meanA", "meanC", "meanB", "dA", "dC")
if (traceS == T) parameters <- c(parameters, subjectparameters)
if (traceP == T) parameters <- c(parameters, "P")
parnames <- tolower(subjectparameters)

dataList = list(k = Fdata, 
                N = N,
                nsubj = nsubj,
                nconds = nconds,
                nsetsize = length(unique(aggdat$setsize)),
                SSidx = SSidx,
                RSSidx = RSSidx,
                Setsizes = Setsizes,
                RSSlist = RSSlist,
                RSSnpl = RSSnpl
              )
if (model == 2) dataList[["Csetsize"]] <- sort(unique(aggdat$setsize))-5 # centered set sizes

################################# RUN THE CHAINS  ##########################

jags2Model <- jags.parallel(data=dataList, inits=NULL, parameters.to.save=parameters, model.file = modelname,
                            n.chains = 3, n.iter = 10000, n.burnin = 5000,
                            n.thin = 1, n.cluster= 3)
codaSamples <- as.mcmc(jags2Model)
mcmcChain = as.matrix( codaSamples )


# Posterior group mean parameters

if (model == 2) {
  
  indivAC <- matrix(NA, nsubj, 2)
  for (id in 1:nsubj) {
    for (p in 1:2) {
      indivAC[id,p] <- mean(mcmcChain[,paste0(muparameters[p], "[", id, "]")])
    }
  }
  x11()
  plot(indivAC[,1], indivAC[,2], xlab="muC", ylab="muA")
    
  x11()
  layout(matrix(1:4, 2, 2, byrow=T))
  plotPostKO(mcmcChain[,"meanA"], xlab=bquote(paste(mu, "A")), breaks=30)
  plotPostKO(mcmcChain[,"meanB"], xlab=bquote(paste(mu, "B")), breaks=30)
  plotPostKO(mcmcChain[,"dA"], xlab=bquote(paste(delta, "A")), breaks=30)
  plotPostKO(mcmcChain[,"dC"], xlab=bquote(paste(delta, "C")), breaks=30)
  
  x11()
  layout(matrix(1:4, 2, 2, byrow=T))
  plotPostKO(mcmcChain[,"sgA"], xlab=bquote(paste(sigma, "A")), breaks=30)
  plotPostKO(mcmcChain[,"sgB"], xlab=bquote(paste(sigma, "B")), breaks=30)
  plotPostKO(mcmcChain[,"sgC"], xlab=bquote(paste(sigma, "C")), breaks=30)
}

#### Prepare data for joint plot of both experiments

if (model == 2) {
  
  parChainFile <- "PairsRSS.MMMparChains.RData"
  chainfilethere <- file.access(parChainFile, mode=0)
  if (chainfilethere==0) {
    load(parChainFile)
  } else {
    MuA <- array(0, dim=c(length(Setsizes), 2, dim(mcmcChain)[1]))
    MuC <- array(0, dim=c(length(Setsizes), 2, dim(mcmcChain)[1]))
    MuB <- array(0, dim=c(1, 2, dim(mcmcChain)[1]))
    meanA <- array(0, dim=c(length(Setsizes), 2, dim(mcmcChain)[1]))
    meanC <- array(0, dim=c(length(Setsizes), 2, dim(mcmcChain)[1]))
    DeltaA <- matrix(0, dim(mcmcChain)[1], 2)
    DeltaC <- matrix(0, dim(mcmcChain)[1], 2)
  }
  
  CSetsizes <- Setsizes - 5
  for (s in 1:length(Setsizes)) {
    MuA[s,experiment,] <- mcmcChain[, "meanA"] + mcmcChain[, "dA"] * CSetsizes[s]
    MuC[s,experiment,] <- mcmcChain[, "meanC"] + mcmcChain[, "dC"] * CSetsizes[s]
    # Posterior samples of group-level means of A and C for each set size, predicted by linear model 
    aIdx <- NULL
    for (j in 1:nsubj) aIdx <- c(aIdx, which(names(as.data.frame(mcmcChain)) == paste0("A[", j, ",", s, "]")))
    cIdx <- NULL
    for (j in 1:nsubj) cIdx <- c(cIdx, which(names(as.data.frame(mcmcChain)) == paste0("C[", j, ",", s, "]")))    
    meanA[s,experiment,] <- rowMeans(mcmcChain[, aIdx])
    meanC[s,experiment,] <- rowMeans(mcmcChain[, cIdx])
  }
  MuB[1,experiment,] <- mcmcChain[, "meanB"]  # Posterior samples of group-level mean of B
  DeltaA[,experiment] <- mcmcChain[, "dA"]    # Posterior samples of group-level mean of regression slope on A
  DeltaC[,experiment] <- mcmcChain[, "dC"]
  
  save(MuA, MuC, MuB, meanA, meanC, DeltaA, DeltaC, file=parChainFile)


  
}



############### Prepare Plot of Data and Predictions #######################################3

predFile <- paste0("PairsRSS.postpred", experiment, ".RData")
predfilethere <- file.access(predFile, mode=0)
if (predfilethere==0) {
  load(predFile)
} else {
  meanPropResponse <- list()
  CIupperPropResponse <- list()
  CIlowerPropResponse <- list()
  meanPostPred <- list()
  CIupperPostPred <- list()
  CIlowerPostPred <- list()
}

predP <- array(0,dim=c(nsubj,nconds,3))
predN <- array(0,dim=c(nsubj,nconds,3))
predProb <- array(0,dim=c(nsubj,nconds,3))
datP <- array(0,dim=c(nsubj,nconds,3))
datProb <- array(0,dim=c(nsubj,nconds,3))

for (idx in 1:nsubj) {
  for (c in 1:nconds) {
    for (r in 1:3) {
      condition <- 1
      Pchain <- mcmcChain[, paste("P[", as.character(c), ",", as.character(idx), ",", as.character(r), "]", sep="")]
      predP[idx,c,r] <- mean(Pchain)
      predN[idx,c,r] <- predP[idx,c,r]*N
      datP[idx,c,r] <- Fdata[c,idx,r]/N
    }
  }
}

respCatName <- c("Correct", "Other", "New")
for (respCat in 1:3) {
  
  # Bakeman-McArthur 95% confidence intervals

  pdat <- datP[,,respCat]   # pick out the response category (1=correct, 2=other, 3=NPL)  
  rmean = rowMeans(pdat)  # each subject's mean across all conditions
  pcdat <- pdat - rmean + mean(rmean)  
  ppred <- predP[,,respCat]
  rmean = rowMeans(ppred)
  pcpred <- ppred - rmean + mean(rmean)
  
  M.data <- array(0,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  M.pred <- array(NA,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  CI.data.u <- array(NA,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  CI.data.l <- array(NA,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  CI.pred.u <- array(NA,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  CI.pred.l <- array(NA,dim=c(4,5,4))  # RSS-npl, RSS-list, setsize
  
  for (c in 1:nconds) {
    M.data[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcdat[,c])
    std <- sd(pcdat[,c])
    CI.data.u[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcdat[,c]) + 1.96*std/sqrt(nsubj)    
    CI.data.l[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcdat[,c]) - 1.96*std/sqrt(nsubj)   
    M.pred[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcpred[,c])
    std <- sd(pcpred[,c])
    CI.pred.u[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcpred[,c]) + 1.96*std/sqrt(nsubj)    
    CI.pred.l[RSSidx[c,2], RSSidx[c,1], SSidx[c]] <- mean(pcpred[,c]) - 1.96*std/sqrt(nsubj)   
  }
  
  meanPropResponse[[respCat]] <- M.data
  CIupperPropResponse[[respCat]] <- CI.data.u
  CIlowerPropResponse[[respCat]] <- CI.data.l
  meanPostPred[[3+respCat]] <- M.pred
  CIupperPostPred[[3+respCat]] <- CI.pred.u
  CIlowerPostPred[[3+respCat]] <- CI.pred.l
  save(meanPropResponse, CIupperPropResponse, CIlowerPropResponse, meanPostPred, CIupperPostPred, CIlowerPostPred, file=predFile)

}


#### Prediction of recall accuracy

recallfile <- "PairsRSS.Recallpred.Rdata"
rfilethere <- file.access(recallfile, mode=0)
if (rfilethere==0) {
  load(recallfile)
} else {
  RecallPredList <- list()
}


nNPL <- 5
CSetsizes <- Setsizes - 5
PredRecall <- matrix(NA, nsubj, length(Setsizes))
PredMeanRecall <- matrix(NA, 3, length(Setsizes))
for (s in 1:length(Setsizes)) {
  for (subj in 1:nsubj) {
    A <- mcmcChain[, paste0("A[", subj, ",", s, "]")]
    C <- mcmcChain[, paste0("C[", subj, ",", s, "]")]
    if (model == 1) B <- 0.1
    if (model == 2) B <- mcmcChain[, paste0("muB[", subj, "]")]
    actCorr <- B + A + C
    actOther <- (B + A)*(Setsizes[s]-1)  # assuming all list elements are in effective response set  
    actNPL <- B * nNPL
    pcorrect <- actCorr/(actCorr + actOther + actNPL)
    PredRecall[subj, s] <- mean(pcorrect)
  }
  if (model == 1) {
    muA <- mcmcChain[, paste0("muA[", s, "]")]
    muC <- mcmcChain[, paste0("muC[", s, "]")]
  }
  if (model == 2) {
    muA <- mcmcChain[, paste0("meanA")] + mcmcChain[, paste0("dA")] * CSetsizes[s]
    muC <- mcmcChain[, paste0("meanC")] + mcmcChain[, paste0("dC")] * CSetsizes[s]
    muB <- mcmcChain[, "meanB"]
  }
  ActCorr <- muB + muA + muC
  ActOther <- (muB + muA)*(Setsizes[s]-1)  # assuming all list elements are in effective response set  
  ActNPL <- muB * nNPL
  Pcorrect <- ActCorr/(ActCorr + ActOther + ActNPL)
  PredMeanRecall[1,s] <- mean(Pcorrect)
  PredMeanRecall[2:3,s] <- HDIofMCMC(Pcorrect)
}

PredRecallCorr <- BakemanWide(PredRecall)
PredMeanCI <- ConfintWide(PredRecallCorr)
RecallPredList[[experiment]] <- PredMeanCI
save(RecallPredList, file=recallfile)

### Read recall data, plot data + predictions

if (plotoption == 1) x11()
if (plotoption == 2) pdf(paste0("RSS.MMM.pdf"))
layout(matrix(1:2,1,2))
for (experiment in 1:2) {
  pathname <- paste0("C:\\Daten\\Rohdaten\\Pairs.Binding\\PairsRSS", experiment)
  dat <- read.table(paste0(pathname, "\\PairsRSS", experiment, "_all.dat"), header=F)
  names(dat) <- c("id", "session", "block", "trial", "setsize", 
                  "rsizeList", "rsizeNPL", "tested", "response", "rcat", "rt")
  dat$rsize <- dat$rsizeList + dat$rsizeNPL
  dat$correct <- dat$rcat==1
  dat$intrusion <- dat$rcat==2
  dat$new <- dat$rcat==3
  
  recogdat <- subset(dat, rsize>0)
  recalldat <- subset(dat, rsize==0)
  recallagg <- aggregate(correct ~ id+setsize, data=recalldat, FUN=mean)
  PredMeanCI.MMM <- RecallPredList[[experiment]]  # prediction from MMM
  PredMeanCI.MPT <- RecallPredList[[2+experiment]]  # prediction from MPT
  
  lineplot.ci3(recallagg, dv="correct", iv="setsize", id="id", x=1, off=-0.1+0.03*7, xlim=c(0.5,8.5), ylim=c(0,1), 
               xlab="Set Size", ylab="P(correct)", pt=21, ptcol="black", cex=1.5)
  # errbar(x=sort(unique(recallagg$setsize))+0.2, y=PredMeanRecall[1,], yplus=PredMeanRecall[2,], yminus=PredMeanRecall[3,],
  #        add=T, col="red", errbar.col="red")
  errbar(x=sort(unique(recallagg$setsize))+0.1, y=PredMeanCI.MMM[1,], yplus=PredMeanCI.MMM[2,], yminus=PredMeanCI.MMM[3,],
         add=T, col="red", errbar.col="red")
  errbar(x=sort(unique(recallagg$setsize))+0.2, y=PredMeanCI.MPT[1,], yplus=PredMeanCI.MPT[2,], yminus=PredMeanCI.MPT[3,],
         add=T, col="blue", errbar.col="blue")
}
if (plotoption==2) dev.off()


