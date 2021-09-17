# Kinase Activity
setwd("~/Dropbox/Columbia/GitHub/MAKINA/")

load("RData/CPTAC_GBM_netReg_kinPPtab.RData", verbose = T)
netReg <- netReg[sapply(netReg, length) >= 10]
kinPPtab <- kinPPtab[kinPPtab$kinase %in% names(netReg), ]

load("RData/CPTAC_GBM_baggingNet.RData", verbose = T)
load("RData/CPTAC_GBM_phosphoProteomics.RData", verbose = T)

load("RData/CPTAC_GBM_groups.RData", verbose = T)
cs <- intersect(names(groups), colnames(Phospho_Proteome))
groups <- groups[cs]
Phospho_Proteome <- Phospho_Proteome[, cs]
Phospho_Proteome <- t(apply(Phospho_Proteome, 1, function(x) (x - mean(x))/sd(x)))

means <- rowMeans(Phospho_Proteome)
means <- sort(means, decreasing = T)

library(Hmisc)
nBins <- 25
(m <- round(nrow(Phospho_Proteome)/nBins))
cuts <- as.numeric(cut2(means, m = m))
names(cuts) <- names(means)

nNullGenes <- 100
set.seed(1234)
nullGenes <- lapply(netReg, function(kin) unique(unlist(lapply(kin, function(gg) sample(names(cuts)[cuts == cuts[gg]], size = nNullGenes)))))

activityMat <- matrix(0, nrow = length(netReg), ncol = ncol(Phospho_Proteome))
rownames(activityMat) <- names(netReg)
colnames(activityMat) <- colnames(Phospho_Proteome)
scoresMat <- nullScoresMat <- activityMat

for(ss in 1:ncol(Phospho_Proteome)){
  scores <- sapply(1:length(netReg), function(x){
    kin <- names(netReg)[x]
    ta <- kinPPtab[kinPPtab$kinase == kin, ]
    rownames(ta) <- ta$siteID
    res <- weighted.mean(Phospho_Proteome[netReg[[kin]], ss], ta[netReg[[kin]], "score"])
    return(res)
  })
  
  nullScores <- sapply(1:length(nullGenes), function(x){
    kin <- names(nullGenes)[x]
    res <- weighted.mean(Phospho_Proteome[nullGenes[[kin]], ss], baggingNet[kin, nullGenes[[kin]]])
    return(res)
  })
  
  scoresMat[, ss] <- scores
  nullScoresMat[, ss] <- nullScores
  activityMat[, ss] <- scores - nullScores
  
  print(ss)
}
# save(activityMat, scoresMat, nullScoresMat, file = "RData/CPTAC_GBM_activityMat.RData")

load("RData/CPTAC_GBM_activityMat.RData", verbose = T)

DEAgroups <- function(ddata, groups){
  gg <- sort(unique(groups))
  
  ans <- vector("list", length(gg))
  names(ans) <- gg
  for(i in 1:length(gg)){
    whichOfInterest <- names(groups)[groups == gg[i]]
    theOthers <- setdiff(colnames(ddata), whichOfInterest)
    diffActivity <- apply(ddata, 1, function(x){
      suppressWarnings(a <- wilcox.test(x[whichOfInterest], x[theOthers]))
      xx <- as.numeric(a$statistic)/length(whichOfInterest)/length(theOthers)
      xx <- log2(xx/(1-xx))
      a <- c(xx, a$p.value)
      return(a)
    })
    diffActivity <- t(diffActivity)
    colnames(diffActivity) <- c("MWWscore", "pValue")
    fc <- rowMeans(ddata[, whichOfInterest]) - rowMeans(ddata[, theOthers])
    qValue <- p.adjust(diffActivity[, "pValue"], method = "fdr")
    
    diffActivity <- data.frame(statistic = diffActivity[, "MWWscore"], dm = fc, p.value = diffActivity[, "pValue"], fdr = qValue)
    colnames(diffActivity) <- c("MWWscore", "logFC", "pValue", "qValue")
    ans[[i]] <- diffActivity
  }
  return(ans)
}
aDEA <- DEAgroups(ddata = activityMat, groups = groups)

ttab <- table(kinPPtab$kinase, kinPPtab$Int)

aDEA <- lapply(aDEA, function(x){
  tmp <- ttab[rownames(x), , drop = F]
  tmp <- as.data.frame.matrix(tmp)
  colnames(tmp) <- c("TotInteractions", "KnownInteractions")
  tmp$TotInteractions <- rowSums(tmp)
  x <- cbind(x, tmp)
  return(x)
})

# save(aDEA, file = "RData/CPTAC_GBM_aDEA_activityMat.RData")

load("RData/CPTAC_GBM_aDEA_activityMat.RData", verbose = T)
aDEAfilt <- lapply(aDEA, function(x){
  x <- x[x$logFC > 0.3 & x$pValue < 0.01, ]
  x <- x[order(x$pValue), ]
  return(x)
})

toRemove <- table(unlist(lapply(aDEAfilt, rownames)))
toRemove <- toRemove[toRemove > 1]
toRemove
# named integer(0)