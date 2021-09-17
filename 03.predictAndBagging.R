setwd("~/Dropbox/Columbia/GitHub/MAKINA/")

load("RData/CPTAC_GBM_proteomics.RData", verbose = T)
load("RData/CPTAC_GBM_phosphoProteomics.RData", verbose = T)
load("RData/known_kinases_substrates.RData", verbose = T)

cs <- intersect(colnames(Proteome), colnames(Phospho_Proteome))
Proteome <- Proteome[, cs]
Phospho_Proteome <- Phospho_Proteome[, cs]
identical(colnames(Proteome), colnames(Phospho_Proteome))
# [1] TRUE

ck <- intersect(rownames(Proteome), kinaseTab$KINASE)
kinaseTab <- kinaseTab[kinaseTab$KINASE %in% ck,]
Proteome <- Proteome[ck, ]

sitesToCheck <- paste(kinaseTab$SUB_GENE, kinaseTab$SUB_MOD_RSD, sep = "-")

toChange <- unlist(lapply(strsplit(Phospho_Proteome_fData$site, "s"), function(x) x[[1]]))
toChange <- unlist(lapply(strsplit(toChange, "t"), function(x) x[[1]]))
toChange <- unlist(lapply(strsplit(toChange, "y"), function(x) x[[1]]))
Phospho_Proteome_fData <- cbind(Phospho_Proteome_fData, site2 = toChange, stringsAsFactors = F)

kinaseTab <- cbind(kinaseTab, site2 = paste0(kinaseTab$SUB_GENE, "-", kinaseTab$SUB_MOD_RSD), stringsAsFactors = F)

# Predict unknown interactions for each boot
library(doMC)
library(kernlab)
dir.create("RData/ksvmFit")

boot <- 100

lf <- paste0("RData/ksvmFit/negativeKinaseTab_ksvmFit_", 1:boot, ".RData")
for(nBoot in 1:length(lf)){
  load(lf[nBoot])

  nCores <- 100
  registerDoMC(nCores)
  kinPPnet <- foreach(i = 1:nrow(Proteome), .combine = "rbind") %dopar% {
    PP <- Phospho_Proteome
    colnames(PP) <- paste0(colnames(PP), "_PP")

    Pr <- Proteome[rep(rownames(Proteome)[i], nrow(PP)), ]
    colnames(Pr) <- paste0(colnames(Pr), "_Pr")

    testMat <- cbind(PP, Pr)

    tmp <- predict(ksvmFit, testMat, type = "probabilities")
    tmp <- tmp[, "1"]

    # print(paste(nBoot, i))
    return(tmp)
  }
  rownames(kinPPnet) <- rownames(Proteome)
  colnames(kinPPnet) <- rownames(Phospho_Proteome)
  ffile <- paste0("RData/kinPPnet/kinPPnet_", nBoot, ".RData")
  save(kinPPnet, file = ffile)
  message(nBoot)
}

# create bagging net
boot <- 100

lf <- paste0("RData/kinPPnet/kinPPnet_", 1:boot, ".RData")
load(lf[1], verbose = T)
baggingNet <- kinPPnet

for(i in 2:length(lf)){
  load(lf[i])
  baggingNet <- baggingNet + kinPPnet
  print(i)
}
baggingNet <- baggingNet/boot
# save(baggingNet, file = "RData/CPTAC_GBM_baggingNet.RData")