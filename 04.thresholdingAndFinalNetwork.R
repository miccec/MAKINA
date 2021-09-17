# Threshold 50% known interactions kept (weakly supervised learning)
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

load("RData/kinaseTab_fromPSP.RData", verbose = T)
Proteome <- Proteome[sort(unique(kinaseTab$KINASE)), ]

load("RData/CPTAC_GBM_baggingNet.RData", verbose = T)

kinaseTab <- cbind(kinaseTab, score = NA, stringsAsFactors = F)
for(i in 1:nrow(kinaseTab)) kinaseTab[i, "score"] <- baggingNet[kinaseTab$KINASE[i], kinaseTab$siteID[i]]

plot(density(as.vector(baggingNet)))
summary(as.vector(baggingNet))
rug(kinaseTab$score)
summary(kinaseTab$score)
(thresh <- quantile(kinaseTab$score, probs = 0.50))
abline(v=thresh, col = "red")

kinPPnet <- matrix(0, nrow = nrow(baggingNet), ncol = ncol(baggingNet))
rownames(kinPPnet) <- rownames(baggingNet)
colnames(kinPPnet) <- colnames(baggingNet)
kinPPnet[baggingNet >= thresh] <- 1

rowSums(kinPPnet != 0)
range(rowSums(kinPPnet != 0))
# save(kinPPnet, file = "RData/CPTAC_GBM_kinPPnet_50perc.RData")

# Create regulons with kinase-phosphosite cor > 0
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

load("RData/CPTAC_GBM_baggingNet.RData", verbose = T)
load("RData/CPTAC_GBM_kinPPnet_50perc.RData", verbose = T)
load("RData/kinaseTab_fromPSP.RData", verbose = T)

kinPPnet <- kinPPnet[rowSums(kinPPnet != 0) > 0, ]
baggingNet <- baggingNet[rownames(kinPPnet), ]

Proteome <- Proteome[rownames(kinPPnet), ]
Phospho_Proteome <- Phospho_Proteome[colnames(kinPPnet), ]
Phospho_Proteome_fData <- Phospho_Proteome_fData[colnames(kinPPnet), ]

# library(psych)
# tmp <- corr.test(t(Proteome), t(Phospho_Proteome), method = "spearman", adjust = "none", ci = FALSE)
# ccor <- tmp$r
# ccor_pVal <- tmp$p
# # save(ccor, ccor_pVal, file = "RData/CPTAC_GBM_kinPPnet_50perc_ccor.RData")
load("RData/CPTAC_GBM_kinPPnet_50perc_ccor.RData", verbose = T)

netReg <- vector("list", nrow(kinPPnet))
names(netReg) <- rownames(kinPPnet)
kinPPtab <- NULL
for(i in 1:nrow(kinPPnet)){
  toTake <- kinPPnet[i, ]
  toTake <- toTake[toTake != 0]
  toTake <- names(toTake)

  score <- baggingNet[i, toTake]
  score <- sort(score, decreasing = T)

  if(length(score) == 1) siteID <- toTake else
    siteID <- names(score)
  Phoshosite <- Phospho_Proteome_fData[siteID, "site"]
  kinase <- rep(rownames(kinPPnet)[i], length(siteID))
  # score
  rho <- ccor[i, siteID]

  Int <- rep(0, length(siteID))
  tab <- data.frame(siteID = siteID, Phoshosite = Phoshosite, kinase =  kinase, score = score, rho = rho, Int = Int, stringsAsFactors = F)
  rownames(tab) <- tab$siteID

  tmp <- kinaseTab[kinaseTab$KINASE == rownames(kinPPnet)[i], ]
  tmp <- tmp[tmp$siteID %in% siteID, ]
  if(nrow(tmp) > 0) tab[unique(tmp$siteID), "Int"] <- 1
  rownames(tab) <- NULL
  tab <- tab[tab$rho > 0, ]

  kinPPtab <- rbind(kinPPtab, tab)

  netReg[[i]] <- tab$siteID

  print(i)
}

sapply(netReg, length)
sort(sapply(netReg, length))

# save(netReg, kinPPtab, file = "RData/CPTAC_GBM_netReg_kinPPtab.RData")