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

# Train SVM classifier from known positive and random negative interactions for each boot
library(doMC)
library(kernlab)
dir.create("RData/ksvmFit")

boot <- 100
nCores <- 100
registerDoMC(nCores)
foreach(nBoot = 1:boot) %dopar% {
  negativeKinaseTab <- NULL
  kinase <- sort(unique(kinaseTab$KINASE))
  for(i in 1:length(kinase)){
    tmpKT <- kinaseTab[kinaseTab$KINASE == kinase[i], ]
    nNeg <- nrow(tmpKT)
    if(nNeg < 10) nNeg <- 10
    toChoose <- Phospho_Proteome_fData[!Phospho_Proteome_fData$siteID %in% tmpKT$siteID, ]
    toChoose <- toChoose[sample(rownames(toChoose), size = nNeg), ]

    KINASE <- rep(kinase[i], nNeg)
    SUB_GENE <- sapply(strsplit(toChoose$site2, "-"), function(x) x[1])
    SUB_MOD_RSD <- sapply(strsplit(toChoose$site2, "-"), function(x) x[2])
    site2 <- toChoose$site2
    site <- toChoose$site
    siteID <- toChoose$siteID
    toAdd <- data.frame(KINASE = KINASE, SUB_GENE = SUB_GENE, SUB_MOD_RSD = SUB_MOD_RSD,
                        site2 = site2, site = site, siteID = siteID, stringsAsFactors = F)

    negativeKinaseTab <- rbind(negativeKinaseTab, toAdd)
  }
  negativeKinaseTab <- cbind(negativeKinaseTab, Int = 0)

  kt <- rbind(kinaseTab, negativeKinaseTab)
  rownames(kt) <- paste0("Int_", 1:nrow(kt))

  PP <- Phospho_Proteome[kt$siteID, ]
  rownames(PP) <- rownames(kt)
  colnames(PP) <- paste0(colnames(PP), "_PP")

  Pr <- Proteome[kt$KINASE, ]
  rownames(Pr) <- rownames(kt)
  colnames(Pr) <- paste0(colnames(Pr), "_Pr")

  trainMat <- cbind(PP, Pr)

  groups <- kt$Int
  names(groups) <- rownames(kt)
  groups <- factor(groups)

  ksvmFit <- ksvm(x = trainMat, y = groups, prob.model = T, kernel = "rbfdot", C = 8)

  ffile <- paste0("RData/ksvmFit/negativeKinaseTab_ksvmFit_", nBoot, ".RData")
  save(negativeKinaseTab, ksvmFit, file = ffile)
  print(nBoot)
}