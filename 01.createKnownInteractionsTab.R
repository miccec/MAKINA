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

x <- sapply(kinaseTab$site2, function(x) sum(Phospho_Proteome_fData$site2 %in% x))
table(x)
#     0     1     2     3     4     5     6     7     8    10    11
# 16214  1863   548   148    44    22     4     5     4     2     2

# Create known interactions table
kinaseTabNew <- NULL
for(i in 1:nrow(kinaseTab)){
  idx <- which(Phospho_Proteome_fData$site2 == kinaseTab$site2[i])
  if(length(idx) == 0) next
  print(paste(i, "of", nrow(kinaseTab)))
  tmp <- Phospho_Proteome_fData[idx, ]
  KINASE <- rep(kinaseTab$KINASE[i], nrow(tmp))
  SUB_GENE <- rep(kinaseTab$SUB_GENE[i], nrow(tmp))
  SUB_MOD_RSD <- rep(kinaseTab$SUB_MOD_RSD[i], nrow(tmp))
  site2 <- rep(kinaseTab$site2[i], nrow(tmp))
  site <- tmp$site
  siteID <- tmp$siteID
  toAdd <- data.frame(KINASE = KINASE, SUB_GENE = SUB_GENE, SUB_MOD_RSD = SUB_MOD_RSD,
                      site2 = site2, site = site, siteID = siteID, stringsAsFactors = F)

  kinaseTabNew <- rbind(kinaseTabNew, toAdd)
}
kinaseTab <- kinaseTabNew
kinaseTab <- cbind(kinaseTab, Int = 1)
# save(kinaseTab, file = "RData/kinaseTab_fromPSP.RData")