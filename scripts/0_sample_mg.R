#--------------------------------LOAD PACKAGES REQUIRED FOR ANALYSIS----------------------------------------------------
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfi)
library(Gviz)
library(DMRcate)
library(stringr)
library(dplyr)
library("IlluminaHumanMethylationEPICv2manifest")
library("IlluminaHumanMethylationEPICv2anno.20a1.hg38")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")


#-----------------------------------------READ ANNOTATIONS----------------------------------------------------------

annEPIC = getAnnotation(jokergoo / IlluminaHumanMethylationEPICv2anno.20a1.hg38)
head(annEPIC)

#---------------------------------------READ SAMPLES SHEET----------------------------------------------------------

BaseDir <- "C:/Users/n9652906/OneDrive - Queensland University of Technology/Desktop/ReNcell Arrays/ReNcell_Array_1/205676390096"
targets <- read.csv(
  file.path(
    "C:/Users/n9652906/OneDrive - Queensland University of Technology/Desktop/ReNcell Arrays/SampleSheet.csv"
  ),
  stringsAsFactors = FALSE,
  skip = 7
)
targets$Basename <- file.path(
  "C:/Users/n9652906/OneDrive - Queensland University of Technology/Desktop/ReNcell Arrays/SampleSheet.csv",
  targets$Sentrix_ID,
  paste(targets$Sentrix_ID, targets$Sentrix_Position, sep =
          "_")
)
targets$ID <- paste(targets$Sample_Group)
RGset <- read.metharray.exp(base = BaseDir, targets = targets)


# ### if error:
# Error in read.metharray(basenames = files, extended = extended, verbose = verbose,  :
#                           !anyDuplicated(basenames) is not TRUE
#
#                         Check SampleSheet - chances are the Sentrix ID has been lost to scientific numbering, or there are extra spaces recognised
#
#                         ---------------- NO NEED FOR TARGETS ID BUT JUST SET EASIER NAMES-----------------------------------------------------------------------------------

names <- targets$Sample_Names
celltype <- targets$Sample_Source
dosing <- targets$Sample_Conditions
passage <- targets$Sample_Passage


detP <- detectionP(RGset)  #### Get Detection P values #####


#--------------------------EXCLUDING BAD QUALITY SAMPLES------------------------------------------------------------------------------------------

keep <- colMeans(detP) < 0.01
RGset <- RGset[, keep]
targets <- targets[keep, ]
targets[, 1:5]
detP <- detP[, keep]
dim(detP)

#----------------------------------NORMALISATION-------------------------------------------------
mSetSq <- preprocessQuantile(RGset)
mSetFun <- preprocessFunnorm(RGset)
mSetRaw <- preprocessRaw(RGset)

#### mSetRaw is genuinely there just for plotting
#### Quantile if not so different tissue types, Funnorm if very different tissue types.
#### Sometimes quantile doesnt work for plotting (NA probes in some cases) so redo with Funnorm


###### MAKE SURE TO ADD PALETTES

pal2 <- brewer.pal(10, "Paired")


#------------------------------------------FILTERING SQ-----------------------------------------------------------------
###### CHANGE mSetSq to mSetFun if different normalisation occured ####################

detP <- detP[match(featureNames(mSetFun), rownames(detP)), ]
keep <- rowSums(detP < 0.01) == ncol(mSetFun)
table(keep)
mSetSqFlt <- mSetSq[keep, ]
mSetSqFlt

mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

######## CAN CHANGE MAF IF NEEDED ###############
#This code is from Miles

epic.cross1 <- read.csv(
  file.path(
    "C:/Users/n9652906/OneDrive - Queensland University of Technology/Desktop/ReNcell Arrays/ReNcell_Array_1/205676390096",
    pattern = "13059_2016_1066_MOESM1_ESM.csv"
  )
)
epic.cross2 <- read.csv(
  file.path(
    "C:/Users/n9652906/OneDrive - Queensland University of Technology/Desktop/ReNcell Arrays/ReNcell_Array_1/205676390096",
    pattern = "13059_2016_1066_MOESM2_ESM.csv"
  )
)
epic.cross3 <- read.csv(
  file.path(
    "C:/Users/n9652906/OneDrive - Queensland University of Technology/Desktop/ReNcell Arrays/ReNcell_Array_1/205676390096",
    pattern = "13059_2016_1066_MOESM3_ESM.csv"
  )
)
epic.variants1 <- read.csv(
  file.path(
    "C:/Users/n9652906/OneDrive - Queensland University of Technology/Desktop/ReNcell Arrays/ReNcell_Array_1/205676390096",
    pattern = "13059_2016_1066_MOESM4_ESM.csv"
  )
)
epic.variants2 <- read.csv(
  file.path(
    "C:/Users/n9652906/OneDrive - Queensland University of Technology/Desktop/ReNcell Arrays/ReNcell_Array_1/205676390096",
    pattern = "13059_2016_1066_MOESM5_ESM.csv"
  )
)
epic.variants3 <- read.csv(
  file.path(
    "C:/Users/n9652906/OneDrive - Queensland University of Technology/Desktop/ReNcell Arrays/ReNcell_Array_1/205676390096",
    pattern = "13059_2016_1066_MOESM6_ESM.csv"
  )
)
epic.add.probes <- c(
  as.character(epic.cross1$X),
  as.character(epic.variants1$PROBE),
  as.character(epic.variants2$PROBE),
  as.character(epic.variants3$PROBE)
)
epic.add.probes <- unique(epic.add.probes)
keep <- !(featureNames(mSetSqFlt) %in% epic.add.probes)
table(keep)

mSetSqFlt <- mSetSqFlt[keep, ]
mSetSqFlt

#----------------------------------Getting Beta Values for each CpG---------------------------------------------------------------

mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)