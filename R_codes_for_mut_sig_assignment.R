### R codes for mutational signature assignment (SBS96 and ID)
### Courtesy of Nanhai Jiang for the help on the script
rm(list = ls())

# Install dependency packages and load in data
library(remotes)
library(dplyr)
library(data.table)
library(parallel)

# to install mSigAct R package, please use: 
# remotes::install_github(repo = "steverozen/mSigAct", ref = "v2.3.2-branch")
library(mSigAct) 

# to install ICAMS R package, please use: 
# remotes::install_github(repo = "steverozen/ICAMS", ref = "v3.0.6-branch")
library(ICAMS) 

# to install cosmicsig R package, please use: 
# remotes::install_github(repo = "Rozen-Lab/cosmicsig", ref = "v1.0.7-branch")
library(cosmicsig) 

download.file("https://github.com/CYHuang-Lab/CGMH-OCCC-WES-project/blob/main/CGMH-OCCC-mut-sig-catalogs.rds")
catalogs <- readRDS("CGMH-OCCC-mut-sig-catalogs.rds")
catSBS96 <- catalogs$catSBS96
colnames(catSBS96) <- sub(".intersect", "", colnames(catSBS96))

sigs_sbs96_genome <- cosmicsig::COSMIC_v3.2$signature$GRCh38$SBS96

sigs_sbs96_exome <- ICAMS::TransformCatalog(
  catalog = sigs_sbs96_genome, target.region = "exome", target.ref.genome = "GRCh37"
)

sigs_prop_sbs96 <- mSigAct::ExposureProportions(
  mutation.type = "SBS96", cancer.type = "Ovary-AdenoCA"
) # also use PCAWG7::CancerTypes() to see other cancer types

sbs96.base <- paste0("SBS",c(1,2,5,8,13,18,39,40,41))
sbs96.full <- paste0("SBS",c(1,2,5,8,13,18,39,40,41,3,6,26,44))
sigs_sbs96_exome <- sigs_sbs96_exome[, sbs96.full, drop = FALSE]
my_opts <- mSigAct::DefaultManyOpts(likelihood.dist = "multinom")

## Include SBS3 (HRD), SBS6 (dMMR), SBS26 (dMMR), SBS44 (dMMR)
## into the signature universe by indel burden
sig.table <- data.frame(
  ID = colnames(catSBS96), 
  MH_5 = colSums(catalogs$catID[grep("MH:5",rownames(catalogs$catID)),]),
  Homopolymer_T = colSums(catalogs$catID[grep("T:1:5",rownames(catalogs$catID)),])
)

sigs_to_use <- list()

for(i in 1:nrow(sig.table)){
  sigs_to_use[[i]] <- sbs96.base
  if(sig.table$MH_5[i] > 5){sigs_to_use[[i]] <- c(sigs_to_use[[i]], "SBS3")} 
  if(sig.table$Homopolymer_T[i] > 20){sigs_to_use[[i]] <- c(sigs_to_use[[i]],"SBS6","SBS26","SBS44")} 
}

retval <- lapply(1:nrow(sig.table), function(x){
  mSigAct::SparseAssignActivity(
    spectra = catSBS96[,x,drop=FALSE],
    sigs = sigs_sbs96_exome[,sigs_to_use[[x]]],
    output.dir = "./mutsig_new/SBS_sparse/",
    max.level = length(sigs_to_use[[x]]) - 1,
    p.thresh = 0.05 / length(sigs_to_use[[x]]),
    m.opts = my_opts,
    num.parallel.samples = 1,
    mc.cores.per.sample = 1,
    seed = 3838,
    max.subsets = 1e15,
    drop.low.mut.samples = FALSE
  )
})

saveRDS(retval, "./mutsig_new/SBS_sparse/sparse.out.SBS96.rds")


## mutational signature assignment for indels
catID <- catalogs$catID
colnames(catID) <- sub(".intersect", "", colnames(catID))
sigs.ID <- PCAWG7::signature$genome$ID
sigs.prop.ID <- mSigAct::ExposureProportions(
  mutation.type = "ID", cancer.type = "Ovary-AdenoCA"
)
sigs.ID.to.use <- c(names(sigs.prop.ID), "ID7") #add ID7 for dMMR tumor

ID_out <- mSigAct::SparseAssignActivity(
  spectra = catID,
  sigs = sigs.ID[,sigs.ID.to.use],
  output.dir = "./mutsig/ID",
  max.level = length(sigs.ID.to.use) - 1,
  m.opts = my_opts,
  p.thresh = 0.05 / length(sigs.ID.to.use),
  seed = 8787,
  max.subsets = 1e15,
  num.parallel.samples = 1,
  mc.cores.per.sample = 1,
  drop.low.mut.samples = FALSE
)
saveRDS(ID_out, "./mutsig/ID/sparse.out.ID.rds")
