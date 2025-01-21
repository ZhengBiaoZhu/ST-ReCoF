#!/bin/Rscript
`%!in%` = Negate(`%in%`)

args=commandArgs(TRUE)
outDir = args[1] # The output directory for saving the results
species = args[2] # Species used for the experiment, only for human and mouse
iniCellLoc = args[3] # ScRNA-seq data of initial cells, which can be in Seurat (.rds) or count (.txt/.csv) format
iniCellType = args[4] # The cell type of the initial scRNA-seq data
tarCellLoc = args[5] # ScRNA-seq data of target cells, which can be in Seurat (.rds) or count (.txt/.csv) format
tarCellType = args[6] # The cell type of the target scRNA-seq data

library(writexl)
source("01_InputProcess.R")
source("02_SPsPreSelection.R")
source("03_NetworkSPsFilter.R")
source("04_MapSPs2SM.R")

# 0. Create the output directory
if(!dir.exists(outDir)){
  dir.create(outDir, recursive = TRUE)
}

# 1. Load the scRNA-seq data of initial and target cells 
iniCellData = scRNALoad(iniCellLoc, iniCellType)
tarCellData = scRNALoad(tarCellLoc, tarCellType)
# Identify differentially expressed genes
seu_merge = merge(iniCellData, tarCellData)
tarInitMarker = gainDEGs(seu_merge, iniCellType, tarCellType)
# Identify differentially expressed transcription factors
# The default threshold for avg_log2FC is 1
log2FC_t = 1
DETFs_list = gainDETFs(tarInitMarker, species, log2FC_t)
initCellCount = getCountMatrix(iniCellLoc)
# (supplement)
# seu_merge = SeuratProcess(seu_merge)

# 2. pre-select the signal protein
# First, compute a similarity score between the transcription factors and the perturbation reference database
if(species == 'human'){
  LINCS_data_list_ff = readRDS("../data/drugPertDB/LINCSDB_6h_human_list_ff.rds")
}else if(species == 'mouse'){
  LINCS_data_list_ff = readRDS("../data/drugPertDB/LINCSDB_6h_mouse_list_ff.rds")
}else{
  stop("Supported species: human or mouse only.")
}
repurposing_data_pp = readRDS("../data/drugPertDB/repurposing_data_pp.rds")
repurposing_data_df = readRDS("../data/drugPertDB/repurposing_data_df.rds")

SimTermRank_df = simiSorceTFsPertDB(DETFs_list, LINCS_data_list_ff)
# Map the obtained small-molecule perturbagens to the drug target database to identify potential signal proteins
preSelectSPs = map2PreSelectSPs(SimTermRank_df, repurposing_data_df, species)
preSelectSPs.dir = paste0(outDir, "/", tarCellType, "-", iniCellType, "_preSelectSPs.rds")
saveRDS(preSelectSPs, preSelectSPs.dir)
# 3. Network based Signal proteins selection and clustering
# Select signal proteins that interact strongly and specifically with transcription factors
SelectSPsFinal = NetworkSPsFilter(initCellCount, species, DETFs_list, preSelectSPs)
SelectSPsFinal.dir = paste0(outDir, "/", tarCellType, "-", iniCellType, "_SelectSPsFinal.xlsx")
writexl::write_xlsx(SelectSPsFinal, SelectSPsFinal.dir)
communities = communitiesDivide(initCellCount, species)

# 4. Match these communities to the corresponding small-molecule perturbagens
if(species == 'human'){
  chemTarDB = readRDS("../data/drugPertDB/chemTarDB_human.rds")
}else if(species == 'mouse'){
  chemTarDB = readRDS("../data/drugPertDB/chemTarDB_mouse.rds")
}else{
  stop("Supported species: human or mouse only.")
}

smallMoleculesCategory = readRDS("../data/publicInfo/smallMoleculesCategory.rds")
pathway_info = readRDS("../data/publicInfo/pathwayInfo.rds")
smallMoluculeFunc = readRDS("../data/publicInfo/smallMoluculeFunction.rds")
# Filter out the unnecessary signal proteins
communities_f = communities[names(communities) %in% SelectSPsFinal$Protein]
ChemScoreDf = pertChemScore(communities_f, chemTarDB, SelectSPsFinal, smallMoleculesCategory)
classfy_df = CommuMajorCategory(pathway_info, ChemScoreDf)
SmallMoleculeResult = gainFinalSmallMolecule(classfy_df, ChemScoreDf, pathway_info, 10)
reprogram_drug = drug_function[drug_function$Reprogramming == 1, 'Chemicals']
SmallMoleculeResult_f = SmallMoleculeResult[SmallMoleculeResult$pert %in% reprogram_drug, ]
SmallMoleculeResult_f = SmallMoleculeResult_f[, colnames(SmallMoleculeResult_f) %!in% c('UnexpectEffect')]
writexl::write_xlsx(classfy_df, paste0(outDir, "/", tarCellType, "-", iniCellType, "_classfy_df.xlsx"))
writexl::write_xlsx(SmallMoleculeResult, paste0(outDir, "/", tarCellType, "-", iniCellType, "_SmallMoleculeResult.xlsx"))
writexl::write_xlsx(SmallMoleculeResult_f, paste0(outDir, "/", tarCellType, "-", iniCellType, "SmallMoleculeResult_f.xlsx"))

