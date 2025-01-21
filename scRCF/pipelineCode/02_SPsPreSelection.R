# 01_SPsPreSelection.R
# The first step is to pre-select the signal protein: 
# 1. First, compute a similarity score between the transcription factors and the perturbation reference database
# 2. Then, filter the obtained small-molecule perturbagens based on their similarity scores using z-score
# 3. Map the obtained small-molecule perturbagens to the drug target database to identify potential signal proteins

library(tidyverse)
library(tidydr)
`%!in%` = Negate(`%in%`)
syn = readRDS("../data/publicInfo/mouse2human.rds")

### Compute the Jaccard similarity score between the queried transcription factors and the reference database
## @DETFs: Queried transcription factors
## @referenceDB: Pre- and post-perturbation small-molecule database
MJaccardScore <- function(DETFs, referenceDB){
  UDETFs<-union(names(DETFs), names(referenceDB))
  mat<-matrix(0,nrow = 2,ncol = length(UDETFs))
  colnames(mat) <- UDETFs
  mat[1,match(names(DETFs),colnames(mat))] <- DETFs
  mat[2,match(names(referenceDB),colnames(mat))] <-referenceDB
  MJaccard <- length(which(mat[1,]*mat[2,]==1)) / length(UDETFs)
  return(MJaccard)
}

### 2. compute a similarity score between the transcription factors and the perturbation reference database
## @DETFs_list: The differentially expressed transcription factors between initial and targe cell type
## @reference_database: Pre- and post-perturbation small-molecule database
simiSorceTFsPertDB <- function(DETFs_list, reference_database, zscore_t = 1){
  SimSorce <- sapply(1:length(reference_database),function(i)MJaccardScore(DETFs_list, reference_database[[i]]))
  pert <- sapply(names(reference_database), function(x) strsplit(x, split = '_')[[1]][1])
  SimRank_df <- data.frame("score"=SimSorce,"pert"=pert,stringsAsFactors = F)
  SimRank_df$score <- as.numeric(SimRank_df$score)
  SimRank_df <- SimRank_df[order(SimRank_df$score,decreasing=TRUE),]
  SimRank_df$rank <- seq(length(SimRank_df$score))
  SimRank_df_f <- SimRank_df[SimRank_df$score > 0, ]
  # 2. Then, filter the obtained small-molecule perturbagens based on their similarity scores using z-score
  scoreThreshold <- mean(SimRank_df_f$score) + zscore_t * sd(SimRank_df_f$score)
  SimRank_df_ff = SimRank_df_f[SimRank_df_f$score > scoreThreshold, ]
  return(SimRank_df_ff)
}

SPsEffect <- function(x){
  effect = NA
  x = unique(x)
  if(length(x) == 3){
    effect = 2
  }else if(length(x) == 2){
    tmp = prod(x)
    effect = ifelse(tmp == -1, 2, ifelse(tmp == 2, 1, -1))
  }else if(length(x) == 1){
    effect = x
  }else{
    effect = NA
  }
  return(effect)
}
convert2mouse <- function(df){
  df$mouseGene = syn[match(df$target, syn$HGNC.symbol), "MGI.symbol"]
  df = df[, c("mouseGene", "effect")]
  colnames(df) = c("target", "effect")
  df = df[!is.na(df$target), ]
  return(df)
}
convert2vector <- function(df){
  tmpVector = df$effect
  names(tmpVector) = df$target
  return(tmpVector)
}

getTermTarget<-function(term, repurposing_data_df){
  term <- tolower(term)
  repurposing_data_with_gene = rownames_to_column(repurposing_data_df, 'gene')
  PertTarget = repurposing_data_with_gene[, c('gene', term)]
  PertTarget <- PertTarget[PertTarget[,2]!=0,]
  colnames(PertTarget) <- c("SignalPro","action")
  return(PertTarget)
}

unique_terms <- function(data) {
  unique_df <- data.frame()
  for (term in unique(data$term)) {
    sub_data <- data[data$term == term, ]
    if (nrow(sub_data) == 1) {
      unique_df <- rbind(unique_df, sub_data)
    } else {
      actions <- sub_data$action
      if (all(actions == actions[1])) {
        unique_df <- rbind(unique_df, sub_data[1, ])
      } else {
        if (-1 %in% actions & 1 %in% actions) {
          unique_action <- 2
        } else if (-1 %in% actions & 2 %in% actions) {
          unique_action <- -1
        } else if (1 %in% actions & 2 %in% actions) {
          unique_action <- 1
        } else {
          unique_action <- actions[1] # 默认情况，保留第一个 action
        }
        unique_df <- rbind(unique_df, data.frame(Freq = sub_data$Freq[1], term = term, action = unique_action))
      }
    }
  }
  return(unique_df)
}

### 3. Map the obtained small-molecule perturbagens to the drug target database to identify potential signal proteins
## @SimRank_df_ff: Entries with high similarity scores in the reference dataset
## @Chem2Tar_df: Database of small-molecule perturbagens and their corresponding targets
## @species: species
map2PreSelectSPs <- function(SimRank_df_ff, repurposing_data_df, species){
  SimRank_target = lapply(SimRank_df_ff$pert, function(pert){
    getTermTarget(pert, repurposing_data_df)
  })
  SimRank_target_p <- do.call(rbind, SimRank_target)
  TermProtein <- paste(SimRank_target_p$SignalPro, SimRank_target_p$action, sep="_")
  TermProtein <- table(TermProtein)
  TermProtein <- as.data.frame(TermProtein,stringsAsFactors = F)
  TermProtein$term <- sapply(TermProtein$TermProtein,function(x)strsplit(x,split = "_")[[1]][1])
  TermProtein$action <- sapply(TermProtein$TermProtein,function(x)strsplit(x,split = "_")[[1]][2])
  TermProtein_f1 <- TermProtein[TermProtein$action!=2,]
  TermProtein_f2 <- TermProtein[TermProtein$action==2,]
  TermProtein_final <- join(TermProtein_f1, TermProtein_f2, by=c("term"),type="full",match="all")
  TermProtein_final<-TermProtein_final[order(TermProtein_final$Freq,decreasing=T),]
  TermProtein_final = TermProtein_final[, c("Freq", "term", "action")]
  TermProtein_final = unique_terms(TermProtein_final)
  # TermProtein_final <- TermProtein_final[!duplicated(TermProtein_final$term),]
  TermProtein_final = TermProtein_final[, c("term", "action")]
  colnames(TermProtein_final) = c("target", "effect")
  if(species == 'mouse'){
    TermProtein_final = convert2mouse(TermProtein_final)
  }
  TermProtein_final = convert2vector(TermProtein_final)
  return(TermProtein_final)
}
