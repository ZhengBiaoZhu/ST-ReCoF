# 04_MapSPs2SM.R
# Gathering a list of candidate small molecules 
# Based on the community partitioning results, first calculate the effect score of each small molecule in the small-molecule target database.
# Next, calculate the impact scores for the signaling pathway classifications and retain the top 10 scoring classifications.
# Finally, retain the small molecule with the highest effect score within the sub-pathways of each signaling pathway classification.

# Calculate the effect score of each small molecule in the small-molecule target database
## @pert: Small molecules in the small-molecule target database
## @communities_f: Community division results
## @SelectSPsFinal: Final list of signaling proteins
## @chemTarDB: The small-molecule target database
calPertChemScore <- function(pert, communities_f, SelectSPsFinal, chemTarDB){
  SelectSPsFinal = data.frame(SelectSPsFinal)
  pert_info = chemTarDB[chemTarDB$pert_iname == pert, ]
  pert_effect = nrow(pert_info)
  target_use <- pert_info$target
  community_u <- communities_f[names(communities_f) %in% target_use]
  if(length(community_u) == 0)
    return(0)
  pert_info$community = communities_f[match(pert_info$target, names(communities_f))]
  UnexpectEffect = sum(is.na(pert_info$community))
  JSD_u <- SelectSPsFinal[SelectSPsFinal$Protein %in% target_use, ]
  cal_df <- data.frame(target = names(community_u), 
                       community = community_u, 
                       pert_effect = pert_info[match(names(community_u), pert_info$target), 'effect'], 
                       ehrich_effect = JSD_u[match(names(community_u), JSD_u$Protein), 'Sign'])
  colnames(cal_df) = c('target', 'community', 'pert_effect', 'enrich_effect')
  cal_df$product = NA
  cal_df$product = as.numeric(cal_df$pert_effect) * as.numeric(cal_df$enrich_effect)
  cal_df$effect = ifelse(cal_df$product == 1, 1, ifelse(abs(cal_df$product) == 2, 0.5, ifelse(cal_df$product == 4, 0.25, 0)))
  cal_df = cal_df %>% 
    group_by(community) %>%
    slice_max(order_by = product, n = 1, with_ties = FALSE) %>%
    ungroup()
  pert_final = data.frame(pert = pert, pertTarNum = pert_effect, pertTar = paste0(pert_info$target, collapse = ','), 
                          EffectScore = sum(cal_df$effect), Community = paste0(cal_df$community, collapse = ','), 
                          UnexpectEffect = UnexpectEffect)
  return(pert_final)
}

# Add classification information to the small molecules and filter out those with an EffectScore of 0
## @communities_f: Community division results
## @chemTarDB: The small-molecule target database
## @SelectSPsFinal: Final list of signaling proteins
## @smallMoleculesCategory: Small molecule classification information
pertChemScore <- function(communities_f, chemTarDB, SelectSPsFinal, smallMoleculesCategory){
  calChemScore <- list()
  unique_pert_inames <- unique(chemTarDB$pert_iname)
  for (i in 1:length(unique_pert_inames)) {
    x <- unique_pert_inames[i]
    result <- calPertChemScore(x, communities_f, SelectSPsFinal, chemTarDB)
    if (length(result) == 6) {
      calChemScore[[length(calChemScore) + 1]] <- result
    }
  }
  calChemScore_df <- do.call(rbind, calChemScore)
  
  smallMoleculesCategory$drug_name = tolower(smallMoleculesCategory$drug_name)
  calChemScore_df_f = calChemScore_df[calChemScore_df$pert %in% smallMoleculesCategory$drug_name, ]
  calChemScore_df_f$category1 = smallMoleculesCategory[match(calChemScore_df_f$pert, smallMoleculesCategory$drug_name), 'category1']
  calChemScore_df_f$category2 = smallMoleculesCategory[match(calChemScore_df_f$pert, smallMoleculesCategory$drug_name), 'category2']
  ChemScoreDf = calChemScore_df_f %>%
    separate_rows(category1, sep = ', ')
  ChemScoreDf = ChemScoreDf[ChemScoreDf$EffectScore > 0, ]
  return(ChemScoreDf)
}

# Calculate the impact score for each pathway
## @pathway_info: Signal pathway classification and its sub-pathway
## @ChemScoreDf: Small molecules along with their effect scores and classification information
CommuMajorCategory <- function(pathway_info, ChemScoreDf){
  classfy_df <- data.frame(classfy_main = pathway_info$pathway, 
                           community = NA, 
                           numCommunity = NA, 
                           score = NA)
  for(i in 1:nrow(classfy_df)){
    classfy = classfy_df[i, 'classfy_main']
    tmp_df = ChemScoreDf[ChemScoreDf$category1 == classfy, , drop = FALSE]
    drug_num = nrow(tmp_df)
    tmp_data = data.frame(community = tmp_df$Community) %>% 
      separate_rows(community, sep = ',')
    community_num = tmp_data$community[!is.na(tmp_data$community)]
    community_num = community_num[order(community_num, decreasing = F)]
    community_num = unique(community_num)
    classfy_df[i, 'community'] = paste0(community_num, collapse = ',')
    classfy_df[i, 'numCommunity'] = sum(tmp_df$EffectScore)
    classfy_df[i, 'score'] = classfy_df[i, 'numCommunity'] / drug_num
  }
  return(classfy_df)
}

# Filter small molecule drugs to obtain the highest-scoring small molecule in each sub-pathway
## @classfy_df: Signal pathway impact score results
## @ChemScoreDf: Small molecules along with their effect scores and classification information
## @pathway_info: Signal pathway classification and its sub-pathway
## num: The number of saved signal pathways
gainFinalSmallMolecule <- function(classfy_df, ChemScoreDf, pathway_info, num){
  filter_num = 10
  classfy_df_f = classfy_df %>%
    slice_max(order_by = score, n = num)
  pathway_info_f = pathway_info[pathway_info$pathway %in% classfy_df_f$classfy_main,]
  chemTarDB = chemTarDB %>%
    group_by(pert_iname) %>%
    filter(n() <= filter_num) %>%
    ungroup()
  result_list = list()
  for(i in 1:num){
    pathway = classfy_df_f[i, 'classfy_main']
    tmp_ChemScoreDf = ChemScoreDf[ChemScoreDf$category1 == pathway, ]
    MAX_score = max(tmp_ChemScoreDf$EffectScore)
    tmp_ChemScoreDf = tmp_ChemScoreDf %>%
      separate_rows(category2, sep = ', ')
    # tmp_ChemScoreDf$Score = tmp_ChemScoreDf$EffectScore
    tmp_ChemScoreDf$Score = tmp_ChemScoreDf$EffectScore - (tmp_ChemScoreDf$UnexpectEffect / tmp_ChemScoreDf$pertTarNum * 1e-5)
    tmp_ChemScoreDf$Score = tmp_ChemScoreDf$Score / MAX_score
    tmp_ChemScoreDf_f = tmp_ChemScoreDf %>%
      group_by(category2) %>%
      slice_max(order_by = Score, n = 1, with_ties = T)
    sub_drug = data.frame(pathway_info_f[pathway_info_f$pathway == pathway, 'sub_drug']) %>%
      separate_rows(sub_drug, sep = ', ')
    tmp_ChemScoreDf_f = tmp_ChemScoreDf_f[tmp_ChemScoreDf_f$category2 %in% sub_drug$sub_drug,]
    # tmp_ChemScoreDf_f = tmp_ChemScoreDf_f[tmp_ChemScoreDf_f$UnexpectEffect < 3, ]
    result_list[[i]] = tmp_ChemScoreDf_f
  }
  result_df = do.call(rbind, result_list)
  return(result_df)
}








