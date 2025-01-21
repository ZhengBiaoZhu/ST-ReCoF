# 04_MapSPs2SM.R
# Map the signaling proteins to small-molecule drugs. 
# First, check the number of communities affected by each small-molecule drug. 
# Then, combine this with the category information of the small molecules to obtain the number of communities affected by the major pathways. 
# Next, use a modified Sorensen-Dice index calculation to determine the similarity scores of small-molecule drugs within each major pathway. 
# Select the small-molecule drugs with the highest scores in each sub-pathway as the final result.

### Obtain the similarity score for each small-molecule drug
## @pert: Small molecules
## @communities_f: Community information after filtering out the unnecessary signal proteins
## @SelectSPsFinal: Selected signaling proteins
## @chemTarDB: Small-molecule target database used
calPertChemScore <- function(pert, communities_f, SelectSPsFinal, chemTarDF){
  pert_info = chemTarDF[chemTarDF$pert_iname == pert, ]
  pert_effect = nrow(pert_info)
  # Key point: Obtain the number of affected communities
  # 1. First, identify the intersection of drugs and communities to determine the affected communities
  # 2. In each community, proteins with the same effect score 1; if there are 2, the score is 0.5; if all are 2, the score is 0.25
  # Sum the scores for each community and divide by the number of intersecting proteins
  target_use <- pert_info$target
  # Obtain the corresponding communities
  community_u <- communities_f[names(communities_f) %in% target_use]
  if(length(community_u) == 0)
    return(0)
  JSD_u <- SelectSPsFinal[SelectSPsFinal$Protein %in% target_use, ]
  cal_df <- data.frame(target = names(community_u), 
                       community = community_u, 
                       pert_effect = pert_info[match(names(community_u), pert_info$target), 'effect'], 
                       ehrich_effect = JSD_u[match(names(community_u), JSD_u$Protein), 'Sign'])
  
  # The highest score is used to represent the score of the community
  pert_final = data.frame(pert = pert, pertTarNum = pert_effect, pertTar = paste0(pert_info$target, collapse = ','))
  return(pert_final)
}

### Obtain the communities affected by each small-molecule drug
## @communities_f: Community information after filtering out the unnecessary signal proteins
## @chemTarDB: Small-molecule target database used
## @SelectSPsFinal: Selected signaling proteins
## @smallMoleculesCategory: Small-molecule Category information
pertChemScore <- function(communities_f, chemTarDB, SelectSPsFinal, smallMoleculesCategory){
  calChemScore = lapply(unique(chemTarDB$pert_iname), function(x) calPertChemScore(x, communities_f, SelectSPsFinal, chemTarDB))
  calChemScore = Filter(function(x) length(x) == 3, calChemScore)
  calChemScore_df = do.call(rbind, calChemScore)
  smallMoleculesCategory$drug_name = tolower(smallMoleculesCategory$drug_name)
  calChemScore_df_f = calChemScore_df[calChemScore_df$pert %in% smallMoleculesCategory$drug_name, ]
  calChemScore_df_f$category1 = smallMoleculesCategory[match(calChemScore_df_f$pert, smallMoleculesCategory$drug_name), 'category1']
  calChemScore_df_f$category2 = smallMoleculesCategory[match(calChemScore_df_f$pert, smallMoleculesCategory$drug_name), 'category2']
  ChemScoreDf = calChemScore_df_f %>%
    separate_rows(category1, sep = ', ')
  return(ChemScoreDf)
}

###
### Obtain the number of communities affected by the major pathways
## @pathway_info: Information containing major pathways and sub-pathways
## @ChemScoreDf: A data.frame showing the communities affected by each small-molecule drug
CommuMajorCategory <- function(pathway_info, ChemScoreDf){
  classfy_df <- data.frame(classfy_main = pathway_info$pathway, 
                           community = NA, 
                           numCommunity = NA)
  for(i in 1:nrow(classfy_df)){
    classfy = classfy_df[i, 'classfy_main']
    tmp_df = ChemScoreDf[ChemScoreDf$category1 == classfy, 'pertTar', drop = FALSE]
    tmp_data = tmp_df %>% 
      separate_rows(pertTar, sep = ',')
    tmp_data$community = communities_f[match(tmp_data$pertTar, names(communities_f))]
    community_num = tmp_data$community[!is.na(tmp_data$community)]
    community_num = community_num[order(community_num, decreasing = F)]
    community_num = unique(community_num)
    classfy_df[i, 'community'] = paste0(community_num, collapse = ',')
    classfy_df[i, 'numCommunity'] = length(community_num)
  }
  return(classfy_df)
}



###
### Calculate the similarity scores between drugs in the major drug categories and the selected signaling proteins
## @pert: Small molecules
## @communities_f: Community information after filtering out the unnecessary signal proteins
## @SelectSPsFinal: Selected signaling proteins
## @pert_category: Corresponding major pathway information
## @classfy_df: The number of communities affected by the major pathways
## @chemTarDB: Small-molecule target database used
communityEffectScore <- function(pert, communities_f, SelectSPsFinal, pert_category, classfy_df, chemTarDB){
  pert_info = data.frame(chemTarDB[chemTarDB$pert_iname == pert, ])
  pert_info$target = pert_info$target
  pert_info$effectCom = communities_f[match(pert_info$target, names(communities_f))]
  pert_info$Sign = SelectSPsFinal[match(pert_info$target, SelectSPsFinal$Protein), 'Sign']
  unexpectEffect = length(pert_info$effectCom[is.na(pert_info$effectCom)])
  if(nrow(pert_info) == unexpectEffect)
    return(0)
  pert_info$product <- as.numeric(pert_info$effect) * as.numeric(pert_info$Sign)
  pert_info$prodEffect <- ifelse(pert_info$product == 1, 1, ifelse(abs(pert_info$product) == 2, 0.5, ifelse(pert_info$product == 4, 0.25, 0)))
  pert_info_final <- pert_info %>%
    group_by(effectCom) %>% 
    top_n(1, wt = prodEffect) %>% 
    # distinct(prodEffect, .keep_all = TRUE) %>%
    ungroup()
  orign_Effect = sum(pert_info_final$prodEffect)
  pert_info_final <- pert_info %>%
    group_by(effectCom) %>% 
    top_n(1, wt = prodEffect) %>% 
    distinct(prodEffect, .keep_all = TRUE) %>%
    ungroup()
  Effect = sum(pert_info_final$prodEffect)
  category_Eff = classfy_df[classfy_df$classfy_main == pert_category, 'numCommunity']
  # MJaccardScore = (Effect / category_Eff) * (Effect / (unexpectEffect + Effect))
  MJaccardScore = (Effect * 2) / (Effect + category_Eff) + (1 - (orign_Effect - Effect) / Effect) * 1e-5
  return(MJaccardScore)
}

### Obtain the corresponding small-molecule results
## @classfy_df: The number of communities affected by the major pathways
## @chemTarDB: Small-molecule target database used
## @communities_f: Community information after filtering out the unnecessary signal proteins
## @pathway_info: Information containing major pathways and sub-pathways
## @chemTarDB: Small-molecule target database used
gainFinalSmallMolecule <- function(classfy_df, ChemScoreDf, communities_f, pathway_info, chemTarDB, SelectSPsFinal){
  SelectSPsFinal = data.frame(SelectSPsFinal)
  classfy_main = classfy_df$classfy_main
  result_list = list()
  chemTarDB = chemTarDB %>%
    group_by(pert_iname) %>%
    filter(n() <= 10) %>%
    ungroup()
  for(i in 1:length(classfy_main)){
    print(i)
    classfy = classfy_main[i]
    pert_classfy = ChemScoreDf[ChemScoreDf$category1 == classfy, 'pert']
    pert_classfy_df = ChemScoreDf[ChemScoreDf$category1 == classfy, ]
    # Obtain all drugs in the category
    # a) Drug target information is obtained from repurposing_data_pp
    # b) Protein interaction direction is obtained from JSD_re_final
    # c) Protein grouping is obtained from communities_f
    tmp_result_df = data.frame(pert = pert_classfy$pert, MJScore = NA, category1 = NA, category2 = NA)
    for(j in 1:nrow(pert_classfy)){
      pert = tmp_result_df[j, 'pert']
      tmp_result_df[j, 'MJScore'] = communityEffectScore(pert, communities_f, SelectSPsFinal, classfy, classfy_df, chemTarDB)
    }
    tmp_result_df$category1 = classfy
    tmp_result_df$category2 = pert_classfy_df[match(tmp_result_df$pert, pert_classfy_df$pert), 'category2', drop = T]
    
    pathway_info_f  = pathway_info[pathway_info$pathway == classfy, 'sub_pathway']
    sub_pathway = strsplit(pathway_info_f, split = ', ')[[1]]
    tmp_result_df_f = tmp_result_df %>% 
      separate_rows(category2, sep = ', ')
    tmp_result_df_f = tmp_result_df_f[tmp_result_df_f$category2 %in% sub_pathway, ]
    tmp_result_df_f = tmp_result_df_f[tmp_result_df_f$MJScore != 0, ]
    
    result <- tmp_result_df_f %>%
      group_by(category2) %>%
      top_n(1, MJScore) %>%
      ungroup()
    result_list[[i]] = result
    names(result_list)[i] = classfy
  }
  result_df = do.call(rbind, result_list)
  return(result_df)
}

