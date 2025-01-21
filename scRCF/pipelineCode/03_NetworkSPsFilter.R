# 03_networkSPsFilter
# This part is referenced from:  
# Zheng M, Xie B, Okawa S, Liew SY, Deng H, del Sol A. A single cell-based computational platform to identify chemical compounds targeting desired sets of transcription factors for cellular conversion. Stem Cell Reports. 2023;18(1):131-144. doi:10.1016/j.stemcr.2022.10.013
# In this part, network filtering is performed. 
# First, construct a gene co-expression network using the scRNA count matrix to obtain the weights between the corresponding nodes in the network. 
# Add the obtained differentially expressed transcription factors and pre-selected signaling proteins. 
# Next, calculate the accessibility and specificity scores from the pre-selected signaling proteins to the differentially expressed transcription factors and integrate the results. 
# Then, use JSD (Jensen-Shannon Divergence) to measure the efficiency of each pre-selected signaling protein.
# After obtaining the selected proteins, incorporate the network information and use the random walk algorithm to perform community detection on them

library(igraph)
library(plyr)
library(foreach)
library(doParallel)

numWorkers =4
registerDoParallel(numWorkers)

## calculate Entropy
NodeEntropy<-function(nodeVector){
  idx<-which(nodeVector>0)
  if(length(idx)>1){
    En<- -sum(nodeVector[idx]*log(nodeVector[idx]))
    NEn<- En/log(length(idx)) ## normalized by the maximum entropy
  }else{
    NEn<-0
  }
  return(NEn)
}

## assign entropy and edge weight to the background network based on scRNA-seq of initial cell state
## @input_data: the scRNA-seq data of initial cell state, tEach row means geneand each column means cell
## @species: sepcies can be human, mouse or rat
BN_assign<-function(input_data,species){
  #  input_data<-readRDS(input_data)
  if(species=='mouse'|species=='rat'){
    network<-readRDS('../data/PKN/mouseNrat/mouse_merged_graph.rds')
    iTFs2reTFs<-readRDS('../data/PKN/mouseNrat/iTFs2reTFs_mouse.rds')
  }else if(species=="human"){
    network<-readRDS('../data/PKN/human/human_merged_graph.rds')
    iTFs2reTFs<-readRDS('../data/PKN/human/iTFs2reTFs_human.rds')
  }else {
    stop("Only human, mouse, rat are support!")
  }
  
  ## assgin expression value to sources of interactions
  b=input_data
  b$from<-rownames(b)
  a<-get.data.frame(network,"edges")
  ab=join(a,b,by=c("from"),type="left",match="first")
  
  ## assgin expression value to targets of interactions
  colnames(b)[ncol(b)]<-"to"
  ab1<-join(a,b,by=c("to"),type="left",match="first")
  ab<-ab[,4:ncol(ab),drop=F]
  ab1<-ab1[,4:ncol(ab1),drop=F]
  
  ## elementwise product
  elePro<-ab * ab1
  
  ## Sum of elementwise product
  sum_product_expression=rowSums(elePro, na.rm = FALSE, dims = 1)
  
  ## add sum_product_expression to edge data frame  
  g3=cbind(a,sum_product_expression)
  g3[is.na(g3)]<-0
  ## discard the edges that have zero production expression
  g3<-g3[g3$sum_product_expression!=0,]
  
  ## add transcriptional level network
  iTFs2reTFs_f<-iTFs2reTFs[iTFs2reTFs$from%in%g3$to,]
  colnames(b)[ncol(b)]<-"from"
  ab3=join(iTFs2reTFs_f,b,by=c("from"),type="left",match="first")
  iTFs2reTFs_f$sum_product_expression<-rowSums(ab3[,4:ncol(ab3)])
  
  ## combine two levels network
  g3<-rbind(g3,iTFs2reTFs_f)
  
  ## calculate the out degree strength for each nodes, named as outstrength
  g <- graph.data.frame(as.data.frame(g3))
  E(g)$weight<-g3$sum_product_expression
  V(g)$strength<-strength(g,vids = V(g),mode = "out")
  g3$outstrength<-V(g)$strength[match(g3$from,as_ids(V(g)))]
  
  ## assign probability to each interactions(edges)
  g3$prob<-g3$sum_product_expression/g3$outstrength
  E(g)$prob<-g3$prob
  ## calculate the entropy for each node
  g_adj<-as_adjacency_matrix(g,attr ="prob")
  V(g)$entropy<-sapply(1:nrow(g_adj),function(x)NodeEntropy(g_adj[x,]))
  return(g)
}

## if the query TFs do not locate in the first layer, identify the first layer TF as proxy TF that can target non-first layer TFs as much as possible, while minimizing off target effect.
## contextualize the prior knowledge TF interaction network that is specific to the initial cellular state
## @hairball_g: the species specific (human/moust&rat) prior knowledge TF interaction network
## @initial_scRNA: the scRNA-seq data of initial cell state, tEach row means geneand each column means cell
network_context<-function(hairball_g,initial_scRNA){
  
  hairball_df<-get.data.frame(hairball_g,"edges")
  
  dat_percent<-rowSums(initial_scRNA!=0)/ncol(initial_scRNA)
  bool<-as.data.frame(ifelse(dat_percent<0.5,-1,1))
  bool$to<-rownames(bool)
  colnames(bool)=c("to_effect","to")
  join_obj<-join(hairball_df,bool,by="to",type="left",match="first")
  colnames(bool)=c("from_effect","from")
  join_obj<-join(join_obj,bool,by="from",type="left",match="first")
  join_obj$to_effect[which(is.na(join_obj$to_effect))]=-1
  join_obj$from_effect[which(is.na(join_obj$from_effect))]=-1
  join_obj$pred<-join_obj[,4]*join_obj[,5]
  join_obj_unspe<-join_obj[join_obj$Effect==2,]
  join_obj_spe<-join_obj[join_obj$Effect!=2,]
  join_obj_spe_contex<-join_obj_spe[join_obj_spe$Effect==join_obj_spe$pred,]
  join_obj_contex<-rbind(join_obj_spe_contex,join_obj_unspe)
  g <- graph.data.frame(join_obj_contex[,1:3])
  return(g)
}

## Identify the 1stTF as proxy of query non-1stTFs
## @context_g: the contextualized prior knowledge TF interaction network from function network_context
## @firstTFs: potential 1stTFs that could be proxy of query non-1stTFs
## @nonFirstTFs: query non-1stTFs, the query TFs that have no direct interactions with interface TFs
proxyTFs_iden<-function(context_g,firstTFs,nonFirstTFs){
  
  firstTFs<-firstTFs[firstTFs%in%as_ids(V(context_g))]
  # identify the length from first layer TFs to desired non-first layer TFs
  Non1stTFsProxy_cal <-foreach(perttar=firstTFs, .combine=rbind, .inorder=FALSE, .multicombine=TRUE, .maxcombine = length(firstTFs),.packages="igraph") %dopar% {
    x <- distances(context_g,perttar,nonFirstTFs,mode="out",weights=NA)
  }
  Non1stTFsProxy_cal[is.infinite(Non1stTFsProxy_cal)]<-0
  ReachNum<-rowSums(Non1stTFsProxy_cal!=0)
  Non1stTFsProxy_cal_reach<-data.frame(Non1stTFsProxy_cal[which(ReachNum!=0),],stringsAsFactors = F)
  # check how many nodes are effected with length as desired non-first layer TFs been effected
  allReachNode<-lapply(1:nrow(Non1stTFsProxy_cal_reach), function(x)as_ids(ego(context_g, order = max(Non1stTFsProxy_cal_reach[x,]), nodes = rownames(Non1stTFsProxy_cal_reach)[x], mode = "out")[[1]]))
  allReachNode_matrix<-lapply(allReachNode, function(x){VectorPred<-rep(0,length(as_ids(V(context_g))))
  names(VectorPred)=as_ids(V(context_g))
  VectorPred[x]=1/length(x)
  return(VectorPred)
  })
  allReachNode_matrix<-do.call(rbind,allReachNode_matrix)
  rownames(allReachNode_matrix)<-rownames(Non1stTFsProxy_cal_reach)
  VectorIdeal<-rep(0,length(as_ids(V(context_g))))
  names(VectorIdeal)=as_ids(V(context_g))
  VectorIdeal[nonFirstTFs]=1/length(nonFirstTFs)
  JSD_val<-sapply(rownames(allReachNode_matrix),function(x){
    calScore<-allReachNode_matrix[x,]
    if(sum(calScore)==0){return(1)}
    M<-0.5*(calScore+VectorIdeal)
    DKL1<-sum(VectorIdeal*log2(VectorIdeal/M),na.rm = T)
    DKL2<-sum(calScore*log2(calScore/M),na.rm = T)
    0.5*(DKL1+DKL2)
  })
  return(names(JSD_val[which(JSD_val==min(JSD_val))])[1]) # if multiple proxy with same specificity, randomly select one.
}
## calculate the JSD value for signaling proteins based on shortest path and entropy
## @DETFs: DETF vector with names and value with 1 and -1 which mean activation and inhibition of the TF
## @g: background network igraph calculated by function BN_assign
## @out_dir: output directory
## @candidate: signaling protein candidates selected from NCPC
## @out_fiel: the output file name
## @initial_scRNA: the scRNA-seq data of initial cell state, tEach row means geneand each column means cell
PathEntropy<-function(species,DETFs,g,candidate,initial_scRNA){
  if(species=='mouse'){
    network_1stTFs<-read.table('../data/PKN/mouseNrat/mouse_1stTFs.txt',stringsAsFactors = F,sep = '\t')
    sigmole<-read.table('../data/PKN/mouseNrat/mouse_SigMole.txt',stringsAsFactors = F,sep = '\t')
    hairball_g<-readRDS('../data/PKN/mouseNrat/mouse_TransReg_TFHairball.rds')
    te_pert_DETF<-DETFs
    
  }else if(species=='human'){
    network_1stTFs<-read.table('../data/PKN/human/human_1stTFs.txt',stringsAsFactors = F,sep = '\t')
    sigmole<-read.table('../data/PKN/human/human_SigMole.txt',stringsAsFactors = F,sep = '\t')
    hairball_g<-readRDS('../data/PKN/human/human_TransReg_TFHairball.rds')
    te_pert_DETF<-DETFs
  }else {
    stop("Only human, mouse, rat are support!")
  }
  
  E(g)$logW<-log2(E(g)$sum_product_expression+1)
  E(g)$surprisal<-log(1/(E(g)$logW/max(E(g)$logW))) 
  #te_pert_target<-te_pert_target[te_pert_target%in%as_ids(V(g))]
  te_pert_DETF_f<-names(te_pert_DETF)[names(te_pert_DETF)%in%network_1stTFs$V1&names(te_pert_DETF)%in%as_ids(V(g))]
  pert_SigMole<-names(candidate)[names(candidate)%in%as_ids(V(g))&names(candidate)%in%sigmole$V1]
  firstTFs<-network_1stTFs$V1[network_1stTFs$V1%in%as_ids(V(g))]
  ## if the query TFs do not locate in the first layer, identify proxy TF in the first layer targeting on these non-first layer TFs
  network_cont<-network_context(hairball_g,initial_scRNA)
  network_cont_df<-get.data.frame(network_cont,"edges")
  te_pert_DETF_non1st<-names(te_pert_DETF)[!names(te_pert_DETF)%in%network_1stTFs$V1&names(te_pert_DETF)%in%network_cont_df$to]
  if(length(te_pert_DETF_non1st)!=0){
    proxyTF<-proxyTFs_iden(network_cont,firstTFs,te_pert_DETF_non1st)
    te_pert_DETF_f<-unique(c(te_pert_DETF_f,proxyTF))
    te_pert_DETF_f<-te_pert_DETF_f[te_pert_DETF_f%in%as_ids(V(g))] # to ensure the query TFs being targeted by expressed interface TFs
  }
  SigPathScore2 <-foreach(perttar=pert_SigMole, .combine=rbind, .inorder=FALSE, .multicombine=TRUE, .maxcombine = length(pert_SigMole),.packages="igraph") %dopar% {
    # find "shortest" path from a signalling molecule to first layer TFs
    x <- shortest_paths(g,perttar,firstTFs,mode="out",weights=E(g)$surprisal,output = "epath")
    edges = x$epath
    E=rep(NA,length(edges))
    P=rep(NA,length(edges))
    if(length(which(lengths(edges)!=0))!=0){
      ##calculate entropy,e.g.specifity
      E[which(lengths(edges)!=0)]<- sapply(which(lengths(edges)!=0), function(y)sum(V(g)$entropy[tail_of(g,edges[[y]])]))
      E_re<-1/exp(as.numeric(E))
      ## calculate probablity,e.g.reachability
      P[which(lengths(edges)!=0)]<- sapply(which(lengths(edges)!=0), function(y)sum(E(g)$surprisal[edges[[y]]]))
      P_re<-1/exp(as.numeric(abs(P)))
      Score<-E_re*P_re
    }else{
      Score<-P
    }
    return(Score)
  }
  
  SigPathScore2_s<-SigPathScore2
  rownames(SigPathScore2_s)<-pert_SigMole
  colnames(SigPathScore2_s)<- firstTFs
  
  SigPathScore2_s[is.na(SigPathScore2_s)]<-0
  SigPathScore2_s<-as.data.frame(SigPathScore2_s)
  idealScore<-rep(0,ncol(SigPathScore2_s))
  names(idealScore)<-colnames(SigPathScore2_s)
  idealScore[te_pert_DETF_f]<-1/length(te_pert_DETF_f)
  
  JSD_val<-sapply(rownames(SigPathScore2_s),function(SigMole){
    calScore<-SigPathScore2_s[SigMole,]
    if(sum(calScore)==0){return(1)}
    calScore<-calScore/sum(calScore)
    M<-0.5*(calScore+idealScore)
    DKL1<-sum(idealScore*log2(idealScore/M),na.rm = T)
    DKL2<-sum(calScore*log2(calScore/M),na.rm = T)
    0.5*(DKL1+DKL2)
  })
  
  JSD_enrich_name<-names(JSD_val)[order(JSD_val,decreasing=F)]
  JSD_val_order<-JSD_val[order(JSD_val,decreasing=F)]
  sign<-candidate[match(JSD_enrich_name,names(candidate))]
  JSD_enrich_rank<-data.frame(cbind("Protein"=JSD_enrich_name,"Rank"=seq(length(JSD_enrich_name)),"Sign"=sign,"JSD_val"=JSD_val_order))
  return(JSD_enrich_rank)
}

### Further filter the protein targets using the PPI network
## @initCellCount: Count matrix of the initial cell types
## @species: species
## @DETFs_list: List of differentially expressed transcription factors
## @preSelectSPs: Preliminarily selected signaling proteins
NetworkSPsFilter <- function(initCellCount, species, DETFs_list, preSelectSPs){
  MNetwork <- BN_assign(initCellCount, species)
  JSDSelectSPs <- PathEntropy(species, DETFs_list, MNetwork, preSelectSPs, initCellCount)
  SelectSPsFinal <- JSDSelectSPs[which(JSDSelectSPs$JSD_val != 1), ]
  return(SelectSPsFinal)
}

### Divide the filtered signaling proteins into different communities using the random walk algorithm
## @initCellCount: Count matrix of the initial cell types
## @species: species
communitiesDivide <- function(initCellCount, species){
  MNetwork <- BN_assign(initCellCount, species)
  clusters <- cluster_walktrap(MNetwork, 
                               weights = E(MNetwork)$weight,
                               merges = T)
  membership <- membership(clusters)
  return(membership)
}
