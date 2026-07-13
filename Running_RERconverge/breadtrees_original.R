#PURPOSE: alternate function to readTrees(), can create subsets of the masterTree if given list of desired species  

library(RERconverge)
breadTrees1=function(file, max.read=NA, masterTree=NULL, minTreesAll=20, reestimateBranches=F, minSpecs=NULL,useSpecies=NULL){
  
  
  
  message("Using readTrees2")
  message("Reading data")
  tmp=scan(file, sep="\t", what="character", quiet = T)
  message(paste0("Read ",length(tmp)/2, " items", collapse=""))
  trees=vector(mode = "list", length = min(length(tmp)/2,max.read, na.rm = T))
  keeptrees=rep(TRUE,length(trees))
  
  treenames=character()
  maxsp=0; # maximum number of species
  allnames=NA # unique tip labels in gene trees
  
  #create trees object, get species names and max number of species
  for ( i in 1:min(length(tmp),max.read*2, na.rm = T)){
    if (i %% 2==1){
      treenames=c(treenames, tmp[i])
    }
    else{
      trees[[i/2]]=unroot(read.tree(text=tmp[i]))

      # check master tree
      if (!is.null(masterTree)) {
        trees[[i/2]] = pruneTree(trees[[i/2]],intersect(trees[[i/2]]$tip.label,masterTree$tip.label))
      }
      
      #check useSpecies
      if(!is.null(useSpecies)){
        if(length(intersect(trees[[i/2]]$tip.label,useSpecies)) < 3){
          keeptrees[[i/2]] = FALSE; next;
        }
        trees[[i/2]] = unroot(keep.tip(trees[[i/2]],intersect(trees[[i/2]]$tip.label,useSpecies)))
      }
      
      #check for new species
      if (sum(trees[[i/2]]$tip.label %in% allnames == F) > 0) {
        allnames = unique(c(allnames,trees[[i/2]]$tip.label))
        maxsp = length(allnames) - 1
        
      }
    }
    
  }
  allnames = allnames[!is.na(allnames)]
  trees = trees[keeptrees]
  treenames = treenames[keeptrees]
  names(trees)=treenames
  treesObj=vector(mode = "list")
  treesObj$trees=trees
  treesObj$numTrees=length(trees)
  treesObj$maxSp=maxsp
  
  
  message(paste("max is", maxsp))
  
  
  
  ### report is a binary matrix showing the species membership of each tree
  report=matrix(nrow=treesObj$numTrees, ncol=maxsp)
  colnames(report)=allnames
  
  rownames(report)=treenames
  for ( i in 1:nrow(report)){
    ii=match(allnames, trees[[i]]$tip.label)
    report[i,]=1-is.na(ii)
    
  }
  treesObj$report=report
  
  ############ This line finds indices of trees that have the complete species
  ii=which(rowSums(report)==maxsp)
  
  
  
  ######################################################################
  if(length(ii)==0 & is.null(masterTree)){
    stop("no tree has all species - you must supply a master tree")
  }
  ######################################################################
  
  
  #Create a master tree with no edge lengths
  if (is.null(masterTree)) {
    master=trees[[ii[1]]]
    master$edge.length[]=1
    treesObj$masterTree=master
  } else {
    
    master=pruneTree(masterTree, intersect(masterTree$tip.label,allnames))
    #prune tree to just the species names in the largest gene tree
    master$edge.length[]=1
    
    master=unroot(pruneTree(masterTree, intersect(masterTree$tip.label,allnames)))
    #prune tree to just the species names in the gene trees
    #master$edge.length[]=1
    
    treesObj$masterTree=master
  }
  
  
  treesObj$masterTree=rotateConstr(treesObj$masterTree, sort(treesObj$masterTree$tip.label))
  #this gets the abolute alphabetically constrained order when all branches
  #are present
  tiporder=treeTraverse(treesObj$masterTree)
  
  #treesObj$masterTree=CanonicalForm(treesObj$masterTree)
  message("Rotating trees")
  
  for ( i in 1:treesObj$numTrees){
    
    treesObj$trees[[i]]=rotateConstr(treesObj$trees[[i]], tiporder)
    
  }
  
  
  ap=allPathsTrackBranches(master)
  treesObj$ap=ap
  matAnc=(ap$matIndex>0)+1-1
  matAnc[is.na(matAnc)]=0
  
  
  
  paths=matrix(nrow=treesObj$numTrees, ncol=length(ap$dist))
  for( i in 1:treesObj$numTrees){
    #Make paths all NA if tree topology is discordant
    paths[i,]=allPathMasterRelativeTrackBranches(treesObj$trees[[i]], master, ap,i)
    
    #calls matchAllNodes -> matchNodesInject
  }
  paths=paths+min(paths[paths>0], na.rm=T)
  treesObj$paths=paths
  treesObj$matAnc=matAnc
  treesObj$matIndex=ap$matIndex
  treesObj$lengths=unlist(lapply(treesObj$trees, function(x){sqrt(sum(x$edge.length^2))}))
  
  #require all species and tree compatibility
  #ii=which(rowSums(report)==maxsp)
  ii=intersect(which(rowSums(report)==maxsp),which(is.na(paths[,1])==FALSE))
  
  
  
  #if masterTree is provided by user, must use minSpecs<maxsp
  #if no user supplied tree and not minSpec, calculate branch lengths from trees with all species
  #if minSpecs<maxsp, calculate branch lengths from trees with minSpecs species
  if(is.null(minSpecs)){
    #if minimum species not specified,
    #minimum is all species
    minSpecs=maxsp
  }
  
  if(!is.null(masterTree) && !reestimateBranches){
    message("Using user-specified master tree")
  }
  
  
  if(minSpecs==maxsp){ #if we're using all species
    if (is.null(masterTree)) { #and if the user did not specify a master tree
      if(length(ii)>=minTreesAll){
        message (paste0("estimating master tree branch lengths from ", length(ii), " genes"))
        tmp=lapply( treesObj$trees[ii], function(x){x$edge.length})
        
        allEdge=matrix(unlist(tmp), ncol=2*maxsp-3, byrow = T)
        allEdge=scaleMat(allEdge)
        allEdgeM=apply(allEdge,2,mean)
        treesObj$masterTree$edge.length=allEdgeM
      }else {
        message("Not enough genes with all species present: master tree has no edge.lengths")
      }
    }else{
      message("Must specify minSpecs when supplying a master tree: master tree has no edge.lengths")
    }
  }else{ #if we are not using all species
    #estimating from trees with minimum number of species
    treeinds=which(rowSums(report)>=minSpecs) #which trees have the minimum species
    message (paste0("estimating master tree branch lengths from ", length(treeinds), " genes"))
    
    
    if(length(treeinds)>=minTreesAll){
      pathstouse=treesObj$paths[treeinds,] #get paths for those trees
      
      colnames(pathstouse) = ap$destinNode
      colBranch = vector("integer",0)
      unq.colnames = unique(colnames(pathstouse))
      
      for (i in 1:length(unq.colnames)){
        ind.cols = which(colnames(pathstouse) == unq.colnames[i])
        colBranch = c(colBranch,ind.cols[1])
      }
      
      allEdge = pathstouse[,colBranch]
      allEdgeScaled = allEdge
      for (i in 1:nrow(allEdgeScaled)){
        allEdgeScaled[i,] = scaleDistNa(allEdgeScaled[i,])
      }
      colnames(allEdgeScaled) = unq.colnames
      
      edgelengths = vector("double", ncol(allEdgeScaled))
      
      edge.master = treesObj$masterTree$edge
      
      for (i in 1:nrow(edge.master)){
        destinNode.i = edge.master[i,2]
        col.Node.i = allEdgeScaled[,as.character(destinNode.i)]
        edgelengths[i] = mean(na.omit(col.Node.i))
      }
      
      treesObj$masterTree$edge.length = edgelengths
    }else{
      message("Not enough genes with minSpecs species present: master tree has no edge.lengths")
    }
  }
  
  message("Naming columns of paths matrix")
  colnames(treesObj$paths)=namePathsWSpecies(treesObj$masterTree)
  class(treesObj)=append(class(treesObj), "treesObj")
  treesObj
}

treeTraverse=function(tree, node=NULL){
  if(is.null(node)){
    rt=getRoot(tree)
    ic=getChildren(tree,rt)
    return(c(treeTraverse(tree, ic[1]), treeTraverse(tree, ic[2])))
    
  }
  else{
    if (node<=length(tree$tip)){
      return(tree$tip[node])
    }
    else{
      ic=getChildren(tree,node)
      return(c(treeTraverse(tree, ic[1]), treeTraverse(tree, ic[2])))
      
    }
  }
}
getRoot = function(phy) phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]
getChildren=function(tree, nodeN){
  tree$edge[tree$edge[,1]==nodeN,2]
}
allPathsTrackBranches=function(tree){
  dd=dist.nodes(tree) #### pairwise distances between nodes in the tree
  allD=double()
  nn=matrix(nrow=0, ncol=2)
  nA=length(tree$tip.label)+tree$Nnode ######### Total number of nodes in the tree (internal + tips)
  matIndex=matrix(nrow=nA, ncol=nA)
  index=1
  
  destinNode = vector("integer", 0)
  ancNode = vector("integer",0)
  
  for ( i in 1:nA){
    ia=getAncestors(tree,i) #### Getting the ancestors of each node in the tree
    
    destinNode = c(destinNode, rep(i, length(ia)))
    ancNode = c(ancNode, ia)
    
    if(length(ia)>0){
      allD=c(allD, dd[i, ia])
      nn=rbind(nn,cbind(rep(i, length(ia)), ia))
      for (j in ia){
        matIndex[i,j]=index
        index=index+1
      }
    }
  }
  return(list(dist=allD, nodeId=nn, matIndex=matIndex, destinNode=destinNode, ancNode=ancNode))
}
getAncestors=function(tree, nodeN){
  if(is.character(nodeN)){
    nodeN=which(tree$tip.label==nodeN)
  }
  im=which(tree$edge[,2]==nodeN)
  if(length(im)==0){
    return()
  }
  else{
    anc=tree$edge[im,1]
    return(c(anc, getAncestors(tree, anc)))
  }
  
}
allPathMasterRelativeTrackBranches=function(tree, masterTree, masterTreePaths=NULL,i=NULL){
  if(! is.list(masterTreePaths)){
    masterTreePaths=allPathsTrackBranches(masterTree)
  }
  
  treePaths=allPaths(tree)
  map=matchAllNodes(tree,masterTree)
  
  #remap the nodes
  treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
  treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]
  
  
  ii=masterTreePaths$matIndex[(treePaths$nodeId[,2]-1)*nrow(masterTreePaths$matIndex)+treePaths$nodeId[,1]]
  
  vals=double(length(masterTreePaths$dist))
  vals[]=NA
  if(sum(is.na(ii))>0 & !is.null(i)) {
    message("error: discordant tree topology in tree", i)
    return(vals)
  }
  vals[ii]=treePaths$dist
  vals
}
allPaths=function(tree, categorical = F){
  if (!categorical){
    dd=dist.nodes(tree)
  }
  allD=double()
  nn=matrix(nrow=0, ncol=2)
  nA=length(tree$tip.label)+tree$Nnode
  matIndex=matrix(nrow=nA, ncol=nA)
  index=1
  for ( i in 1:nA){
    ia=getAncestors(tree,i)
    if(length(ia)>0){
      if(categorical) {
        # add the state of node i to allD length(ia) times
        x = which(tree$edge[,2] == i)
        state = tree$edge.length[x]
        allD = c(allD, rep(state, length(ia)))
      }
      else {
        allD=c(allD, dd[i, ia])
      }
      nn=rbind(nn,cbind(rep(i, length(ia)), ia))
      for (j in ia){
        matIndex[i,j]=index
        index=index+1
      }
    }
  }
  return(list(dist=allD, nodeId=nn, matIndex=matIndex))
}
matchAllNodes=function(tree1, tree2){
  map=matchNodesInject(tree1,tree2)
  map=map[order(map[,1]),]
  map
}
matchNodesInject=function (tr1, tr2){
  if(length(tmpsp<-setdiff(tr1$tip.label, tr2$tip.label))>0){
    #stop(paste(paste(tmpsp, ","), "in tree1 do not exist in tree2"))
    stop(c("The following species in tree1 do not exist in tree2: ",paste(tmpsp, ", ")))
  }
  commontiplabels <- intersect(tr1$tip,tr2$tip)
  if(RF.dist(pruneTree(tr1,commontiplabels),pruneTree(tr2,commontiplabels))>0){
    stop("Discordant tree topology detected - gene/trait tree and treesObj$masterTree have irreconcilable topologies")
  }
  #if(RF.dist(tr1,tr2)>0){
  #  stop("Discordant tree topology detected - trait tree and treesObj$masterTree have irreconcilable topologies")
  #}
  
  toRm=setdiff(tr2$tip.label, tr1$tip.label)
  desc.tr1 <- lapply(1:tr1$Nnode + length(tr1$tip), function(x) extract.clade(tr1,
                                                                              x)$tip.label)
  names(desc.tr1) <- 1:tr1$Nnode + length(tr1$tip)
  desc.tr2 <- lapply(1:tr2$Nnode + length(tr2$tip), function(x) extract.clade(tr2,
                                                                              x)$tip.label)
  names(desc.tr2) <- 1:tr2$Nnode + length(tr2$tip)
  Nodes <- matrix(NA, length(desc.tr1), 2, dimnames = list(NULL,
                                                           c("tr1", "tr2")))
  for (i in 1:length(desc.tr1)) {
    Nodes[i, 1] <- as.numeric(names(desc.tr1)[i])
    for (j in 1:length(desc.tr2)) if (all(desc.tr1[[i]] %in%
                                          desc.tr2[[j]]))
      Nodes[i, 2] <- as.numeric(names(desc.tr2)[j])
  }
  
  iim=match(tr1$tip.label, tr2$tip.label)
  Nodes=rbind(cbind(1:length(tr1$tip.label),iim),Nodes)
  if(any(table(Nodes[,2])>1)){
    stop("Incorrect pseudorooting detected - use fixPseudoroot() function to correct trait tree topology")
  }
  
  Nodes
}
scaleMat=function(mat){t(apply(mat,1,scaleDist))}
scaleDist=function(x){
  x/sqrt(sum(x^2))
}

#rerpath = find.package('RERconverge')
#toytreefile = "subsetMammalGeneTrees.txt" 
#tree = paste(rerpath,"/extdata/",toytreefile,sep="")
#
#b = breadTrees1(tree)
#
#subset = b$masterTree$tip.label[7:62]
#b2 = breadTrees1(tree,useSpecies=subset)
