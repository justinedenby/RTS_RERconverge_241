#PURPOSE: helper function for running RERconverge with trees created from the breadTrees1() function -- ensures that the tree is preordered, which is required for RERconverge

library(TreeTools)

patchTreesObjToPreorder <- function(obj) {
  ord <- attr(obj$masterTree, "order")
  if (!is.null(ord) && ord == "preorder") { message("Already preorder."); return(obj) }
  
  master_old <- obj$masterTree
  master_new <- Preorder(master_old)
  
  # descendant tip-sets -> old-node -> new-node correspondence
  descTips <- function(tree) {
    nTip <- length(tree$tip.label)
    tr <- reorder(tree, "postorder")
    res <- vector("list", nTip + tree$Nnode)
    for (t in seq_len(nTip)) res[[t]] <- tree$tip.label[t]
    for (e in seq_len(nrow(tr$edge))) {
      p <- tr$edge[e, 1]; ch <- tr$edge[e, 2]
      res[[p]] <- c(res[[p]], res[[ch]])
    }
    vapply(res, function(x) paste(sort(unique(x)), collapse = "\x1f"), character(1))
  }
  old2new <- match(descTips(master_old), descTips(master_new))
  stopifnot(!anyNA(old2new), !anyDuplicated(old2new))   # errors if rooting differs
  
  # column layouts (deterministic from topology)
  ap_old <- allPathsTrackBranches(master_old)
  ap_new <- allPathsTrackBranches(master_new)
  mi_new <- ap_new$matIndex
  
  dk <- old2new[ap_old$destinNode]
  ak <- old2new[ap_old$ancNode]
  perm <- mi_new[cbind(dk, ak)]                 # perm[k] = new position of old column k
  stopifnot(!anyNA(perm), !anyDuplicated(perm), length(perm) == ncol(obj$paths))
  
  reCol <- function(M) {
    if (is.null(M)) return(NULL)
    stopifnot(ncol(M) == length(perm))
    Mn <- M; Mn[, perm] <- M
    cn <- colnames(M); if (!is.null(cn)) { nc <- character(length(perm)); nc[perm] <- cn; colnames(Mn) <- nc }
    Mn
  }
  
  obj$paths        <- reCol(obj$paths)
  obj$weights      <- reCol(obj$weights)
  obj$pathsImputed <- reCol(obj$pathsImputed)
  obj$masterTree   <- master_new
  obj$ap           <- ap_new
  obj$matIndex     <- mi_new
  ma <- (mi_new > 0) + 1 - 1; ma[is.na(ma)] <- 0; obj$matAnc <- ma
  obj
}
