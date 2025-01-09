makeContrastsAlt <- function(targets, comparisons)
{
  cont.matrix <- matrix(0,nrow = length(unique(targets$condition)), ncol = length(comparisons))
  i <- 1
  for (comparison in comparisons)
  {
    for (j in 1:length(comparison))
    {
      cont.matrix[abs(comparison[j]),i] <- cont.matrix[abs(comparison[j]),i]+(comparison[j]/abs(comparison[j]))
    }
    i <- i + 1
  }
  return(cont.matrix)
}

runLimma <- function(measurements, targets, comparisons = NULL)
{
  if (!is.null(comparisons)){
      cont.matrix <- makeContrastsAlt(targets, comparisons)
      
      cont.matrix <- as.data.frame(cont.matrix)
      row.names(cont.matrix) <- unique(targets$condition)
      cont.matrix <- as.matrix(cont.matrix)
      
      fcond <- factor(targets$condition, levels = unique(targets$condition))
      
      design <- model.matrix(~0+fcond)
      design <- as.data.frame(design)
      names(design) <- unique(targets$condition)
      design <- as.matrix(design)
      
      fit <- lmFit(measurements, design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      
      return(list(fit2, cont.matrix, fit))
  }
  else{
    print("error, no comparison")
    return(NULL)
  }
}

ttopFormatter <- function(ttop)
{
  ttop$ID <- row.names(ttop)
  ttop <- ttop[,c(7,1,2,3,4,5,6)]
  ttop <- ttop[complete.cases(ttop),]
  return(ttop)
}
