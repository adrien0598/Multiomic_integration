runLimma <- function(measurements, targets, comparisons = NULL)
{
  if (!is.null(comparisons)){
      cont.matrix <- matrix(c(-1, 1), ncol = 1)
      
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

plot_top <- function(d, n_top, score){
  colnames(d) = c("ID", "score")
  top_up <- slice_head(d, n = n_top)
  top_down <- slice_tail(d, n = n_top)
  
  p <- bind_rows(list(up = top_up, down = top_down), .id = "status") %>%
    mutate(ID = fct_inorder(ID)) %>%
    ggplot(aes(x = score, y = ID, fill = status)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
    ylab(score)+
    theme_prism()
  
  return(p)
}