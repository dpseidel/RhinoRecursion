pairwise_ks_test <- function(value, group, alternative = "two.sided"){
  

  
  f <- function(x, y){ 
    p <- ks.test(x, y, alternative = alternative, exact = 
                   NULL)$p.value 
    return(p) 
  }
  
  res <- suppressWarnings(lapply(lst, function(x) lapply(lst, function(y) f(x, y)))) 
  res <- unlist(res)
  res <- res[res!=1] %>% unique()
  res <- p.adjust(res, method = "bonferroni")
  res <- matrix(res, nrow = length(lst), ncol = length(lst), byrow = T)
  row.names(res) <- colnames(res) <- names(lst)
  cat("Pairwise Kolmogorov-Smirnov Test with Correction: p-value Matrix","\n","\n")
  return(res)
}


#dumbest way i can figure this out tonight... 
expand.grid(1:4, 1:4) %>% filter(Var1 != Var2) -> x
x[x[,2]>=x[,1],] 

comb <- expand.grid(unique(levels(df$ToD)), unique(levels(df$ToD))) %>% 
  filter(Var1 != Var2) %>% .[x[,2]>=x[,1],] %>% mutate_all(as.character())


ToD_disp <- filter(lag4, !is.na(dt4)) %>%
  mutate(disphr = dt4/time4) %>% 
  select(ToD, disphr) %>% 
  split(.$ToD, ) %>%
  map(~pull(.x, disphr))

results <- pmap_dbl(list(a = comb$Var1, b = comb$Var2), 
     function(a,b){
  suppressWarnings(ks.test(ToD_disp[[a]], ToD_disp[[b]])$p.value)
}) %>% p.adjust(method = "bonferroni") 

names(results) <- paste(comb$Var1, comb$Var2, sep =".")
