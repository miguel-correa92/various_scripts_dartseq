## read and organize STRUCTURE results
## Code modified from Rrunstruct R package
## https://github.com/eriqande/Rrunstruct

library(tidyverse)
library(combinat)

# read_and_process_structure_output

# Read and process logs

read_traces <- function (D, pattern){
  files <- dir(D, pattern = pattern, full.names = TRUE, recursive = TRUE)
  
  lapply(files, function(x) {
    fn <- str_replace(x, "^.*_k", "") %>% str_replace('\\.log', '')
    K <- str_split(fn, "_")[[1]][1]
    rep <- str_split(fn, "_")[[1]][2]
    slurp_traces(x, K = K, rep = rep) %>% 
      pivot_longer(cols = -(K:Sweep), names_to = 'variable', values_to = 'value')
  }) %>% do.call(rbind, .)
}

slurp_traces <- function (file, K, rep) {
  x <- readLines(file)
  x <- x[which(str_detect(x, "Finished initialization; starting MCMC")):
           which(str_detect(x, "MCMC completed"))]
  x2 <- x[str_detect(x, "^ *[0-9]+:")] %>% 
    str_replace_all(":", "") %>% 
    str_trim %>% 
    str_replace_all(" +0$", "")
  
  header <- x[str_detect(x, "Rep#:")] %>% 
    unique %>% .[1] %>% 
    str_replace_all("Ln Like", "LnLike") %>% 
    # str_replace_all('Est Ln P(D)"', 'EstLnP(D)') %>% 
    str_replace_all("Rep#:", "Sweep") %>% 
    str_trim
 
  df <- as_tibble(cbind(K = K, Rep = rep, 
                     read.table(textConnection(c(header, x2)),
                                header = TRUE, na.strings = "--")
                     )
               )
  return(df)
}


#' make the trace data frame tidy and long format
gather_traces <- function(df) {
  gather(df, "variable", "value", -(K:Sweep))
}

# read and process results


read_results <- function (D, pattern) {
  files <- dir(D, pattern = pattern, full.names = TRUE, recursive = TRUE)
  lapply(files, function(x) {
    # x <- files[1]
    fn <- str_replace(x, "^.*_k", "") %>% str_remove('_f')
    K <- str_split(fn, "_")[[1]][1]
    Rep <- str_split(fn, "_")[[1]][2]
    cbind(fn, K, Rep)
    
    slurp_results(x, K = K, Rep = Rep) %>% gather_results
  }) %>% do.call(rbind, .)
}

# a <- slurp_results(file = files, K = K, Rep = Rep)

slurp_results <- function(file, K, Rep) {
  # x <- files[1]
  x <- suppressWarnings(readLines(file))
  start <-which(str_detect(x, "^Inferred ancestry of individuals"))
  end <- which(str_detect(x, "^Estimated Allele Frequencies"))
  
  # get the elements of the header
  header <- x[start+1] %>%
    str_replace_all("[^a-zA-Z ]", "") %>%
    paste("Index", .) %>%
    str_replace("Inferred clusters", "") %>%
    str_trim %>%
    str_split(" +") %>%
    "[["(1)
  
  
  # then extract the output for each individual
  indivs <- x[(start+2):end]
  
  # then get that all in a data frame
  tmp <- indivs[str_detect(indivs, "[0-9]")] %>% # toss blank lines
    str_replace_all("[():]", "")  %>%  # toss the ()'s and the :'s from the output
    str_trim
  
  # read it into a table
  df <- read.table(textConnection(tmp), header = FALSE)
  
  # then make the header as appropriate
  names(df) <- c(header, 1:(ncol(df)-length(header)) )
  
  as_tibble(cbind(K = K, Rep = Rep, df))
}


#' makes a long format data frame of the results
gather_results <- function(df) {
  
  pivot_longer(df, cols = -c(K:Miss), names_to = 'cluster', values_to = 'probability')
  
}


MAP_cluster  <- function(df) {
  df %>%
    group_by(pick(!c(cluster, probability))) %>% # group by everything except cluster and probability
    filter(probability == max(probability)) %>%   #%>% # get the row that has maximum posterior
    summarise(across(everything(), .fns = first)) %>% 
    # summarise_each(funs(first)) #%>% # take only one if the max posterior is shared between clusters
    transmute(map_cluster = cluster) %>%
    ungroup
}

relabel_map_clusters <- function(df) {
  # df <- res2
  # get the map_clusters from rep 1 and attach them to the original data frame
  # and then for each rep > 1 determine which permutation gives the closest correspondence
  # with rep 1, and put that on there
  df2 <- df %>%
    filter(Rep == 1) %>%
    mutate(rep1_map_cluster = map_cluster) %>%
    select(-map_cluster, -Rep) %>%
    ungroup %>%
    left_join(df, .)
  
  df3 <- df2 %>%  # slap those values back onto the original data frame
    group_by(K, Rep) %>%  # now, group it by K and rep and compute the permutation that makes each rep look most like rep 1
    summarise(max_perm = find_max_perm(rep1_map_cluster, map_cluster, K)) %>%
    left_join(df2, .) %>%
    group_by(K, Rep) %>%
    mutate(map_cluster_relabeled = relabel_by_max_perm(rep1_map_cluster, map_cluster, K))
  
  df3
  
  
}

# this is a window function
relabel_by_max_perm <- function(c1, c2, k) {
  k <- as.numeric(as.character(k[1]))
  permy <- permn(k)
  scores <- sapply(permy, function(x) sum(c1 == x[as.numeric(c2)]))
  maxp <- which.max(scores)
  permy[[maxp]][as.numeric(c2)]
  
}

# this is a summary function
find_max_perm <- function(c1, c2, k) {
  k <- as.numeric(as.character(k[1]))
  permy <- permn(k)
  scores <- sapply(permy, function(x) sum(c1 == x[as.numeric(c2)]))
  which.max(scores)
}




fix_fs <- function(v, mp, k) {
  
  mp <- as.numeric(as.character(mp[1]))
  k <- as.numeric(as.character(k[1]))
  
  
  
  perm <- permn(k)[[mp]]
  new.v <- v
  fs <- str_detect(new.v, "^F")
  nums <- as.numeric(str_replace_all(new.v[fs], "F", "") )
  new.nums <- perm[nums]
  new.v[fs] <- paste("F", new.nums, sep = "")
  new.v
}


relabel_traces <- function(tr, max_perms) {
  tr %>%
    inner_join(ungroup(max_perms)) %>%
    group_by(K, Rep) %>%
    mutate(variables_relabeled = fix_fs(as.character(variable), mp, K))
}



# now we have to relabel the results.
fix_clust <- function(cl, mp, k) {
  mp <- as.numeric(as.character(mp[1]))
  k <- as.numeric(as.character(k[1]))
  cl <- as.numeric(as.character(cl))
  perm <- permn(k)[[mp]]
  
  perm[cl]
  
}


relabel_results <- function(res, max_perms) {
  res %>%
    inner_join(ungroup(max_perms)) %>%
    group_by(K, Rep) %>%
    mutate(cluster_relabeled = fix_clust(cluster, mp, K))
}
