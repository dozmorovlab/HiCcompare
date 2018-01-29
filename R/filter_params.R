#' Determine the A quantile cutoff to be used
#' 
#' 
#' 


filter_params <- function(hic.table) {
  # split up table by distances
  temp_list <- S4Vectors::split(hic.table, hic.table$D)
  # randomize IFS
  temp_tables <- lapply(temp_list, randomize_IFs)
  # combine tables back into 1
  new.table <- data.table::rbindlist(temp_tables)
  # order tables
  hic.table <- hic.table[order(start1, start2),]
  new.table <- new.table[order(start1, start2),]
  # create new tables to have differences added
  sparse1 <- cbind(hic.table$start1, hic.table$start2, hic.table$IF1)
  sparse2 <- cbind(new.table$start1, new.table$start2, new.table$IF1)
  sparse3 <- cbind(hic.table$start1, hic.table$start2, hic.table$IF2)
  sparse4 <- cbind(new.table$start1, new.table$start2, new.table$IF2)
  shuffle1 <- create.hic.table(sparse1, sparse2, chr = hic.table$chr1[1])
  shuffle2 <- create.hic.table(sparse3, sparse4, chr = hic.table$chr1[1])
  # add changes
  
  # normalize
  shuffle1 <- hic_loess(shuffle1)
  suffle1 <- hic_compare(shuffle1, Plot = T)
}


# function to shuffle the IFs within a table containing only entries at a single distance
randomize_IFs <- function(dist_table) {
  # shuffle vector of IFs
  newIF1 <- sample(dist_table$IF1, size = length(dist_table$IF1), replace = FALSE)
  newIF2 <- sample(dist_table$IF2, size = length(dist_table$IF2), replace = FALSE)
  # create new hic.table with new IF vectors
  sparse1 <- cbind(dist_table$start1, dist_table$start2, newIF1)
  sparse2 <- cbind(dist_table$start1, dist_table$start2, newIF2)
  temp.table <- create.hic.table(sparse1, sparse2, chr = dist_table$chr1[1])
  return(temp.table)
}
