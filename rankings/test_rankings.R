# testing rankings

library(readr)
library(data.table)
library(HiCcompare)
library(dplyr)

# Set up
#####
centromeres <- read.table("D:/brain/data/hg38_centromeres.bed")

# for dplfc 1 (dplfc folder)
dplfc_1_mat <- read_tsv("D:/brain/data/all_resolutions/dplfc_1000000.matrix", col_names = FALSE) %>% as.data.table()
dplfc_1_bed <- read_tsv("D:/brain/data/all_resolutions/dplfc_1000000_abs.bed", col_names = FALSE) %>% as.data.table()

# for amyg
amyg_mat <- read_tsv("D:/brain/data/all_resolutions/amyg_1000000.matrix", col_names = FALSE) %>% as.data.table()
amyg_bed <- read_tsv("D:/brain/data/all_resolutions/amyg_1000000_abs.bed", col_names = FALSE) %>% as.data.table()


dplfc_1 <- hicpro2bedpe(dplfc_1_mat, dplfc_1_bed)$cis[-c(23,25)] # remove chrM and chrY
amyg <- hicpro2bedpe(amyg_mat, amyg_bed)$cis[-c(23,25)]

amyg_dplfc1 <- mapply(create.hic.table, amyg, dplfc_1, MoreArgs = list(exclude.regions = centromeres, exclude.overlap = 0,
                                                                       scale = FALSE), SIMPLIFY = FALSE)

amyg_dplfc1 <- hic_loess(amyg_dplfc1, Plot = TRUE)

amyg_dplfc1 <- hic_diff(amyg_dplfc1, Plot = TRUE, diff.thresh = 'auto')


#######
# Test rankings

# MD plot of p-values
hic.table <- amyg_dplfc1[[1]]
MD.plot2(hic.table$adj.M, hic.table$D, hic.table$p.value)

# set up
idx <- 1:(nrow(hic.table) * alpha)
alpha = 0.05

# get top 5% of max ranks where max is taken from M, raw difference, and Avg expression
hic.table <- hic.table[order(rnkMax),]
topRanks <- rep(1, nrow(hic.table)) # make indicator for top ranks
topRanks[idx] <- 0 # set top ranking rows to 0 indicator for plotting on MD plot
MD.plot2(hic.table$adj.M, hic.table$D, p.val = topRanks)

# get top 5% of mean ranks
hic.table <- hic.table[order(rnkMean),] # order by mean rank
topRanks <- rep(1, nrow(hic.table)) # make indicator for top ranks
topRanks[idx] <- 0 # set top ranking rows to 0 indicator for plotting on MD plot
MD.plot2(hic.table$adj.M, hic.table$D, p.val = topRanks)


# by top 5% of rnkM and rnkDiff
hic.table <- hic.table[order(rnkM, rnkDiff),] # order by mean rank
topRanks <- rep(1, nrow(hic.table)) # make indicator for top ranks
topRanks[idx] <- 0 # set top ranking rows to 0 indicator for plotting on MD plot
MD.plot2(hic.table$adj.M, hic.table$D, p.val = topRanks)

# by top 5% of rnkDiff
hic.table <- hic.table[order(rnkDiff),] # order by mean rank
topRanks <- rep(1, nrow(hic.table)) # make indicator for top ranks
topRanks[idx] <- 0 # set top ranking rows to 0 indicator for plotting on MD plot
MD.plot2(hic.table$adj.M, hic.table$D, p.val = topRanks)


get_ranks <- function(hic.table, alpha = 0.1) {
  idx <- 1:(nrow(hic.table) * alpha)
  topRanks <- rep(1, nrow(hic.table)) # make indicator for top ranks
  topRanks[idx] <- 0 # set top ranking rows to 0 indicator for plotting on MD plot
  MD.plot2(hic.table$adj.M, hic.table$D, p.val = topRanks)
  
}