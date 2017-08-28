# # figure out excluding regions
#
# # prostate_analysis_final.rmd
# exclude.regions <- exclude.gi[findOverlaps(fish_genes[3,], exclude.gi)@to,]
# iset <- fish_genes
#
# findOverlaps(iset@regions, exclude.regions)
# sum(width(pintersect(exclude.regions, iset@regions[1,])))
# sum(width(pintersect(exclude.regions, iset@regions[4,])))
#
# olaps <- GenomicRanges::findOverlaps(exclude.regions, iset@regions[1,])
# olap_length <- width(ranges(olaps, ranges(exclude.regions), ranges(iset@regions[1,])))
#
# q <- GRanges(seqnames = rep('chr1', 2), ranges = IRanges(start = c(1, 101), end = c(100, 200)))
# exclude <- GRanges(seqnames = rep('chr22', 3), ranges = IRanges(start = c(80, 1, 16500000), end = c(120, 20, 16900000)))
#
# width(pintersect(exclude, q))
#
# # think i need to go through all regions in matrix individual and sum width of pintersect of the excluded regions
#
# data('HMEC.chr22')
# data('NHEK.chr22')
# hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
#
# iset <- make_InteractionSet(hic.table)
# iset@regions
#
# exclude.overlap <- 0.5
#
# # first split iset into list containing each region as an individual GRanges
# target_regions <- split(iset@regions, as.factor(iset@regions))
# # get widths for each target region
# target_widths <- width(iset@regions)
# # get number of basepairs overlapping with excluded regions
# target_overlap <- sapply(target_regions, function(x) {
#                         sum(width(pintersect(exclude.regions, x)))
#                       })
# # get percentage of overlap
# percent_olap <- target_overlap / target_widths
# # get which regions to remove
# regions_to_remove <- which(percent_olap >= exclude.overlap)
# if (length(regions_to_remove) > 0) {
#   new_exclude <- iset@regions[regions_to_remove,]
#   to_remove <- InteractionSet::findOverlaps(new_exclude, iset)
# }
