library('DNAcopy')
library('optparse')
option_list <- list( 
    make_option(c("-s", "--smooth"), 
                action="store_true", 
                default=FALSE, 
                help="Smooth the data before segmentation"),
    make_option(c("-t", "--threshold"), 
                action="store", 
                default=0.01, 
                help="Significance threshold"),
    make_option(c("-n", "--nperm"), 
                action="store", 
                default=10000, 
                help="Number of permutations"),
    make_option(c("-m", "--perm_method"), 
                action="store", 
                default='perm', 
                help="P. value calculation method"),
    make_option(c("-u", "--undo_split"), 
                action="store", 
                default="none", 
                help="How change-points have to be undone"),
    make_option(c("-i", "--min_width"), 
                action="store", 
                default=2, 
                help="Minimum number of bins for a changed segment")
)
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments = 2)

input.path <- args$args[1]
output.path <- args$args[2]
smooth_cbs <- args$options$smooth
threshold <- args$options$threshold
nperm <- args$options$nperm
perm_method <- args$options$perm_method
undo_split <- args$options$undo_split
min.width <- args$options$min_width


write(paste0("Input: ", input.path), stderr())
write(paste0("Output: ", output.path), stderr())
write(paste0("Smoothing: ", smooth_cbs), stderr())
write(paste0("Threshold: ", threshold), stderr())
write(paste0("N. permutations: ", nperm), stderr())
write(paste0("Permutation method: ", perm_method), stderr())
write(paste0("Undo split: ", undo_split), stderr())
write(paste0("Min. width: ", min.width), stderr())

write("Loading probe coverages into a data frame", stderr())
tbl = read.delim(input.path)
if (!is.null(tbl$weight)) {
    # Drop any 0-weight bins
    tbl = tbl[tbl$weight > 0,]
}
cna = CNA(cbind(tbl$ratio), tbl$chr, tbl$start,
          data.type="logratio", presorted=T)

write("Segmenting the probe data", stderr())
set.seed(0xA5EED)

# additional smoothing (if --smooth-cbs provided)
if (smooth_cbs) {
	write("Performing smoothing of the data", stderr())
	cna = smooth.CNA(cna)
}

if (is.null(tbl$weight)) {
    fit = segment(cna, 
                  alpha=threshold, 
                  nperm = nperm, 
                  p.method = perm_method, 
                  undo.splits = undo_split, 
                  min.width = min.width)
} else {
    fit = segment(cna, 
                  weights=tbl$weight, 
                  alpha=threshold, 
                  nperm = nperm, 
                  p.method = perm_method, 
                  undo.splits = undo_split,
                  min.width = min.width)
}

write("Printing the CBS table to standard output", stderr())
write.table(fit$output, output.path, sep='\t', row.names=FALSE)




