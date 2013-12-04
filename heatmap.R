#!/net/gs/vol3/software/modules-sw/R/2.15.1/Linux/RHEL6/x86_64/bin/Rscript

# Make a colored heatmap plot from heatmap.txt.
#
# Josh Burton
# December 2012


library(ggplot2,quietly=TRUE) # ggplot
library(RColorBrewer) # colorRampPalette
#library(plyr) # used by library(reshape)
#library(reshape) # melt

heatmap.file <- 'heatmap.txt'


# Get a ColorBrewer palette.
palette <- colorRampPalette( rev( brewer.pal( 11, "Spectral" ) ) )

# Read the file.
heatmap <- read.table( heatmap.file, header=TRUE )

# Melt the heatmap from an NxN table into a 4x(N^2) list of columns.
#heatmap.m <- melt( heatmap, id.vars=c("Row"), na.rm=TRUE )

# Create a plot of this heatmap.
p <- ggplot( heatmap, aes( x=X, y=Y, fill=log10(Z+1) ) ) # load data
#p <- ggplot( heatmap, aes( x=X, y=Y, fill=Z ) ) # load data
p <- p + geom_tile() # choose a rectangular tiling
p <- p + xlab("") # clear x- and y-axis labels
p <- p + ylab("")
p <- p + scale_fill_gradientn( colours = palette(100), name="log10(N links)" ) # choose colors from the palette

jpeg.file <- "~/public_html/heatmap.jpg"
ggsave( filename=jpeg.file, plot=p, width=7, height=6 )
