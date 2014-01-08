#!/net/gs/vol3/software/modules-sw/R/2.15.1/Linux/RHEL6/x86_64/bin/Rscript

# Make a colored hex plot from heatmap.txt.
#
# Josh Burton
# December 2012


library(ggplot2,quietly=TRUE) # ggplot
library(RColorBrewer) # colorRampPalette
library(grid) # unit

heatmap.file <- 'heatmap.txt'
jpeg.file <- "~/public_html/heatmap.MWAH.jpg"
shuffle <- 0 # TEMP: if 1, randomly permute the heatmap to show the 'pre-Lachesis' dataset.

# Get a ColorBrewer palette.
# Spectrum: blue = no links; yellow = some links; red = lots of links
pal <- c( "#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2" ) # this is brewer.pal( 11,"Spectral" )
pal <- c( "#9E0142","#B91F48","#D53E4F","#E45549","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2" ) # added 2 interpolated colors at the red end of the spectrum
pal <- rev(pal)

# Reds!
pal <- brewer.pal( 9, "Reds" ) # e.g., c( "#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D" )
pal <- c( pal[1:7], pal[8], "#76040F", pal[9], "#500007", "#300000" ) # darken overall palette by interpolating at dark end of spectrum
palette <- colorRampPalette( pal )


# Read the file.
heatmap <- read.table( heatmap.file, header=TRUE )

# Make a random permutation of the rows.
if ( shuffle == 1 ) {
    N <- max( heatmap[,1] ) + 1
    shuffle.order <- order( runif(N) )
}

# Read the chrom_borders file made by the function MakeWholeAssemblyHeatmap().
chrom.borders <- as.numeric( scan( 'heatmap.chrom_breaks.txt', quiet=TRUE ) )

# Create a plot of this heatmap.
p <- ggplot( heatmap, aes( x=X, y=Y, fill=log10(Z+1) ) ) # load data
if ( shuffle == 1 ) {
    # Load shuffled data.  Note that the shuffle order is 1-indexed so we'll need to twiddle the data to use it properly.
    p <- ggplot( heatmap, aes( x=shuffle.order[X+1]-1, y=shuffle.order[Y+1]-1, fill=log10(Z+1) ) )
}
p <- p + theme_bw()
p <- p + geom_tile() # choose a rectangular tiling
p <- p + xlab("") # clear x- and y-axis labels
p <- p + ylab("")
p <- p + scale_fill_gradientn( colours = palette(100), name="log link density" )

# Stub code for variable tile sizes (note: this requires dim(heatmap) == N*N, so there can't be any absent data points e.g. on the main diagonal)
#cell.size <- rep(1,N)
#cell.size[300] = 20
#cell.size[500] = 40
#cell.boundary <- rep(cell.size,N)
#p <- p + geom_tile( aes( width=cell.boundary, height=cell.boundary ) )


# Add vertical lines to delineate chromosomes.
if ( shuffle == 0 ) {
    p <- p + geom_vline( xintercept = chrom.borders, col='black', linetype='dashed', size=.4 )
    p <- p + geom_hline( yintercept = chrom.borders, col='black', linetype='dashed', size=.4 )
}


# Tweak the theme for better readability on a slide (November 2013 Research Report).
#p <- p + theme( axis.text.x = element_text( size=0 ), axis.text.y = element_text( size=0 ) ) # hide axis labels
#p <- p + theme( axis.ticks = element_blank() ) # hide tick marks
#p <- p + theme( legend.text=element_text(size=20), legend.key.size=unit(1,"cm"), legend.key.width=unit(1.5,"cm") )
#p <- p + scale_fill_gradientn( colours = palette(100), name="" ) # same colors, but remove legend name

ggsave( filename=jpeg.file, plot=p, width=11.45, height=10 )
