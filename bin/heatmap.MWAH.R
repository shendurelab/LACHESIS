#!/usr/bin/env Rscript
#///////////////////////////////////////////////////////////////////////////////
#//                                                                           //
#// This software and its documentation are copyright (c) 2014-2015 by Joshua //
#// N. Burton and the University of Washington.  All rights are reserved.     //
#//                                                                           //
#// THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS  //
#// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF                //
#// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  //
#// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY      //
#// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT //
#// OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR  //
#// THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                //
#//                                                                           //
#///////////////////////////////////////////////////////////////////////////////


# Make a colored hex plot from heatmap.txt.
#
# Josh Burton
# December 2012


library(ggplot2,quietly=TRUE) # ggplot
library(RColorBrewer) # colorRampPalette
library(grid) # unit

heatmap.file <- 'heatmap.txt'
jpeg.file <- "out/HiC_heatmap.jpg"
shuffle <- 0 # if 1, randomly permute the heatmap to show the 'pre-Lachesis' dataset.


########## Get a ColorBrewer palette.

# Orange color scheme for S. stipitis image in MetaPhase fig. 2c
#pal <- colorRampPalette( c('white', '#FC4E2A', 'black' ), bias=3 )(9)
# Magenta color scheme for K. wickerhamii image in MetaPhase fig. 2d
#pal <- colorRampPalette( c( brewer.pal( 5, 'PuRd' ), 'black'), bias=5 )(9)

# Reds.  Decomment the first line, then one of the following three lines.
pal <- brewer.pal( 9, "Reds" ) # i.e., c( "#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D" )
# Decomment ONLY ONE of the following 3 lines.  These palettes are (progressively darker) darkened versions that interpolate at the dark end of the spectrum.
#pal <- c( "white", "#FFEBE1", pal[2:8], "#860811", pal[9], "#400000", "#100000" )
pal <- c( pal[1:8], "#950B13", "#860811", "#76040F", pal[9], "#500000", "#400000", "#300000", "#200000", "#100000", "#000000" )
#pal <- c( pal[1], pal[3], pal[5], pal[7], pal[8], "#950B13", "#860811", pal[9], "#500000", "#400000", "#300000", "#200000", "#100000", "#000000" ) # darkest red: for MetaPhase Fig. 2c,d

# BEST: Yellow-to-orange-to-red palette.  Try bias = 1,2,3 for progressively darkened versions.
pal <- colorRampPalette( c( brewer.pal( 9, 'YlOrRd' ), 'black'), bias=2 )(9)

##########


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

# Create a plot of this heatmap.  If the data values are all non-negative, log=transform them.
if ( min( heatmap$Z ) >= 0 ) {
    p <- ggplot( heatmap, aes( x=X, y=Y, fill=log10(Z+1) ) )
} else {
    p <- ggplot( heatmap, aes( x=X, y=Y, fill=Z ) )
}


if ( shuffle == 1 ) {
    # Load shuffled data.  Note that the shuffle order is 1-indexed so we'll need to twiddle the data to use it properly.
    p <- ggplot( heatmap, aes( x=shuffle.order[X+1]-1, y=shuffle.order[Y+1]-1, fill=log10(Z+1) ) )
}


# Set the theme for this plot.
p <- p + theme_bw()
p <- p + geom_tile() # choose a rectangular tiling
p <- p + xlab("") # clear x- and y-axis labels
p <- p + ylab("")
#p <- p + scale_fill_gradientn( colours = palette(100), name="log link density" )

# Stub code for variable tile sizes (note: this requires dim(heatmap) == N*N, so there can't be any absent data points e.g. on the main diagonal)
#cell.size <- rep(1,N)
#cell.size[300] = 20
#cell.size[500] = 40
#cell.boundary <- rep(cell.size,N)
#p <- p + geom_tile( aes( width=cell.boundary, height=cell.boundary ) )


# Add vertical lines to delineate chromosomes.
if ( shuffle == 0 ) {
    #p <- p + geom_vline( xintercept = chrom.borders, col='black', linetype='dashed', size=.5 )
    #p <- p + geom_hline( yintercept = chrom.borders, col='black', linetype='dashed', size=.5 )
    p <- p + geom_vline( xintercept = chrom.borders, col='black', linetype='dashed', size=1 )
    p <- p + geom_hline( yintercept = chrom.borders, col='black', linetype='dashed', size=1 )
}


# Tweak the theme for better readability on a slide (November 2013 Research Report).
p <- p + theme( axis.text.x = element_text( size=0 ), axis.text.y = element_text( size=0 ) ) # hide axis labels
p <- p + theme( axis.ticks = element_blank() ) # hide tick marks

# Tweak the theme to make a high-quality image for publication
#p <- p + theme( axis.text.x = element_text( size=20 ), axis.text.y = element_text( size=20 ) )
p <- p + theme( legend.text=element_text(size=20), legend.key.size=unit(1,"cm"), legend.key.width=unit(1.5,"cm") )
p <- p + scale_fill_gradientn( colours = palette(100), name="" ) # same colors, but remove legend name

# Set width and height of output file.  The aspect ratio should be tweaked so that the heatmap portion of the image is an exact square.
w <- 11.45
h <- 10
ggsave( filename=jpeg.file, plot=p, width=w, height=h )
