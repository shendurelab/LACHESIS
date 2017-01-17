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

jpeg.file <- "out/heatmap.jpg"
ggsave( filename=jpeg.file, plot=p, width=7, height=6 )
