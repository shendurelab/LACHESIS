#!/usr/bin/env Rscript
# The above "shebang" allows this file to be self-executing
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



# QuickDotplot.R


library( ggplot2, quietly=TRUE ) # ggplot
library( scales, quietly=TRUE) # labels=comma

arglist <- commandArgs(trailingOnly=TRUE)
data.file <- arglist[[1]]
jpheg.file <- 'out/dotplot.renumbered.jpg'



# theme_black: Copied from http://stackoverflow.com/questions/13999103/control-color-of-legend-elements-that-are-not-colour-guides-in-ggplot
theme_black <- function (base_size = 16, base_family = ""){
    theme_minimal() %+replace%
        theme(
              line = element_line(colour = "white", size = 0.5, linetype = 1,
                        lineend = "butt"),
              rect = element_rect(fill = "white",
                        colour = "white", size = 0.5, linetype = 1),
              text = element_text(family = base_family,
                        face = "plain", colour = "white", size = base_size,
                        angle = 0, lineheight = 0.9, hjust = 0.5, vjust = 0.5),
              plot.background = element_rect(colour = 'black', fill = 'black'),
              plot.title = element_text(size = rel(1.2)),
              panel.border = element_rect(fill = NA, colour = "white"),
              panel.grid.major = element_line(colour = "grey20", size = 0.2),
              panel.grid.minor = element_line(colour = "grey5", size = 0.5),
              strip.background = element_rect(fill = "grey30", colour = "grey30")
             )
    }





# Read the data file.
data <- read.table( data.file, header=FALSE )

# Re-order the chromosomes.  This factoring should work for chromosomes whether or not they're prepended with 'chr', though it hasn't been test for with 'chr'.
#data$V3 <- factor( data$V3, levels = c( "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY" ) )


# Load the data into a ggplot object.
p <- ggplot( data, aes( x=V1, y=V2 ) )
p <- p + theme_black()
p <- p + theme_bw()
#p <- p + theme( axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) # remove some lines for Nat. Biotech

# Plot the points, and load in the colors, if any were given.
p <- p + geom_point( size = 2.5, aes( col=data[,3] ) )

# Apply the color scheme for SKY.
# These RGB values were taken from a close-up of an actual SKY image, sampled in MS Paint like a boss: http://www.ncbi.nlm.nih.gov/sky/ccap_helper.cgi?tsc=2
# Chr9 has been changed from white to black because I'm using a white background.
p <- p + scale_colour_manual( values=rev( c('#FFE800','#900101','#8F8090','#8DE3F6','#8A6911','#9C4D94','#F0A0B0','#F66824','#000000','#049618','#0495D1','#EE04EF','#F70914','#9386EE','#7FE7A0','#F8A80C','#0028F8','#C01C68','#48F41B','#B064EF') ), name = "Chromosome" )



# Increase font size on axis tick labels.
p <- p + theme( axis.text.x  = element_text(size=16) )
p <- p + theme( axis.text.y  = element_text(size=16) )
#p <- p + theme( axis.title.x = element_text(size=26,vjust=-0.1) )
#p <- p + theme( axis.title.y = element_text(size=26,vjust=0.3) )

# Remove text labels.  We don't want any text on a JPEG image.
p <- p + ggtitle("") + xlab("") + ylab("")
p <- p + theme( legend.position="none" )
p <- p + scale_x_continuous( breaks=seq(0,8000,1000), labels = comma )

# Save the plot to the jpeg file.  Adjust the width to make room for a legend, if there is one.
suppressWarnings( ggsave( filename=jpeg.file, plot=p, width=8, height=5 ) )
