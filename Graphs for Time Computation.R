####################################################################
# SCRIPT FOR GRAPHICS OF COMPUTING TIMES OF SCCs
####################################################################

install.packages("plotly")
install.packages("gapminder")
library(ggplot2)
library(plotly)
library(gapminder)

setwd("~/MEGA/PhD/12. ICCSA Technical Publication/ICCSA Selected Papers '22/R Code")

grids <- function(axis = c("xy", "x", "y"), color = "grey92", size = NULL, linetype = NULL)
{
  axis <- match.arg(axis)
  grid.major <- element_line(color = color, size = size,
                             linetype = linetype)
  grid.minor <- element_line(color = color, size = 0.25,
                             linetype = linetype)

  switch(axis,
         xy = theme(panel.grid.major = grid.major, panel.grid.minor = grid.minor),
         x = theme(panel.grid.major.x = grid.major, panel.grid.minor.x = grid.minor),
         y = theme(panel.grid.major.y = grid.major, panel.grid.minor.y = grid.minor)
         )
}

###
# For saving:

tiff("test.tiff", units="in", width=8, height=5, res=300)
ggplot(data_triangulation, aes(n, t)) +
       geom_point() + geom_smooth() +
       theme_classic() + grids(linetype = "dashed") +
       xlab("Fineness Degree (n)") + ylab("Computing Time (s)") +
       ggtitle("Computing times for Delaunay Triangulations") + 
       theme(plot.title=element_text(hjust=1, vjust=0.5, face='bold'))

dev.off()

# For interactive plot:
triang <- ggplot(data_triangulation, aes(n, t)) +
       geom_point() + geom_smooth() +
       theme_classic() + grids(linetype = "dashed") +
       xlab("Fineness Degree (n)") + ylab("Computing Time (s)") +
       ggtitle("Computing times for Delaunay Triangulations") + 
       theme(plot.title=element_text(hjust=1, vjust=0.5, face='bold'))
ggplotly(triang)

###

# ONE-SAMPLE TIMES:


data <- read.csv2("one_sample_times.csv", sep = ";", dec = ",", header = T)

one_sample <- ggplot(data = data, aes(x = n,  y = t)) +
  geom_point() +
  geom_smooth() +
  theme_classic() + grids(linetype = "dashed") +
  xlab("Triangulation Fineness Degree (n)") + ylab("Computing Time (s)") +
  ggtitle("Computing times for One-Sample SCC Estimation") + 
  theme(plot.title=element_text(hjust=1, vjust=0.5, face='bold')) +
  theme(legend.position=c(.9,.9))

one_sample
ggplotly(one_sample)
  
## this is just for saving
tiff("one_sample.tiff", units="in", width=8, height=5, res=300)
one_sample <- ggplot(data = data, aes(x = n,  y = t)) +
  geom_point() +
  geom_smooth() +
  theme_classic() + grids(linetype = "dashed") +
  xlab("Triangulation Fineness Degree (n)") + ylab("Computing Time (s)") +
  ggtitle("Computing times for One-Sample SCC Estimation") + 
  theme(plot.title=element_text(hjust=1, vjust=0.5, face='bold')) +
  theme(legend.position=c(.9,.9))

dev.off()

###

# TWO-SAMPLE TIMES:

data <- read.csv2("two_sample_times.csv", sep = ";", dec = ",", header = T)

two_sample <- ggplot(data = data, aes(x = n,  y = t)) +
  geom_point() +
  geom_smooth() +
  theme_classic() + grids(linetype = "dashed") + ylim(0,200000) + xlim(5,25)+
  xlab("Triangulation Fineness Degree (n)") + ylab("Computing Time (s)") +
  ggtitle("Computing times for Two-Sample SCC Estimation") + 
  theme(plot.title=element_text(hjust=1, vjust=0.5, face='bold')) +
  theme(legend.position=c(.9,.9))

two_sample
ggplotly(two_sample)
  
## this is just for saving
tiff("two_sample.tiff", units="in", width=8, height=5, res=300)
ggplot(data = data, aes(x = n,  y = t, color = N)) +
  geom_point() +
  geom_smooth() +
  theme_classic() + grids(linetype = "dashed") + ylim(200,600) + xlim(50,500)+
  xlab("Nº Simulated Cases (n)") + ylab("Computing Time (s)") +
  ggtitle("Computing times for Two-Sample SCC Estimation") + 
  theme(plot.title=element_text(hjust=1, vjust=0.5, face='bold')) +
  theme(legend.position=c(.9,.9))

dev.off()



