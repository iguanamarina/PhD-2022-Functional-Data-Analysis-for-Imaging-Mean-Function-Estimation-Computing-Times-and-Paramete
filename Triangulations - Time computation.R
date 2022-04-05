######################################################
#  SCRIPT FOR TRIANGULATION TIME COMPUTATION
######################################################

# FIRST WE LOAD SAMPLE DATA FOR A PET STUDY TO FIND CONTOURS FOR DATA TO BE SIMULATED

##########################################################
###             f.clean  (z,x,y,pet)                   ###
##########################################################

library(oro.nifti); library(memisc)

## Set as working directory -> directory with sample NIfTI files

setwd("~/MEGA/PhD/12. ICCSA Technical Publication/R Code")
load("f.clean.RData")

#Example:
"samplePET" <- f.clean("PETdataEXAMPLE")
head(samplePET)

# We will need the data in Functional Format:

dfPET <- subset(samplePET, samplePET$z == 40) # Z=40 to work in two-dimensional spaces
dfPET <- dfPET[1:9919, "pet"] # Make sure to keep only 9919 and PET
dfPET <- as.matrix(dfPET)
dfPET <- t(dfPET) 
dfPET[is.nan(dfPET)] <- 0
View(dfPET)


##########################################################
###            LIBRARIES FOR TRIANGULATION             ###
##########################################################

install.packages(c("devtools", "remotes", "ggplot2", "contoureR"))

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE) # Forces installation to finish even with warning messages
remotes::install_github("funstatpackages/Triangulation")

library(devtools); library(remotes); library(ggplot2); library(contoureR); library(Triangulation)

# Z are the coordinates where data is measures in the given grid:
# our sample data is: 91x109

x <-rep(1:91, each=109, length.out = 9919) 
y <-rep(1:109, length.out = 9919)
Z <- cbind(as.matrix(x),as.matrix(y))
head(Z)

# Bondaries calculation:

dat <- cbind(Z, t(dfPET)) # takes a sample to calculate boundaries
dat <- as.data.frame(dat)
sum(is.na(dat$pet)) # should be = 0
rownames(dat) <- NULL

memory.size(max = TRUE) # Usually necessary
df = getContourLines(dat[1:9919,], 
                     levels=c(0)) # Place where 0 turns to other values (boundaries)

ggplot(df,aes(x,y,colour=z)) + geom_path() # Display of brain boundaries for Z=30

contour = df[, c("x", "y")] # These are the contours we need

save(contour, file = "contour.rda")

###################################
#  DELAUNAY TRIANGULATIONS PER SE:
###################################

# "N": An integer parameter controlling the fineness of the triangulation 
# and subsequent triangulation. As n increases the fineness increases. 

# TESTS:

VT10 <- TriMesh(contour, n = 10)
VT25 <- TriMesh(contour, n = 25)
VT50 <- TriMesh(contour, n = 50) 

TriPlot(VT10$V, VT10$Tr, col = 1, lwd = 1)
TriPlot(VT25$V, VT25$Tr, col = 1, lwd = 1)
TriPlot(VT50$V, VT50$Tr, col = 1, lwd = 1)


n <- rep(2:100, each=3) # three times for averaging

triangulation_times <- data.frame(n = integer(), t = numeric()) 

for (i in 1:length(n)){
    t <- system.time(TriMesh(contour, n = n[i]))
    dat <- data.frame(n = n[i], t = t[3])
    triangulation_times <- rbind(triangulation_times, dat)
}

View(triangulation_times)
write.csv2(triangulation_times, file = "triangulation_times.csv", col.names = c(n, time))
