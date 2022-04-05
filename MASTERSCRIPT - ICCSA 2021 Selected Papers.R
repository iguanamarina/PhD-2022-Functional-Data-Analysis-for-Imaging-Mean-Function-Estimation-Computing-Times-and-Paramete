############################## ################### ############## ### 
##
## Script name: MASTERSCRIPT - ICCSA 2021 Selected Papers
##
## Purpose of script: I've been requested with an extended article building upon 
## Arias-López, J. A., Cadarso-Suárez, C., & Aguiar-Fernández, P. (2021). Computational 
## Issues in the Application of Functional Data Analysis to Imaging Data (pp. 630–638). 
## https://doi.org/10.1007/978-3-030-86960-1_46
## Again, this is a methodological paper evaluating computing times for this FDA methodology
## but, this time, with a real practical case and with the service of a highly-specialized 
## computer. I aim to evaluate: optimal parameters for triangulations, computing times for
## one-sample and two-sample practical cases, along with optimal parameters to optimize 
## sensibility of this process given the TRUE pixels in a PET image which change between groups.
##
## Date Created: 2022-02-15
##
## Author: Juan A. Arias (M.Sc.)
## Email: juanantonio.arias.lopez@usc.es
## Webpage: https://messy-dataset.xyz
##
## Notes: Just run and have fun.
##       
##   
############################## ################### ############## ### 


####  
# PREAMBLE: ----
#### 


#* Set working directory: ----

setwd("~/MEGA/PhD/12. ICCSA Technical Publication/ICCSA Selected Papers '22/R Code")

#* Tune Options: ----
options(scipen = 6, digits = 4) # View outputs in non-scientific notation
memory.limit(30000000)     # This is needed on some PCs to increase memory allowance

#* Load up packages: ---- 

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE) # Forces installation to finish even with warning messages
install.packages(c("devtools","remotes","readr","imager","itsadug","ggplot2","contoureR","fields","plotly","gapminder", "ggplotly"))

remotes::install_github("funstatpackages/BPST")
remotes::install_github("funstatpackages/Triangulation")
remotes::install_github("funstatpackages/ImageSCC")

library(oro.nifti); library(memisc); library(remotes); library(ggplot2); library(contoureR); library(Triangulation); library(ggplot2);library(plotly);library(gapminder);library(readr);library(imager);library(itsadug);library(fields); library(BPST);library(Triangulation);library(ImageSCC); library(ggplotly)

library(devtools)

#* Load up functions: ----

load("complementary_materials/f.clean.RData") # Function to transform PET into a data.frame

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


####  
# PART 1: TRIANGULATION TIMES ----------
####  


#* Section 1: Get Contours ----------

param.z = 30 # Z=35 bc two-dimensional spaces
samplePET <- f.clean("complementary_materials/PETdataEXAMPLE.nii")
dfPET <- subset(samplePET, samplePET$z == param.z) 
dfPET <- dfPET[1:9919, "pet"] # Make sure to keep only 9919 and PET
dfPET <- as.matrix(dfPET)
dfPET <- t(dfPET) # Because we need the data in functional format
dfPET[is.nan(dfPET)] <- 0

# Our sample data is: 91x109

x <- rep(1:91, each = 109, length.out = 9919) 
y <- rep(1:109, length.out = 9919)
Z <- cbind(as.matrix(x),as.matrix(y))
rm(x); rm(y)
dat <- cbind(Z, t(dfPET)) # takes a sample to calculate boundaries
dat <- as.data.frame(dat)
sum(is.na(dat$pet)) # should be = 0
rownames(dat) <- NULL
df = getContourLines(dat[1:9919,], 
                     levels = c(0)) 

ggplot(df,aes(x,y,colour = z)) + geom_path() # Display of brain boundaries for Z=35
contour = df[, c("x", "y")]

save(contour, file = paste0("contour", as.numeric(param.z), ".rda"))  
  
 
#* Section 2: Get Sample of Triangulations ----------
 

# "N": An integer parameter controlling the fineness of the triangulation 
# and subsequent triangulation. As n increases the fineness increases and so does computing times.
# This is what I want to know: original authors of the methodology suggest N=8, but I suggest
# that a server can run much higher N values with relatively small changes in computing times. 

# TESTS:

load(paste0("~/MEGA/PhD/12. ICCSA Technical Publication/ICCSA Selected Papers '22/R Code/", "contour", as.numeric(param.z), ".rda"))

VT10 <- TriMesh(contour, n = 10)
VT25 <- TriMesh(contour, n = 25)
VT50 <- TriMesh(contour, n = 50) 

TriPlot(VT10$V, VT10$Tr, col = 1, lwd = 1)
TriPlot(VT25$V, VT25$Tr, col = 1, lwd = 1)
TriPlot(VT50$V, VT50$Tr, col = 1, lwd = 1)

# Save these plots for the article


#* Section 3: Get Triangulation Times ----------


n <- rep(2:100, each = 3) # three times for averaging

triangulation_times <- data.frame(n = integer(), t = numeric()) 

for (i in 1:length(n)) {
    t <- system.time(TriMesh(contour, n = n[i]))
    dat <- data.frame(n = n[i], t = t[3])
    triangulation_times <- rbind(triangulation_times, dat)
}

write.csv2(triangulation_times, file = "triangulation_times.csv", col.names = c(n, time))


#* RESULTS ----------


data_triangulation <- read.csv2("triangulation_times.csv", header = T, sep = ";")
data_triangulation <- data_triangulation[,2:3]

# For saving:
tiff("test.tiff", units = "in", width = 8, height = 5, res = 300)
ggplot(data_triangulation, aes(n, t)) +
       geom_point() + geom_smooth() +
       theme_classic() + grids(linetype = "dashed") +
       xlab("Fineness Degree (n)") + ylab("Computing Time (s)") +
       ggtitle("Computing times for Delaunay Triangulations") + 
       theme(plot.title = element_text(hjust = 1, vjust = 0.5, face = 'bold'))

dev.off()

# For interactive plot:
triang <- ggplot(data_triangulation, aes(n, t)) +
       geom_point() + geom_smooth() +
       theme_classic() + grids(linetype = "dashed") +
       xlab("Fineness Degree (n)") + ylab("Computing Time (s)") +
       ggtitle("Computing times for Delaunay Triangulations") + 
       theme(plot.title = element_text(hjust = 1, vjust = 0.5, face = 'bold'))
plotly::ggplotly(triang)


####
# PART 2: ONE-GROUP SCC ----
####

#* Section 1: Prepare the Data ----------

load("C:/Users/Juan A. Arias/Desktop/Preprocessed Data/Database_masked.RData")
Data <- Data[,-8]; attach(Data); rm(database)

# Names of PPT's for two groups (AD and CN) 

names_CN <- as.vector(Data[,1])[Data$group == "CN"]
names_CN <- names_CN[duplicated(names_CN) != TRUE]; names_CN

names_AD <- as.vector(Data[,1])[Data$group == "AD"]
names_AD <- names_AD[duplicated(names_AD) != TRUE]; names_AD

YcreatorCN <- function(i){
  Y <- subset(Data, Data$PPT == names_CN[i] & Data$z == param.z)
  Y <- Y[1:nrow(Z), 8] 
  Y <- as.matrix(Y)
  Y = t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

    SCC_matrix_CN <- matrix(ncol = nrow(Z), nrow = 0)
    for (i in 1:length(names_CN)) { 
      temp <- YcreatorCN(i)
      SCC_matrix_CN <- rbind(SCC_matrix_CN,temp)    
    }

YcreatorAD <- function(i){
  Y <- subset(Data, Data$PPT == names_AD[i] & Data$z == param.z) 
  Y <- Y[1:nrow(Z),8] 
  Y <- as.matrix(Y)
  Y = t(Y) 
  Y[is.nan(Y)] <- 0
  print(Y)
}

    SCC_matrix_AD <- matrix(ncol = nrow(Z),nrow = 0)
    for (i in 1:length(names_AD)) {
      temp <- YcreatorAD(i)
      SCC_matrix_AD <- rbind(SCC_matrix_AD,temp)    
    }


#* Section 2: Prepare Parameters ----         
    
VT8 <- TriMesh(contour, n = 8)        
VT15 <- TriMesh(contour, n = 15)
VT25 <- TriMesh(contour, n = 25)
VT35 <- TriMesh(contour, n = 35)
VT45 <- TriMesh(contour, n = 45)
VT55 <- TriMesh(contour, n = 55)
VT65 <- TriMesh(contour, n = 65)
VT75 <- TriMesh(contour, n = 75)
VT85 <- TriMesh(contour, n = 85)
VT95 <- TriMesh(contour, n = 95)

# Simulation parameters:

d.est = 5; d.band = 2; r = 1;
lambda = 10^{seq(-6,3,0.5)}
alpha.grid = c(0.10,0.05,0.01)
n <- rep(c(8, 15, 25, 35, 45, 55, 65, 75, 85, 95), each = 1) # three times for averaging
n <- c(5, 30)
VT5 <- TriMesh(contour, n = 5)
VT10 <- TriMesh(contour, n = 10)        
VT20 <- TriMesh(contour, n = 20)
VT30 <- TriMesh(contour, n = 30)

#* Section 3: Calculations ----

load("~/MEGA/PhD/12. ICCSA Technical Publication/ICCSA Selected Papers '22/R Code/SCC_CN.RData")
load("~/MEGA/PhD/12. ICCSA Technical Publication/ICCSA Selected Papers '22/R Code/SCC_roiAD_4.RData")

# Computing times:

one_sample_times <- data.frame(n = integer(), t = numeric()) 
list <- list()

for (i in 1:length(n)) {
    
    VT <- paste0("VT", as.numeric(n[i]))
    V.est <- as.matrix(get(VT)$V)
    Tr.est <- as.matrix(get(VT)$Tr)
    V.band <- as.matrix(get(VT)$V)
    Tr.band <- as.matrix(get(VT)$Tr)
    t <- system.time(out   <-  scc.image(Ya = SCC_matrix, Z = Z, d.est = d.est, d.band = d.band, r = r,
                               V.est.a = V.est, Tr.est.a = Tr.est,
                               V.band.a = V.band, Tr.band.a = Tr.band,
                               penalty = TRUE, lambda = lambda, alpha.grid = alpha.grid,
                               adjust.sigma = TRUE))
    list[[i]] <- out 
    datos <- data.frame(n = n[i], t = t[3])
    one_sample_times <- rbind(one_sample_times, datos)
    print(i)
}

write.csv2(one_sample_times, file = "one_sample_times.csv", col.names = c(n, time)) 

#* RESULTS ----

# Visualization of Images:
 
plot(list[[1]],   # Triangulation N=5
     breaks = seq(from = 0, to = 2, length.out = 65),
     col = ,
     xlab = "Longitudinal (1-109)",
     ylab = "Transversal (1-91)",
     sub = "Triangulation N=5",
     col.sub = "black",
     family = "serif")

plot(list[[30]],  # Triangulation N=95
     breaks = seq(from = 0, to = 2, length.out = 65),
     col = ,
     xlab = "Longitudinal (1-109)",
     ylab = "Transversal (1-91)",
     sub = "Triangulation N=95",
     col.sub = "black",
     family = "serif")


# Visualization of Computing Times:

one_sample <- ggplot(data = one_sample_times, aes(x = n,  y = t)) +
  geom_point() +
  geom_smooth() +
  theme_classic() + grids(linetype = "dashed") +
  xlab("Triangulation Degree of Fineness") + ylab("Computing Time (s)") +
  ggtitle("Computing times for One-Sample SCC Estimation") + 
  theme(plot.title = element_text(hjust = 1, vjust = 0.5, face = 'bold')) +
  theme(legend.position = c(.9,.9))

one_sample
ggplotly(one_sample)

# For saving:

ggplot(data_triangulation, aes(n, t)) +
       geom_point() + geom_smooth() +
       theme_classic() + grids(linetype = "dashed") +
       xlab("Fineness Degree (n)") + ylab("Computing Time (s)") +
       ggtitle("Computing times for One-Sample SCC") + 
       theme(plot.title = element_text(hjust = 1, vjust = 0.5, face = 'bold'))
tiff("one-sample.tiff", units = "in", width = 8, height = 5, res = 300)
dev.off()


####
# PART 3: TWO-GROUP SCC ----
####

#* Section 1: Prepare the Data ----

# Computing times:

two_sample_times <- data.frame(n = integer(), t = numeric()) 
list2 <- list()

for (i in 1:length(n)) {
    
    VT <- paste0("VT", as.numeric(n[i]))
    t <- system.time(out   <-  scc.image( Ya = SCC_matrix_AD, Yb = SCC_matrix_CN, Z = Z, d.est = d.est, d.band = d.band,
                                          r = r, V.est.a = get(VT)$V, Tr.est.a = get(VT)$Tr,
                                          V.band.a = get(VT)$V, Tr.band.a = get(VT)$Tr,
                                          V.est.b = get(VT)$V, Tr.est.b = get(VT)$Tr,
                                          V.band.b = get(VT)$V, Tr.band.b = get(VT)$Tr,
                                          penalty = TRUE, lambda = lambda, alpha.grid = alpha.grid,
                                          adjust.sigma = TRUE))
    list2[[i]] <- out
    datos <- data.frame(n = n[i], t = t[3])
    two_sample_times <- rbind(two_sample_times, datos)
    print(i)
}


#* RESULTS ----

# Visualization of Images:

load("complementary_materials/my_points.RData")

points <- my_points(list2[[1]], 2)

plot(list2[[1]],   # Triangulation N=5
     breaks = seq(from = 0, to = 2, length.out = 65),
     col = ,
     xlab = "Longitudinal (1-109)",
     ylab = "Transversal (1-91)",
     sub = "Triangulation N=5",
     col.sub = "black",
     family = "serif")

plot(out_two,
     breaks = c(0,100),
     col = "turquoise",
     # breaks=seq(from=-1,to=6, length.out = 65),
     # xlab="Longitudinal (1-95)",
     # ylab="Transversal (1-79)",
     # sub="Difference between estimated mean functions: CNs - ADs",
     # col.sub="red",
     family = "serif")

    points(points[1],
           type = "p",
           pch = ".",
           col = "red",
           cex = 12)  
    
    points(points[2],
           type = "p",
           pch = ".",
           col = "blue",
           cex = 12) 


# Visualization of Computing Times:

two_sample <- ggplot(two_sample_times = datos, aes(x = n,  y = t)) +
  geom_point() +
  geom_smooth() +
  theme_classic() + grids(linetype = "dashed") +
  xlab("Triangulation Degree of Fineness") + ylab("Computing Time (s)") +
  ggtitle("Computing times for Two-Sample SCC Estimation") + 
  theme(plot.title = element_text(hjust = 1, vjust = 0.5, face = 'bold')) +
  theme(legend.position = c(.9,.9))

two_sample
ggplotly(two_sample)

# For saving:

ggplot(two_sample_times, aes(n, t)) +
       geom_point() + geom_smooth() +
       theme_classic() + grids(linetype = "dashed") +
       xlab("Fineness Degree (n)") + ylab("Computing Time (s)") +
       ggtitle("Computing times for Two-Sample SCC") + 
       theme(plot.title = element_text(hjust = 1, vjust = 0.5, face = 'bold'))
tiff("one-sample.tiff", units = "in", width = 8, height = 5, res = 300)
dev.off()