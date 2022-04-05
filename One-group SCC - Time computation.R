####################################################################
# SCRIPT FOR COMPUTING TIMES OF SCC IN ONE-GROUP SETUP
####################################################################

setwd("~/MEGA/PhD/12. ICCSA Technical Publication/R Code")
load("contour.rda")


##########################################################
###            LIBRARIES
##########################################################

# install.packages(c("devtools", "remotes"))

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE) # Forces installation to finish even with warning messages
remotes::install_github("funstatpackages/ImageSCC")
remotes::install_github("funstatpackages/BPST")

library(devtools); library(remotes); library(ImageSCC); library(BPST); library(fields)


##########################################################
###               SIMULATION  OF  DATA                 ###
##########################################################

VT8 <- TriMesh(contour, n = 8)
VT10 <- TriMesh(contour, n = 10)
VT25 <- TriMesh(contour, n = 25)

# Triangulation Information  ## Data sets in package ‘ImageSCC’


TriPlot(Brain.V1, Brain.Tr1, col = 1, lwd = 1)
TriPlot(Brain.V2, Brain.Tr2, col = 1, lwd = 1)
TriPlot(Brain.V3, Brain.Tr3, col = 1, lwd = 1)


V.est1=as.matrix(Brain.V1)
Tr.est1=as.matrix(Brain.Tr1)
V.band1=as.matrix(Brain.V1)
Tr.band1=as.matrix(Brain.Tr1)

V.est2=as.matrix(Brain.V2)
Tr.est2=as.matrix(Brain.Tr2)
V.band2=as.matrix(Brain.V2)
Tr.band2=as.matrix(Brain.Tr2)

V.est3=as.matrix(Brain.V3)
Tr.est3=as.matrix(Brain.Tr3)
V.band3=as.matrix(Brain.V3)
Tr.band3=as.matrix(Brain.Tr3)

# Sample location information

n1=79; n2=79;
npix=n1*n2
u1=seq(0,1,length.out=n1)
v1=seq(0,1,length.out=n2)
uu=rep(u1,each=n2)
vv=rep(v1,times=n1)
uu.mtx=matrix(uu,n2,n1)
vv.mtx=matrix(vv,n2,n1)
Z=as.matrix(cbind(uu,vv))
ind.inside=inVT(V.est1,Tr.est1,Z[,1],Z[,2])$ind.inside


# Simulation parameters

d.est=5; d.band=2; r=1;
# n=50; # n=50/100/200
lam1=0.5; lam2=0.2;
lambda=10^{seq(-6,3,0.5)}
mu.func=3  # mu.func=1/2/3/4
alpha.grid=c(0.10,0.05,0.01)

# SIMULATION IN LOOP:

n <- rep(c(50, 100, 150, 200, 250, 300, 500), each=3) # three times for averaging

one_sample_times_N8 <- data.frame(n = integer(), t = numeric()) 
list <- list()

for (i in 1:length(n)){
 
    dat <-  data1g.image(n=n[i],Z,ind.inside,mu.func,noise.type='Func',lam1,lam2)
    Y <- dat$Y
    t <- system.time(out<-scc.image(Ya=Y,Z=Z,d.est=d.est,d.band=d.band,r=r,
                               V.est.a=V.est1,Tr.est.a=Tr.est1,
                               V.band.a=V.band1,Tr.band.a=Tr.band1,
                               penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,
                               adjust.sigma=TRUE))
    list[[i]] <- out
    datos <- data.frame(n = n[i], t = t[3])
    one_sample_times_N8 <- rbind(one_sample_times_N8, datos)
}

plot(out)

    plot(out,
        breaks=seq(from=0,to=2,length.out = 65),
        col=,
        xlab="Longitudinal (1-95)",
        ylab="Transversal (1-79)",
        sub="Control Group",
        col.sub="black",
        family ="serif")

########

one_sample_times_N15 <- data.frame(n = integer(), t = numeric()) 
list15 <- list()

for (i in 1:length(n)){
 
    dat <-  data1g.image(n=n[i],Z,ind.inside,mu.func,noise.type='Func',lam1,lam2)
    Y <- dat$Y
    t <- system.time(out<-scc.image(Ya=Y,Z=Z,d.est=d.est,d.band=d.band,r=r,
                               V.est.a=V.est2,Tr.est.a=Tr.est2,
                               V.band.a=V.band2,Tr.band.a=Tr.band2,
                               penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,
                               adjust.sigma=TRUE))
    list15[[i]] <- out
    datos <- data.frame(n = n[i], t = t[3])
    one_sample_times_N15 <- rbind(one_sample_times_N15, datos)
}

########

one_sample_times_N25 <- data.frame(n = integer(), t = numeric()) 
list25 <- list()

for (i in 1:length(n)){
 
    dat <-  data1g.image(n=n[i],Z,ind.inside,mu.func,noise.type='Func',lam1,lam2)
    Y <- dat$Y
    t <- system.time(out<-scc.image(Ya=Y,Z=Z,d.est=d.est,d.band=d.band,r=r,
                               V.est.a=V.est3,Tr.est.a=Tr.est3,
                               V.band.a=V.band3,Tr.band.a=Tr.band3,
                               penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,
                               adjust.sigma=TRUE))
    list25[[i]] <- out
    datos <- data.frame(n = n[i], t = t[3])
    one_sample_times_N25 <- rbind(one_sample_times_N25, datos)
}

#####
# SIMULATION FOR TWO-SAMPLE:
#####

two_sample_times_N8 <- data.frame(n = integer(), t = numeric()) 
list_two8 <- list()

for (i in 1:length(n)){
 
    dat <-  data1g.image(n=n[i],Z,ind.inside,mu.func,noise.type='Func',lam1,lam2)
    Y <- dat$Y
    t <- system.time(out<-scc.image(Ya=Y,Z=Z,d.est=d.est,d.band=d.band,r=r,
                               V.est.a=V.est1,Tr.est.a=Tr.est1,
                               V.band.a=V.band1,Tr.band.a=Tr.band1,
                               penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,
                               adjust.sigma=TRUE))
    list[[i]] <- out
    datos <- data.frame(n = n[i], t = t[3])
    one_sample_times_N8 <- rbind(one_sample_times_N8, datos)
}

plot(out)

    plot(out,
        breaks=seq(from=0,to=2,length.out = 65),
        col=,
        xlab="Longitudinal (1-95)",
        ylab="Transversal (1-79)",
        sub="Control Group",
        col.sub="black",
        family ="serif")        