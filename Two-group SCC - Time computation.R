####################################################################
# SCRIPT FOR COMPUTING TIMES OF SCC IN TWO-GROUP SETUP
####################################################################


if(!require(devtools)){install.packages('devtools');library(devtools)}
if(!require(BPST)){install_github("funstatpackages/BPST");library(BPST)}
if(!require(Triangulation)){install_github("funstatpackages/Triangulation");library(Triangulation)}
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE) # Forces installation to finish even with warning messages
remotes::install_github("funstatpackages/ImageSCC")

library(devtools); library(remotes); library(ImageSCC); library(BPST); library(fields)


# Triangulation Information
V.est.a=as.matrix(Brain.V1)
Tr.est.a=as.matrix(Brain.Tr1)
V.est.b=as.matrix(Brain.V1)
Tr.est.b=as.matrix(Brain.Tr1)
V.band.a=as.matrix(Brain.V1)
Tr.band.a=as.matrix(Brain.Tr1)
V.band.b=as.matrix(Brain.V1)
Tr.band.b=as.matrix(Brain.Tr1)

V.est.a2=as.matrix(Brain.V2)
Tr.est.a2=as.matrix(Brain.Tr2)
V.est.b2=as.matrix(Brain.V2)
Tr.est.b2=as.matrix(Brain.Tr2)
V.band.a2=as.matrix(Brain.V2)
Tr.band.a2=as.matrix(Brain.Tr2)
V.band.b2=as.matrix(Brain.V2)
Tr.band.b2=as.matrix(Brain.Tr2)

V.est.a3=as.matrix(Brain.V3)
Tr.est.a3=as.matrix(Brain.Tr3)
V.est.b3=as.matrix(Brain.V3)
Tr.est.b3=as.matrix(Brain.Tr3)
V.band.a3=as.matrix(Brain.V3)
Tr.band.a3=as.matrix(Brain.Tr3)
V.band.b3=as.matrix(Brain.V3)
Tr.band.b3=as.matrix(Brain.Tr3)

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
ind.inside=inVT(V.est.a,Tr.est.a,Z[,1],Z[,2])$ind.inside

# Simulation parameters
d.est=5; d.band=2; r=1;
n <- rep(c(50, 100, 150, 200, 250, 300, 500), each=3) # three times for averaging
lam1=0.5; lam2=0.2;
lambda=10^{seq(-6,3,0.5)}
mu1.func=3
delta=0.8   # delta=0,...,0.8
alpha.grid=c(0.10,0.05,0.01)

# Example: Run two sample SCC construction one time:

two_sample_times_N8 <- data.frame(n = integer(), t = numeric()) 
list_two_8 <- list()

for (i in 1:length(n)){
 
    dat=data2g.image(na=n[i],nb=n[i],Z,ind.inside,mu1.func,noise.type='Const',lam1,lam2,delta)
    Ya=dat$Ya
    Yb=dat$Yb
    t <- system.time(out_two <- scc.image(Ya=Ya,Yb=Yb,Z=Z,d.est=d.est,d.band=d.band,r=r,
                                          V.est.a=V.est.a,Tr.est.a=Tr.est.a,
                                          V.band.a=V.band.a,Tr.band.a=Tr.band.a,
                                          V.est.b=V.est.b,Tr.est.b=Tr.est.b,
                                          V.band.b=V.band.b,Tr.band.b=Tr.band.b,
              penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE))
    list_two_8[[i]] <- out_two
    datos <- data.frame(n = n[i], t = t[3])
    two_sample_times_N8 <- rbind(two_sample_times_N8, datos)
}

###

two_sample_times_N15 <- data.frame(n = integer(), t = numeric()) 
list_two_15 <- list()

for (i in 1:length(n)){
 
    dat=data2g.image(na=n[i],nb=n[i],Z,ind.inside,mu1.func,noise.type='Const',lam1,lam2,delta)
    Ya=dat$Ya
    Yb=dat$Yb
    t <- system.time(out_two <- scc.image(Ya=Ya,Yb=Yb,Z=Z,d.est=d.est,d.band=d.band,r=r,
                                          V.est.a=V.est.a2,Tr.est.a=Tr.est.a2,
                                          V.band.a=V.band.a2,Tr.band.a=Tr.band.a2,
                                          V.est.b=V.est.b2,Tr.est.b=Tr.est.b2,
                                          V.band.b=V.band.b2,Tr.band.b=Tr.band.b2,
              penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE))
    list_two_15[[i]] <- out_two
    datos <- data.frame(n = n[i], t = t[3])
    two_sample_times_N15 <- rbind(two_sample_times_N15, datos)
}

###

two_sample_times_N25 <- data.frame(n = integer(), t = numeric()) 
list_two_25 <- list()

for (i in 1:length(n)){
 
    dat=data2g.image(na=n[i],nb=n[i],Z,ind.inside,mu1.func,noise.type='Const',lam1,lam2,delta)
    Ya=dat$Ya
    Yb=dat$Yb
    t <- system.time(out_two <- scc.image(Ya=Ya,Yb=Yb,Z=Z,d.est=d.est,d.band=d.band,r=r,
                                          V.est.a=V.est.a3,Tr.est.a=Tr.est.a3,
                                          V.band.a=V.band.a3,Tr.band.a=Tr.band.a3,
                                          V.est.b=V.est.b3,Tr.est.b=Tr.est.b3,
                                          V.band.b=V.band.b3,Tr.band.b=Tr.band.b3,
              penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE))
    list_two_25[[i]] <- out_two
    datos <- data.frame(n = n[i], t = t[3])
    two_sample_times_N25 <- rbind(two_sample_times_N25, datos)
}

### This is for visualization of areas with values outside SCC

aa = out_two

Z.band <- matrix(aa$Z.band,ncol=2) # Posiciones
z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Posiciones por separado
n1 <- length(z1); n2 <- length(z2) # Longitud de dichas posiciones
scc <- matrix(NA,n1*n2,2) # Se crea la matriz donde irá el valor de SCC para cada posicion
ind.inside.band <- aa$ind.inside.cover # Solo las zonas cubiertas por la triangulacion
scc[ind.inside.band,] <- aa$scc[,,2] # se asigna el SCC a esas zonas (1,2,3 para los diferentes alpha)
scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # Limites: minimo de la inferior y máximo de la superior

scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
scc.l.mtx[scc.l.mtx<0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Lower SCC que son positivos se quedan
scc.u.mtx[scc.u.mtx>0]=NA # Los que cumplen lo esperable se ponen NA para no representarlos, solo los Upper SCC que son negativos se quedan


image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es positiva (cae encima del 0), o sea, que la imagen 1 es más fuerte que la imagen 2 en esas áreas
image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regiones en que la diferencia de medias es negativa (cae debajo del 0), o sea, que la imagen 2 es más fuerte que la imagen 1 en esas áreas

points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Puntos con diferencia de medias positiva (primera más fuerte)
points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Puntos con diferencia de medias negativa (segunda más fuerte)

points.P <- (points.P/79) # Re-escalation
points.N <- (points.N/79)


plot(out_two,
     breaks=c(0,100),
     col="turquoise",
     # breaks=seq(from=-1,to=6, length.out = 65),
     # xlab="Longitudinal (1-95)",
     # ylab="Transversal (1-79)",
     # sub="Difference between estimated mean functions: CNs - ADs",
     # col.sub="red",
     family ="serif")

    points(points.P,
           type="p",
           pch=".",
           col="red",
           cex=12)  
    
    points(points.N,
           type="p",
           pch=".",
           col="blue",
           cex=12) 
