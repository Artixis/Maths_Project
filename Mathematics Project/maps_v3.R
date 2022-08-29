require(rgdal)
require(raster)

require(gmt)
require(rgeos)
require(forecast)


# Weighting for kernel smoother
w=15


# GAussian kernel smooth weights
wt <- function(d,w) {
  exp(-0.5*(d/w)^2)
}


mig2 <- read.csv("./data/RIME_REGION_22112016094642536.csv")


ogrDrivers()

#au= readOGR("./data/1270055001_sa2_2016_aust_shape", "SA2_2016_AUST")
#au$SA2_MAIN16 <- as.numeric(levels(au$SA2_MAIN16))[au$SA2_MAIN16]
#syd <- crop(au,extent(150,152.5,-34.5,-32.5))
#writeOGR(syd,"./data/1270055001_sa2_2016_aust_shape","syd","ESRI Shapefile")


mll <- c(143.5,146.5,-39,-36.5)

#mel <- crop(au,extent(mll))
#writeOGR(mel,"./data/1270055001_sa2_2016_aust_shape","mel","ESRI Shapefile")

syd = readOGR("./data/1270055001_sa2_2016_aust_shape", "syd")
mel = readOGR("./data/1270055001_sa2_2016_aust_shape", "mel")


temp <- ecdf( abs(mig2$Value) )
mig2$col <- ceiling(temp(mig2$Value)*5)
mig2$col <- ifelse( mig2$col>0 , mig2$col+5 , -ceiling(temp(-mig2$Value)*5)+5 )






#### SYDNEY #########

png("./output/figs/Syd_Migration.png",width=595,height=842)
par(mfrow=c(3,3))
par(mar=c(0.1,0.1,4,0.1))
for(year in levels(factor(mig2$TIME_FY))) {
  
  mig2012 <- mig2[mig2$TIME_FY==year,] 
  plot(syd,col=cm.colors(10)
     [mig2012$col[match(syd$SA2_MAIN16,mig2012$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]]
     ,xlim=c(150,152.5),ylim=c(-34.5,-32.5),main=year,lwd=0.01)

  points(151.7817,-32.9283,cex=2,lwd=3)
  text(151.7817+0.05,-32.9283,cex=1,"Newcastle",adj=c(0,0))

  points(150.8931,-34.4278,cex=2,lwd=3)
  text(150.8931+0.05,-34.4278,cex=1,"Wollongong",adj=c(0,0))

  points(150.3119,-33.7125,cex=2,lwd=3)
  text(150.3119+0.05,-33.7125,cex=1,"Katoomba",adj=c(0,0))
 
  
  points(151.2093,-33.8688,cex=2,lwd=3)
  text(151.2093+0.05,-33.8688,cex=1,"Sydney",adj=c(0,0))
  
  
}
legend("bottomright",pch=22,pt.bg=cm.colors(2)[c(2,1)],c("positive","negative"),cex=2)

dev.off()


# Centroid of each Region
centroids <- getSpPPolygonsLabptSlots(syd)

# kernel smoothing
#close <- dd[i,] < c/1000
#weights <- exp(-0.5*(dd[!close,i]/w)^2)
#proresid2[i] <- sum( proresid[!close] * weights / sum(weights) )			

n=dim(centroids)[1]
dd <- matrix(nrow=n,ncol=n)
for(i in 1:n) {
  for(j in 1:n) {
    dd[i,j] <- geodist(centroids[i,2],centroids[i,1],centroids[j,2],centroids[j,1]) 
  }
}

# grid=> 150,152.5,-34.5,-32.5
stp = 0.02 # deg
lat_pts <- 1+(-32.5--34.5)/stp
lon_pts <-  1+(152.5-150)/stp
points <- data.frame(lon=rep(seq(150,152.5,stp),each=lat_pts)
                     , lat=rep(seq(-34.5,-32.5,stp),times=lon_pts),land=NA)
coordinates(points) <- c("lon","lat")
proj4string(points) <- proj4string(syd)

inside.park <- !is.na(over(points, as(syd, "SpatialPolygons")))
points@data$land <- inside.park

for(year in levels(factor(mig2$TIME_FY))) {
  mig <- mig2[mig2$TIME_FY==year,] 
  v <- mig$Value[match(syd@data$SA2_MAIN16,mig$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]
  npts <- dim(points)[1]
  for(i in 1:npts) {
    if(points@data$land[i]) {
      d <- sapply(1:n,function(j) geodist(centroids[j,2],centroids[j,1]
               ,points@coords[i,2],points@coords[i,1]) )
      wts <- wt(d,w)
      points@data[i,paste0("wm",year)] <- mean( v * wts / sum(wts) , na.rm=T)
    } else {
      points@data[i,paste0("wm",year)] <- NA
    }
  }
}

mm <- max(abs(range(points@data[,paste0("wm",2007:2015)],na.rm=T)))


        
png("./output/figs/Syd_Migration_smoothed.png",width=842,height=842)
par(mfrow=c(3,3))
par(mar=c(0.4,0.4,4,0.4))
for(year in levels(factor(mig2$TIME_FY))) {
  col = floor(100*points@data[,paste0("wm",year)]/(mm+1e-6))+101 # a number between 1 and 200
  
  plot(points,col=cm.colors(200)[col]
       ,xlim=c(150,152.5),ylim=c(-34.5,-32.5),main=year,pch=15,cex.main=2)
  box()
  points(151.7817,-32.9283,cex=2,lwd=3)
  text(151.7817+0.1,-32.9283,cex=1.5,"Newcastle",adj=c(0,0))
  
  points(150.8931,-34.4278,cex=2,lwd=3)
  text(150.8931+0.1,-34.4278,cex=1.5,"Wollongong",adj=c(0,0))
  
  points(150.3119,-33.7125,cex=2,lwd=3)
  text(150.3119+0.1,-33.7125,cex=1.5,"Katoomba",adj=c(0,0))
  
  
  points(151.2093,-33.8688,cex=2,lwd=3)
  text(151.2093+0.1,-33.8688,cex=1.6,"Sydney",adj=c(0,0))
  
  
}
legend("bottomright",pch=22,pt.bg=cm.colors(2)[c(2,1)],c("positive","negative"),cex=2)

dev.off()


wg=30

png("./output/figs/Temporal_syd_migration.png")
years <- levels(factor(mig2$TIME_FY))
ny <- paste0("wm",years)

#sydney
d <- sapply(1:n,function(j) geodist(centroids[j,2],centroids[j,1]
                                    ,-33.8688,151.2093) )
wts <- wt(d,wg)
mig60 <- data.frame(year=as.numeric(levels(factor(mig2$TIME_FY))),syd=NA)
for(year in levels(factor(mig2$TIME_FY))) {
  mig <- mig2[mig2$TIME_FY==year,] 
  v <- mig$Value[match(syd@data$SA2_MAIN16,mig$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]
  mig60$syd[mig60$year==year] <- mean( v * wts / sum(wts) , na.rm=T)
}

plot(mig60$year,mig60$syd,type="b",xlim=c(2007,2020),ylim=c(-0.3,0.3)
     ,pch=25,bg="grey",ylab="Average Net Migration",xlab="year")
fit <- lm(syd~year,data=mig60)
abline(fit,lty=2,col="grey")


#Wollongong
d <- sapply(1:n,function(j) geodist(centroids[j,2],centroids[j,1]
,-34.4278,150.8931) )
wts <- wt(d,wg)
for(year in levels(factor(mig2$TIME_FY))) {
  mig <- mig2[mig2$TIME_FY==year,] 
  v <- mig$Value[match(syd@data$SA2_MAIN16,mig$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]
  mig60$wol[mig60$year==year] <- mean( v * wts / sum(wts) , na.rm=T)
}

lines(mig60$year,mig60$wol,type="b",pch=21,bg="red")
fit <- lm(wol~year,data=mig60)
abline(fit,lty=2,col="red")

#Newcastle
d <- sapply(1:n,function(j) geodist(centroids[j,2],centroids[j,1]
                                    ,-32.9283,151.7817) )
wts <- wt(d,wg)
for(year in levels(factor(mig2$TIME_FY))) {
  mig <- mig2[mig2$TIME_FY==year,] 
  v <- mig$Value[match(syd@data$SA2_MAIN16,mig$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]
  mig60$new[mig60$year==year] <- mean( v * wts / sum(wts) , na.rm=T)
}

lines(mig60$year,mig60$new,type="b",pch=22,bg="blue")
fit <- lm(new~year,data=mig60)
abline(fit,lty=2,col="blue")

#abline(lm(diff~as.numeric(years)),lty=2,col="blue")

#Katoomba
d <- sapply(1:n,function(j) geodist(centroids[j,2],centroids[j,1]
                                    ,-33.7125,150.3119) )
wts <- wt(d,wg)
for(year in levels(factor(mig2$TIME_FY))) {
  mig <- mig2[mig2$TIME_FY==year,] 
  v <- mig$Value[match(syd@data$SA2_MAIN16,mig$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]
  mig60$kat[mig60$year==year] <- mean( v * wts / sum(wts) , na.rm=T)
}

lines(mig60$year,mig60$kat,type="b",pch=24,bg="green")
fit <- lm(kat~year,data=mig60)
abline(fit,lty=2,col="green")

abline(h=0,lty="dotted")
legend("bottomright",legend=c("Wollongong","Newcastle","Katoomba","Sydney")
       ,pch=c(21,22,24,25),pt.bg=c("red","blue","green","grey"))
dev.off()

points_syd <- points










###MELBOURNE###

tiff("./output/figs/MEL_Migration.tiff",width=595,height=842,compression="none")
par(mfrow=c(3,3))
par(mar=c(0.1,0.1,4,0.1))
for(year in levels(factor(mig2$TIME_FY))) {
  
  mig2012 <- mig2[mig2$TIME_FY==year,] 
  plot(mel,col=cm.colors(10)
       [mig2012$col[match(mel$SA2_MAIN16,mig2012$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]]
       ,xlim=mll[1:2],ylim=mll[3:4],main=year,lwd=0.01,cex.main=1.5)
  
  points(145.9298,-38.1538,cex=2,lwd=3)
  text(145.9298+0.05,-38.1538,cex=1.3,"Warragul",adj=c(0,1))
  
  points(144.3617,-38.1499,cex=2,lwd=3)
  text(144.3617+0.05,-38.1499,cex=1.3,"Geelong",adj=c(0,0))
  
  points(143.8503,-37.5622,cex=2,lwd=3)
  text(143.8503+0.05,-37.5622,cex=1.3,"Ballarat",adj=c(0,0))
  
  points(144.2794,-36.7570,cex=2,lwd=3)
  text(144.2794+0.05,-36.7570,cex=1.3,"Bendigo",adj=c(0,0))
  
  points(144.9631,-37.8136,cex=2,lwd=3)
  text(144.9631+0.05,-37.8136,cex=1.3,"Melbourne",adj=c(0,0))
  
  
}
legend("bottomright",pch=22,pt.bg=cm.colors(2)[c(2,1)],c("positive","negative"),cex=1.5)

dev.off()


# Centroid of each Region
centroids <- getSpPPolygonsLabptSlots(mel)

# kernel smoothing
#close <- dd[i,] < c/1000
#weights <- exp(-0.5*(dd[!close,i]/w)^2)
#proresid2[i] <- sum( proresid[!close] * weights / sum(weights) )			

n=dim(centroids)[1]
dd <- matrix(nrow=n,ncol=n)
for(i in 1:n) {
  for(j in 1:n) {
    dd[i,j] <- geodist(centroids[i,2],centroids[i,1],centroids[j,2],centroids[j,1]) 
  }
}

#> mll
#[1] 143.5 146.5 -39.0 -36.5
stp = 0.04 # deg
lat_pts <- 1+(mll[4]-mll[3])/stp
lon_pts <-  1+(mll[2]-mll[1])/stp
points <- data.frame(lon=rep(seq(mll[1],mll[2],stp),each=lat_pts)
                     , lat=rep(seq(mll[3],mll[4],stp),times=lon_pts),land=NA)
coordinates(points) <- c("lon","lat")
proj4string(points) <- proj4string(mel)

inside.park <- !is.na(over(points, as(mel, "SpatialPolygons")))
points@data$land <- inside.park

wt <- function(d,w) {
  exp(-0.5*(d/w)^2)
}


for(year in levels(factor(mig2$TIME_FY))) {
  mig <- mig2[mig2$TIME_FY==year,] 
  v <- mig$Value[match(mel@data$SA2_MAIN16,mig$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]
  mel@data[,paste0("mig",year)] <- v
  npts <- dim(points)[1]
  for(i in 1:npts) {
    if(points@data$land[i]) {
      d <- sapply(1:n,function(j) geodist(centroids[j,2],centroids[j,1]
                                          ,points@coords[i,2],points@coords[i,1]) )
      wts <- wt(d,w)
      points@data[i,paste0("wm",year)] <- mean( v * wts / sum(wts) , na.rm=T)
    } else {
      points@data[i,paste0("wm",year)] <- NA
    }
  }
}

mm <- max(abs(range(points@data[,paste0("wm",2007:2015)],na.rm=T)))



png("./output/figs/MEL_Migration_smoothed.png",width=842,height=842)
par(mfrow=c(3,3))
par(mar=c(0.4,0.4,4,0.4))
for(year in levels(factor(mig2$TIME_FY))) {
  col = floor(100*points@data[,paste0("wm",year)]/(mm+1e-6))+101 # a number between 1 and 200
  
  plot(points,col=cm.colors(200)[col]
       ,xlim=mll[1:2],ylim=mll[3:4],main=year,pch=15,cex.main=2)
  box()
  
  points(145.9298,-38.1538,cex=2,lwd=3)
  text(145.9298+0.1,-38.1538,"Warragul",adj=c(0,1),cex=1.5)
  
  points(144.3617,-38.1499,cex=2,lwd=3)
  text(144.3617+0.1,-38.1499,"Geelong",adj=c(0,0),cex=1.5)
  
  points(143.8503,-37.5622,cex=2,lwd=3)
  text(143.8503+0.1,-37.5622,"Ballarat",adj=c(0,0),cex=1.5)
  
  points(144.2794,-36.7570,cex=2,lwd=3)
  text(144.2794+0.1,-36.7570,"Bendigo",adj=c(0,0),cex=1.5)
  
  points(144.9631,-37.8136,cex=2,lwd=3)
  text(144.9631+0.1,-37.8136,"Melbourne",adj=c(0,0),cex=1.6)
  
  
}
legend("bottomright",pch=22,pt.bg=cm.colors(2)[c(2,1)],c("positive","negative"),cex=2)

dev.off()



wg = 30
png("./output/figs/Temporal_mel_migration.png")
years <- levels(factor(mig2$TIME_FY))
ny <- paste0("wm",years)

#melbourne
d <- sapply(1:n,function(j) geodist(centroids[j,2],centroids[j,1]
                                    ,-37.8136,144.9631) )
wts <- wt(d,wg)
mig60 <- data.frame(year=as.numeric(levels(factor(mig2$TIME_FY))),mel=NA)
for(year in levels(factor(mig2$TIME_FY))) {
  mig <- mig2[mig2$TIME_FY==year,] 
  v <- mig$Value[match(mel@data$SA2_MAIN16,mig$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]
  mig60$mel[mig60$year==year] <- mean( v * wts / sum(wts) , na.rm=T)
}
plot(mig60$year,mig60$mel,type="b",xlim=c(2007,2020),ylim=c(-0.5,0.5),pch=25,bg="grey"
     ,ylab="Average Net Migration",xlab="year")
fit <- lm(mel~year,data=mig60)
abline(fit,lty=2,col="grey")


#Geelong
d=sapply(1:n,function(j) geodist(-38.1499,144.3617
                                    ,centroids[j,2],centroids[j,1]) )
wts <- wt(d,wg)
mig60 <- data.frame(year=as.numeric(levels(factor(mig2$TIME_FY))),mel=NA)
for(year in levels(factor(mig2$TIME_FY))) {
  mig <- mig2[mig2$TIME_FY==year,] 
  v <- mig$Value[match(mel@data$SA2_MAIN16,mig$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]
  mig60$mel[mig60$year==year] <- mean( v * wts / sum(wts) , na.rm=T)
}
lines(mig60$year,mig60$mel,type="b",pch=21,bg="lightskyblue")
fit <- lm(mel~year,data=mig60)
abline(fit,lty=2,col="lightskyblue")

#Bendigo
d=sapply(1:n,function(j) geodist(-36.7570,144.2794
                                    ,centroids[j,2],centroids[j,1]) )
wts <- wt(d,wg)
mig60 <- data.frame(year=as.numeric(levels(factor(mig2$TIME_FY))),mel=NA)
for(year in levels(factor(mig2$TIME_FY))) {
  mig <- mig2[mig2$TIME_FY==year,] 
  v <- mig$Value[match(mel@data$SA2_MAIN16,mig$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]
  mig60$mel[mig60$year==year] <- mean( v * wts / sum(wts) , na.rm=T)
}
lines(mig60$year,mig60$mel,type="b",pch=22,bg="maroon")
fit <- lm(mel~year,data=mig60)
abline(fit,lty=2,col="maroon")

#Ballarat
d=sapply(1:n,function(j) geodist(-37.5622,143.8503
                                 ,centroids[j,2],centroids[j,1]) )
wts <- wt(d,wg)
mig60 <- data.frame(year=as.numeric(levels(factor(mig2$TIME_FY))),mel=NA)
for(year in levels(factor(mig2$TIME_FY))) {
  mig <- mig2[mig2$TIME_FY==year,] 
  v <- mig$Value[match(mel@data$SA2_MAIN16,mig$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]
  mig60$mel[mig60$year==year] <- mean( v * wts / sum(wts) , na.rm=T)
}
lines(mig60$year,mig60$mel,type="b",pch=24,bg="green")
fit <- lm(mel~year,data=mig60)
abline(fit,lty=2,col="green")


#Warragul
d=sapply(1:n,function(j) geodist(-38.1538,145.9298
                                 ,centroids[j,2],centroids[j,1]) )
wts <- wt(d,wg)
mig60 <- data.frame(year=as.numeric(levels(factor(mig2$TIME_FY))),mel=NA)
for(year in levels(factor(mig2$TIME_FY))) {
  mig <- mig2[mig2$TIME_FY==year,] 
  v <- mig$Value[match(mel@data$SA2_MAIN16,mig$ASGS_2011_STATE_GCCSA_SA4_SA3_SA2)]
  mig60$mel[mig60$year==year] <- mean( v * wts / sum(wts) , na.rm=T)
}
lines(mig60$year,mig60$mel,type="b",pch=23,bg="yellow")
fit <- lm(mel~year,data=mig60)
abline(fit,lty=2,col="yellow")

abline(h=0,lty="dotted")
legend("bottomright",legend=c("Geelong","Bendigo","Ballarat","Warragul","Melbourne")
       ,pch=c(21,22,24),pt.bg=c("lightskyblue","maroon","green","yellow","grey"))
dev.off()



