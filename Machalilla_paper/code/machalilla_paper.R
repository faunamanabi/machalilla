

library(lubridate)
# library(dplyr)
library(maptools)
# Load 'rgdal' package, which is used to read/write shapefiles and rasters
library(rgdal)
source("C:/Users/CaracterizaciónD/Documents/GitHub/pacoche/code/TEAM_code.R")
# source("C:/Users/Diego/Documents/CodigoR/ULEAM/Infor_Caract/code/calendar.R")
load(file = "data/machalilla_fixed.RData")



#machalilla.fixed<-read.csv("C:/Users/Diego/Documents/CodigoR/Machalilla_paper/data/machalilla_fixed.csv") 

#machalilla.fixed$photo_date2<-as.Date(as.character(machalilla.fixed$camera_trap_start_date),  "%Y-%m-%d")
#machalilla.fixed$Sampling.Period<-as.numeric(machalilla.fixed$Sampling.Period)
#machalilla.fixed$binomial<-as.character(machalilla.fixed$binomial)
#machalilla.fixed$year<-as.numeric(machalilla.fixed$year)
#machalilla.fixed$month<-as.numeric(machalilla.fixed$month)

#machalilla.fixed$camera_trap_start_date<-as.Date(as.character(machalilla.fixed$camera_trap_start_date),  "%Y-%m-%d")
#machalilla.fixed$camera_trap_end_date<-as.Date(as.character(machalilla.fixed$camera_trap_end_date),  "%Y-%m-%d")


# library(dplyr)
mat.per.sp<-f.matrix.creator2(data = machalilla.fixed,year = 2014)
sp.names<-names(mat.per.sp) # species names

# counting how many (total) records per species by all days
cont.per.sp<-data.frame(row.names = sp.names)
row.per.sp<-as.data.frame(matrix(nrow = length(sp.names), ncol=c(60)))
col.per.sp<-as.data.frame(matrix(nrow = length(sp.names), ncol=c(161)))
rownames(row.per.sp)<-sp.names
rownames(col.per.sp)<-sp.names

for (i in 1:length(mat.per.sp)){
  cont.per.sp[i,1]<-sum(apply(as.data.frame(mat.per.sp [[i]]),FUN=sum,na.rm=T, MARGIN = 1))
  
  row.per.sp[i,]<-apply(mat.per.sp[[i]],1, function(x) sum(x, na.rm=T))
  # row.per.sp[i,which(row.per.sp[i,]>0)]<-1 # convert to presence absence  1 and 0
  
  col.per.sp[i,]<-apply(mat.per.sp[[i]],2, function(x) sum(x, na.rm=T))
  # col.per.sp[i,which(col.per.sp[i,]>0)]<-1 # convert to presence absence  1 and 0
  
}


#### se full mat
row.per.sp

# colnames(row.per.sp)<-rep
# cont.per.sp$especie<-rownames(cont.per.sp)
# colnames(cont.per.sp)<-c("Numero_de_fotos","especie")

##################################################
### Just wild mammals not including domesticated
##################################################

# Chk names in sp.names Myotis myotis should be Artibeus sp
row.per.sp<-row.per.sp[-37,] # elimina Pheucticus chrysogaster
row.per.sp<-row.per.sp[-36,] # elimina Cathartes aura
row.per.sp<-row.per.sp[-35,] # elimina Gallus gallus
row.per.sp<-row.per.sp[-34,] # elimina Zenaida auriculata
row.per.sp<-row.per.sp[-33,] # elimina Momotus momota

row.per.sp<-row.per.sp[-30,] # elimina Heliomaster longirostris 
row.per.sp<-row.per.sp[-29,] # elimina Pipistrellus pipistrellus 
row.per.sp<-row.per.sp[-26,] # elimina Ortalis vetula 
row.per.sp<-row.per.sp[-25,] # elimina Homo sapiens
row.per.sp<-row.per.sp[-21,] # elimina Tinamus major 

row.per.sp<-row.per.sp[-16,] # elimina Equus africanus
row.per.sp<-row.per.sp[-15,] # elimina Sus scrofa
row.per.sp<-row.per.sp[-14,] # elimina Rattus rattus
row.per.sp<-row.per.sp[-12,] # elimina Equus ferus
row.per.sp<-row.per.sp[-10,] # elimina Bos primigenius

row.per.sp<-row.per.sp[-07,] # elimina Buteogallus urubitinga
row.per.sp<-row.per.sp[-04,] # elimina Leptotila verreauxi

row.per.sp<-row.per.sp[-03,] # elimina Capra aegagrus
row.per.sp<-row.per.sp[-02,] # elimina Canis lupus

row.per.sp<-row.per.sp[-01,] # elimina 1a sp
# row.per.sp<-row.per.sp[-10,] # elimina Leptotila verreauxi
# row.per.sp<-row.per.sp[-08,] # elimina Tinamus major
# row.per.sp<-row.per.sp[-01,] # elimina 1era especie vacia


######### CHK names
rownames(row.per.sp)


# Numero de especies observadas
length(row.per.sp[,1])

############################################################
## first Table 
############################################################
library(knitr)
presence<-row.per.sp
presence[presence > 0]<-1 # replace to 1
naiveoccu<-apply(X = presence,FUN = sum, MARGIN = 1 ) / 55 #divided number of sites
events<-apply(X = row.per.sp ,FUN = sum, MARGIN = 1 )
RAI<-(events/(37.8 * 55)) * 100 #37.8 is average sampling effort, in days per camera
Table1<-cbind(events, RAI, naiveoccu)
kable(Table1, format = "rst")

############################################################
## Distribucion posterior de la riqueza de especies
############################################################

# Riqueza de especies y acumulación, modelando la ocurrencia y la detectabilidad. 
# Este análisis sigue el método de Dorazio et al. (2006).

source("C:/Users/CaracterizaciónD/Documents/GitHub/pacoche/code/MultiSpeciesSiteOcc.R")

X1 = as.matrix(row.per.sp) # col.per.sp por dias y row.per.sp por sitios (camaras)
nrepls = 130 #dias 
especies = MultiSpeciesSiteOcc(nrepls, X1)

# summary(especies$fit$sims.matrix)

alpha.post = especies$fit$sims.matrix[,"alpha"]
sigmaU.post = especies$fit$sims.matrix[,"sigma.u"]
N.post = especies$fit$sims.matrix[,"N"]

nsites = 60 
cum_sp<-CumNumSpeciesPresent(nsites, alpha.post, sigmaU.post, N.post)

#histogram of posteriors
hist(especies$fit$sims.matrix[,"N"],breaks = c(16:30), xlab="Number of mammal species", ylab="Relative frecuency", main="")
abline(v=length(row.per.sp[,1]),col="blue", lty = 2) # -> lines.histogram(*) observadas
abline(v=median(especies$fit$sims.matrix[,"N"]),col="red", lty = 2) # -> esperadas por detectabilidad


mean(especies$fit$sims.matrix[,"N"])
median(especies$fit$sims.matrix[,"N"])

#####################
#### compare to vegan
#####################
library(vegan)
sac <- specaccum(t(row.per.sp))
plot(sac)

sp1 <- specaccum(t(row.per.sp), method = "random")
sp2 <- specaccum(t(row.per.sp),  method = "rarefaction")
sp3 <- specaccum(t(row.per.sp),  method = "collector")

summary(sp1)
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sp1, col="yellow", add=TRUE, pch="+")

plot(sp2, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue") # rarefaction
#plot(sp2)
# boxplot(sp1, col="yellow", add=TRUE, pch="+")

plot(sp3, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
#plot(sp2)

H <- diversity(t(row.per.sp))
simp <- diversity(t(row.per.sp), "simpson")
invsimp <- diversity(t(row.per.sp), "inv")
r.2 <- rarefy(t(row.per.sp), 2)
alpha <- fisher.alpha(t(row.per.sp))
pairs(cbind(H, simp, invsimp, r.2, alpha), pch="+", col="blue")

## Species richness (S) and Pielou's evenness (J):
S <- specnumber(t(row.per.sp)) ## rowSums(BCI > 0) does the same...
J <- H/log(S)

rarecurve(row.per.sp)#,step = 20, sample = raremax)
## Rarefaction
(raremax <- max(rowSums(row.per.sp)))
Srare <- rarefy(t(row.per.sp), raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(row.per.sp)#, step = 2, sample = raremax, col = "blue", cex = 0.6)



# filter mat.per.sp to sp of wild mammals
mammal.per.sp<- mat.per.sp[-c(37,36,35,34,33,30,29,26,25,21,16,15,14,12,10,7,4,3,2,1)]
full.mammal<-ldply(mammal.per.sp, data.frame)

sp.abund.count<-col.per.sp[-c(37,36,35,34,33,30,29,26,25,21,16,15,14,12,10,7,4,3,2,1),]

# 
# # convert list of species to unmarked 
# # fullmat2<-f.convert.to.unmarked(mammal.per.sp[2])
# spcies<-as.data.frame( full.mammal[,1])
# colnames(spcies)<-"sp"
# obs<-as.data.frame(full.mammal[,2:162])
# 
# #  covariates of detection and occupancy in that order.
# fullmat<-unmarkedFrameOccu(y=obs,siteCovs=spcies)
# 
# fm0 <- occu(~ 1 ~ 1, fullmat) 
# fm1 <- occu(~ 1 ~ sp, fullmat)
# fm2 <- occu(~ sp ~ 1, fullmat)
# fm3 <- occu(~ sp ~ sp, fullmat)
# 
# 
# models <- fitList(
#   'p(.)psi(.)' = fm0,
#   'p(.)psi(sp)' = fm1,
#   'p(sp)psi(.)' = fm2,
#   'p(sp)psi(sp)' = fm3)
# 
# ms <- modSel(models)
# ############### someting weird in number of parameters compared to Rovero

obs2<-row.per.sp #sp.abund.count
spcies2<-as.data.frame(rownames(row.per.sp)) #sp.abund.count
colnames(spcies2)<-"sp"


spcies2$guild<-c("Insectivore", "Carnivore", "Hervibivore",
                 "Carnivore",  "Hervibivore", "Hervibivore",
                 "Carnivore", "Carnivore", "Hervibivore",
                 "Onmivore", "Insectivore", "Hervibivore",
                 "Carnivore", "Hervibivore", "Carnivore",
                 "Onmivore", "Insectivore")
spcies2$mass<-c(4209, 3910, 55508, 3249, 949, 8000, 11900,
                6750, 2674, 21266, 4203, 433, 6950, 22800,
                4030, 1091, 42)
habitat<-c("dry","humid","dry","humid","dry",
           "dry","dry","dry","dry","dry",
           "dry", "humid", "humid", "transition", "humid",
           "transition", "humid",  "transition", "humid",
           "humid", "dry", "dry", "dry", "dry", "humid",
           "dry", "humid", "dry", "dry", "dry", "dry", "humid",
           "dry", "transition", "dry", "transition", "dry", "dry",
           "dry", "dry", "dry", "dry", "dry", "dry", "transition",
           "dry", "transition", "dry", "dry",  "dry", "dry", "transition",
           "dry", "transition",  "dry", "dry", "transition", "transition",
           "humid", "dry")

habitat2<-matrix(rep(habitat, 17), nrow = 17, ncol = 60, byrow = TRUE)

library(unmarked)
fullmat2<-unmarkedFrameOccu(y=obs2,siteCovs=spcies2)

#covariates of detection and occupancy in that order.
fm0 <- occu(~ 1 ~ 1, fullmat2) 
fm1 <- occu(~ 1 ~ guild, fullmat2)
fm2 <- occu(~ guild ~ 1, fullmat2)
fm3 <- occu(~ guild ~ guild, fullmat2)
fm4 <- occu(~ mass~ 1, fullmat2)
fm5 <- occu(~ mass ~ guild, fullmat2)
fm6 <- occu(~ mass + guild ~ 1, fullmat2)
fm7 <- occu(~ sp ~ guild, fullmat2)

models1 <- fitList(
  'p(.)psi(.)' = fm0,
  'p(.)psi(guild)' = fm1,
  'p(guild)psi(.)' = fm2,
  'p(guild)psi(guild)' = fm3,
  'p(mass)psi(.)' = fm4,
  'p(mass)psi(guild)' = fm5,
  'p(mass, guild)psi(.)' = fm6,
  'p(sp)psi(guild)' = fm7)

ms1 <- modSel(models1)


##############################

obs3<-t(row.per.sp) #sp.abund.count
species3<-(rownames(row.per.sp)) #sp.abund.count
# (spcies2)<-"sp"
species.obs<-matrix(rep(species3,60), # the data elements 
     nrow=60,              # number of rows 
     ncol=17,              # number of columns 
     byrow = TRUE)        # fill matrix by rows


guild<-c("Insectivore", "Carnivore", "Hervibivore",
                 "Carnivore",  "Hervibivore", "Hervibivore",
                 "Carnivore", "Carnivore", "Hervibivore",
                 "Onmivore", "Insectivore", "Hervibivore",
                 "Carnivore", "Hervibivore", "Carnivore",
                 "Onmivore", "Hervibivore")
guild.obs<-matrix(rep(guild,60), # the data elements 
                    nrow=60,              # number of rows 
                    ncol=17,              # number of columns 
                    byrow = TRUE)

mass<-c(4209, 3910, 55508, 3249, 949, 8000, 11900,
                6750, 2674, 21266, 4203, 433, 6950, 22800,
                4030, 1091, 42)
mass.obs<-matrix(rep(mass,60), # the data elements 
                  nrow=60,              # number of rows 
                  ncol=17,              # number of columns 
                  byrow = TRUE)

habitat3<-c("dry","humid","dry","humid","dry",
           "dry","dry","dry","dry","dry",
           "dry", "humid", "humid", "transition", "humid",
           "transition", "humid",  "transition", "humid",
           "humid", "dry", "dry", "dry", "dry", "humid",
           "dry", "humid", "dry", "dry", "dry", "dry", "humid",
           "dry", "transition", "dry", "transition", "dry", "dry",
           "dry", "dry", "dry", "dry", "dry", "dry", "transition",
           "dry", "transition", "dry", "dry",  "dry", "dry", "transition",
           "dry", "transition",  "dry", "dry", "transition", "transition",
           "humid", "dry")


obs.covs <- list(
  sp=species.obs,
  guild= guild.obs,
  mass= mass.obs)


fullmat3<-unmarkedFrameOccu(y=obs3, 
                            siteCovs=as.data.frame(habitat3), 
                            obsCovs=obs.covs)

#covariates of detection and occupancy in that order.
fm0 <- occu(~ 1 ~ 1, fullmat3) 
fm1 <- occu(~ guild ~ 1, fullmat3)
fm2 <- occu(~ mass ~ 1, fullmat3)
fm3 <- occu(~ guild + mass ~ 1, fullmat3)
fm4 <- occu(~ 1~ habitat3, fullmat3)
fm5 <- occu(~ guild ~ habitat3, fullmat3)
fm6 <- occu(~ mass ~ habitat3, fullmat3)
fm7 <- occu(~ guild + mass ~ habitat3, fullmat3)
fm8 <- occu(~ habitat3 ~ habitat3, fullmat3)
fm9 <- occu(~ sp ~ habitat3, fullmat3)



models2 <- fitList(
  'p(.)psi(.)' = fm0,
  'p(guild)psi(.)' = fm1,
  'p(mass)psi(.)' = fm2,
  'p(guild, mass)psi(.)' = fm3,
  'p(.)psi(habitat)' = fm4,
  'p(guild)psi(habitat)' = fm5,
  'p(mass)psi(habitat)' = fm6,
  'p(guild, mass)psi(habitat)' = fm7,
  'p(habitat)psi(habitat)' = fm8,
  'p(sp)psi(habitat)' = fm9)  

ms2 <- modSel(models2)


###############################################
##### Mammal tree for pyloenetic diversity
##### 
##### tree from: 
##### Fritz, S. A., Bininda-Emonds, O. R. P. & Purvis, A. 
##### 2009 Geographical variation in predictors of mammalian 
##### extinction risk: big is bad, but only in the tropics. 
##### Ecol. Lett. 12, 538–49. (doi:10.1111/j.1461-0248.2009.01307.x)
###############################################

library(picante)

fullcomm<-(as.matrix(t(row.per.sp)))
spname_rayita<-sub( " ", "_", colnames(fullcomm))
#### fix Puma yagouaroundi and Myotis 
spname_rayita[8]<-"Puma_concolor"
spname_rayita[17]<-"Artibeus_jamaicensis"

colnames(fullcomm)<-spname_rayita

mammal.tree<-read.nexus("data/MammalTree_nexus2.txt")
prunedphy <- prune.sample(fullcomm, mammal.tree[[1]])

# prunedphy
plot(prunedphy)

#############################
# community by forest type
# cam.habitat<-read.csv("data/cam_habitat.csv", header = T)
abun.sp<-as.data.frame(obs3)
abun.sp$habitat<-habitat3

humid.cam<-dplyr::filter(abun.sp, habitat3=="humid")
humid.com<-apply(humid.cam[,1:17],2,sum)

transition.cam<-dplyr::filter(abun.sp, habitat3=="transition")
transition.com<-apply(transition.cam[,1:17],2,sum)

dry.cam<-dplyr::filter(abun.sp, habitat3=="dry")
dry.com<-apply(dry.cam[,1:17],2,sum)

communi<-matrix(
          c(humid.com, transition.com, dry.com),
          nrow=3,              # number of rows 
          ncol=17,              # number of columns 
          byrow = TRUE)

colnames(communi)<-spname_rayita
rownames(communi)<-c("humid", "transition", "dry")



#############################
# traits mass guild

sp.traits<-cbind(guild, mass)
rownames(sp.traits)<-spname_rayita


# see species present in each community
par(mfrow = c(1, 3))
  for (i in row.names(communi)) {
  plot(prunedphy, show.tip.label = FALSE, main = i)
  tiplabels(tip = which(prunedphy$tip.label %in% names(which(communi [i, ] > 0))) , pch=19, cex=2)
  }


pd.result <- pd(communi, mammal.tree[[1]], include.root=TRUE)
pd.result

phydist <- cophenetic(mammal.tree[[1]]) #Cophenetic Distances for mammal Clustering
ses.mpd.result <- ses.mpd(communi, phydist, null.model="taxa.labels", abundance.weighted=TRUE, runs = 1000) #pairwise distances in communities
ses.mpd.result

ses.mntd.result <- ses.mntd(communi, phydist, null.model="taxa.labels", abundance.weighted=TRUE, runs = 1000) #pairwise distances in communities
ses.mntd.result
                          

############################################
### Geographic covariates
############################################

library(raster)
library(rgdal)
# library(dismo)
# library(biomod2)
library(spatstat)
library(sp)
library(dplyr)
library(maptools)


long<-unique(machalilla.fixed$longitude)
lati<-unique(machalilla.fixed$latitude)
centercoord<-c(mean(subset(long, long<=1)),mean(unique(subset(lati, lati<=1))))
coordsubset<-subset(machalilla.fixed,select = c(camera_trap,longitude,latitude,first_name_set_camera))

#################################
# get elevation
################################

# elevation<-getData('SRTM',lon=centercoord[1], lat=centercoord[2])

# read elevation from disk
elevation<- raster("C://Users//Diego//Documents//CodigoR//ULEAM//Infor_Caract//code//srtm_20_13.tif")

elevation2 <- readGDAL("C://Users//Diego//Documents//CodigoR//ULEAM//Infor_Caract//code//srtm_20_13.tif")
image(elevation, col= grey(1:99/100), axes=TRUE)
elevation<-raster(elevation2)


cam.cords<-as.data.frame(distinct(coordsubset))
coordinates(cam.cords) <- ~longitude+latitude #make sppatial data fram
geo <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") #def cord
proj4string(cam.cords)<-geo # set cords

e<-extent (-80.9, -80.5, -1.75,-1.35)
elevation.crop<-crop(elevation, e)

e = extent(raster(xmn = -80.9, xmx = -80.5, ymn = -1.75, ymx = -1.35))
cr1 = crop(elevation, e)



# plot(elevation.crop)
# plot(cam.cords, add=T, col="red")
# title(main="Altitud", sub="en color rojo se muestra donde se instalaron las camaras")

slope<-terrain(elevation.crop, opt='slope', unit='degrees', neighbors=4)
# plot(slope)
# plot(cam.cords, add=T, col="red")
# title(main="Pendiente", sub="en color rojo se muestra donde se instalaron las camaras")

cam.cords.sp<-SpatialPoints(cam.cords)
proj4string(cam.cords.sp)<-geo 
# etract values
elev.ovr <- extract(elevation.crop, cam.cords, method='bilinear')
slope.ovr <- extract(slope, cam.cords, method='bilinear')

# add to table
cam.cords$elev<-elev.ovr
cam.cords$slope<-slope.ovr

############################
# road
############################

roadpol <- readShapeSpatial("C:/Users/Diego/Documents/CodigoR/ULEAM/Infor_Caract/shp/machalilla_roadsclip.shp")
names(roadpol)<-"dist_rd"
proj4string(roadpol)<-geo
dist_rd<-over(x = cam.cords, y = roadpol)
# add to table
cam.cords$dist_rd<-as.numeric(dist_rd[,1])

roadpol.ow<-as(as(roadpol, "SpatialPolygons"), "owin") # make owin
cam.and.covs<-as.data.frame(cam.cords)

# plot(roadpol, col=topo.colors(65))
# plot(cam.cords, add=T, col="red")
# title(main="Distancia a las carreteras", sub="en color rojo se muestra donde se instalaron las camaras")


#####################################################
## Deforestation
####################################################

dist.def<-raster("C:/Users/Diego/Documents/CodigoR/ULEAM/Infor_Caract/Data/dist_def.tif")

# plot(dist.def)
# plot(cam.cords, add=T, col="red")
# title(main="Distancia a Deforestacion", sub="en color rojo se muestra donde se instalaron las camaras")

# plot(deforestado, col="red", add=T)
# etract values
dist.def.ovr <- extract(dist.def, cam.cords, method='bilinear')
index<-which(is.na(dist.def.ovr)) # detect NA. Means cam is in deforested
dist.def.ovr [index]<-0 # camera in distance cero to deforested


# add to table
cam.and.covs$dist_def<-as.numeric(dist.def.ovr)

cam.cords$dist_def<-dist.def.ovr
# am.and.covs<-as.data.frame(cam.cords)


##################################
## Estructura Veget
##################################

est.veget<-read.csv("C:/Users/Diego/Documents/CodigoR/ULEAM/Infor_Caract/Data/estructVeget.csv") # read table

cam.and.covs<-cbind(cam.and.covs, est.veget) #Paste covs


################################
cam.and.covs.scaled <-cbind(scale(cam.and.covs[5],center = T,scale = T),
                            scale(cam.and.covs[6],center = T,scale = T),
                            scale(cam.and.covs[7],center = T,scale = T),
                            scale(cam.and.covs[8],center = T,scale = T),
                            scale(cam.and.covs[10],center = T,scale = T),
                            scale(cam.and.covs[11],center = T,scale = T),
                            scale(cam.and.covs[12],center = T,scale = T))





