

library(lubridate)
# library(dplyr)
library(maptools)
# Load 'rgdal' package, which is used to read/write shapefiles and rasters
library(rgdal)
source("code/TEAM_code.R")
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

source("code/MultiSpeciesSiteOcc.R")

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

plot(sp2, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
plot(sp2)
# boxplot(sp1, col="yellow", add=TRUE, pch="+")

H <- diversity(t(row.per.sp))
simp <- diversity(t(row.per.sp), "simpson")
invsimp <- diversity(t(row.per.sp), "inv")
r.2 <- rarefy(t(row.per.sp), 2)
alpha <- fisher.alpha(t(row.per.sp))
pairs(cbind(H, simp, invsimp, r.2, alpha), pch="+", col="blue")

## Species richness (S) and Pielou's evenness (J):
S <- specnumber(t(row.per.sp)) ## rowSums(BCI > 0) does the same...
J <- H/log(S)

rarecurve(t(row.per.sp),step = 20, sample = raremax)

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

# community by forest type
cam.habitat<-read.csv("data/cam_habitat.csv", header = T)

# filter mat.per.sp to sp of wild mammals
mammal.per.sp<- mat.per.sp[-c(37,36,35,34,33,30,29,26,25,21,16,15,14,12,10,7,4,3,2,1)]
full.mammal<-ldply(mammal.per.sp, data.frame)



# convert list of species to unmarked 
# fullmat<-f.convert.to.unmarked(mammal.per.sp)
spcies<-as.data.frame( full.mammal[,1])
colnames(spcies)<-"sp"
obs<-as.data.frame(full.mammal[,2:162])

#  covariates of detection and occupancy in that order.
fullmat<-unmarkedFrameOccu(y=obs,siteCovs=spcies)

fm0 <- occu(~ 1 ~ 1, fullmat) 
fm1 <- occu(~ 1 ~ sp, fullmat)
fm2 <- occu(~ sp ~ 1, fullmat)
fm3 <- occu(~ sp ~ sp, fullmat)


models <- fitList(
  'p(.)psi(.)' = fm0,
  'p(.)psi(sp)' = fm1,
  'p(sp)psi(.)' = fm2,
  'p(sp)psi(sp)' = fm3)

ms <- modSel(models)
