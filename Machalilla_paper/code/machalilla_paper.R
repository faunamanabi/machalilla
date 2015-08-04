

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
row.per.sp<-as.data.frame(matrix(nrow = length(sp.names), ncol=c(59)))
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
naiveoccu<-apply(X = presence,FUN = sum, MARGIN = 1 ) / 59 #divided number of sites
events<-apply(X = row.per.sp ,FUN = sum, MARGIN = 1 )
RAI<-(events/(114 * 59)) * 100 #114 is total sampling effort in days
Table1<-cbind(events, RAI, naiveoccu)
kable(Table1, format = "markdown")
############################################################
## Distribucion posterior de la riqueza de especies
############################################################

# Riqueza de especies y acumulación, modelando la ocurrencia y la detectabilidad. 
# Este análisis sigue el método de Dorazio et al. (2006).

source("code/MultiSpeciesSiteOcc.R")

X1 = as.matrix(row.per.sp) # col.per.sp por dias y row.per.sp por sitios (camaras)
nrepls = 90 #dias 
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




