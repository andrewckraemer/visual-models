### Visual Model for the Blue Tit
# Written by Andrew C. Kraemer (procedure adapted from Maan and Cummings 2009,2012; Siddiqi et al 2004; Vorobyev et al 1998,2001); most of the terminology follows the work of Maan and Cummings

rm(list=ls()) # this helps to clear R (I like to start with a blank slate)

##### Input data
# Weber Fraction: v
# fraction of photoreceptor types: nu, ns, nme, nl 
# photoreceptor noise: wu, ws, wm, wl = v/n

refTar<-read.csv("refTspec.csv") # target reflectance spectra (animals): refTar
refBac<-read.csv("refBspec.csv") # background reflectance spectra (leaves, soil, etc.): refBac
irr<-read.csv("irr.csv") # habitat irradiance (the amount of light available to receiver): irr
btit<-read.csv("btit.csv") # taxon-specific absorptance spectra for each photoreceptor class in each species: abUbtit,abSbtit,abMbtit,abLbtit for the Blue Tit
abUbtit<-as.matrix(btit$UV)
abSbtit<-as.matrix(btit$short)
abMbtit<-as.matrix(btit$medium)
abLbtit<-as.matrix(btit$long)
v<-0.1
nu<-0.076
ns<-0.146
nme<-0.204
nl<-0.574
wu<-v/nu
ws<-v/ns
wm<-v/nme
wl<-v/nl

#Quantum Catch
QuTar<-refTar%*%as.matrix(irr*abUbtit)
QuBac<-refBac%*%as.matrix(irr*abUbtit)

QsTar<-refTar%*%as.matrix(irr*abSbtit)
QsBac<-refBac%*%as.matrix(irr*abSbtit)

QmTar<-refTar%*%as.matrix(irr*abMbtit)
QmBac<-refBac%*%as.matrix(irr*abMbtit)

QlTar<-refTar%*%as.matrix(irr*abLbtit)
QlBac<-refBac%*%as.matrix(irr*abLbtit)

#Von Kries Transformation
QuIrr<-t(irr)%*%abUbtit
QsIrr<-t(irr)%*%abSbtit
QmIrr<-t(irr)%*%abMbtit
QlIrr<-t(irr)%*%abLbtit
nspec<-length(QuTar)
quTar<-NULL
for(i in 1:nspec){
	quTar.temp<-QuTar[i]/QuIrr
	quTar<-rbind(quTar,quTar.temp)
}
qsTar<-NULL
for(i in 1:nspec){
	qsTar.temp<-QsTar[i]/QsIrr
	qsTar<-rbind(qsTar,qsTar.temp)
}
qmTar<-NULL
for(i in 1:nspec){
	qmTar.temp<-QmTar[i]/QmIrr
	qmTar<-rbind(qmTar,qmTar.temp)
}
qlTar<-NULL
for(i in 1:nspec){
	qlTar.temp<-QlTar[i]/QlIrr
	qlTar<-rbind(qlTar,qlTar.temp)
}

quBac<-rep(QuBac/QuIrr,length(quTar))
qsBac<-rep(QsBac/QsIrr,length(quTar))
qmBac<-rep(QmBac/QmIrr,length(quTar))
qlBac<-rep(QlBac/QlIrr,length(quTar))

#Contrast between target and background (photoreceptor-specific)
fu<-log(abs(quTar-quBac))
fs<-log(abs(qsTar-qsBac))
fm<-log(abs(qmTar-qmBac))
fl<-log(abs(qlTar-qlBac))

#Luminosity Contrast
LconBtit<-abs(fl/wl)

#Spectral Contrast
SconBtit<-sqrt(
((wu*ws)^2*(fl-fm)^2
+(wu*wm)^2*(fl-fs)^2
+(wu*wl)^2*(fm-fs)^2
+(ws*wm)^2*(fl-fu)^2
+(ws*wl)^2*(fm-fu)^2
+(wm*wl)^2*(fs-fu)^2)/
((wu*ws*wm)^2
+(wu*ws*wl)^2
+(wu*wm*wl)^2
+(ws*wm*wl)^2))