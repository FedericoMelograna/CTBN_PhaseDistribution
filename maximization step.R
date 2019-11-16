
# adesso creiamo la funzine per lo step di maximization step --------------

mpr<-MEtotale
diag(mpr)<-0
diag(mpr)

Mx<-rowSums(MEtotale)
# Mx<-rowSums(mpr)
Mx


(qx<-Mx/ET)
##probelmi, vengono degli infiniti!!!
(qx<-Mx/(ET+0.01 ))



# theta_xxi<-Metotalexxi/Mx

thetaxxprimo<-MEtotale/Mx #teoricamente dovrebbero anche esserci le condizioni giuste

round(thetaxxprimo,2)
#questa parte qua è giusta, sono sbagliate le funzioni!!!
#prove
diag(x)
x=matrix(c(2,5,23,45),ncol=2)
# (diag(x)<-0)
x
x-diag(x)
rowSums(x)
xx<-c(25,50)
x/xx
