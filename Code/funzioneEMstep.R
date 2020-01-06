
# funzione che fa un passo di EM ------------------------------------------

EMstep<-function(matricepartenza,storiaoss,Nv,k,N){
print("aa")
EEsemplice=funzioneETtotalesigma(matrice=matricepartenza,storia=storiaoss,Nv=Nv,k=k,N=N)
print("ciao")
metotalesemplice=funzMEtotale(matrice=matricepartenza,storia=storiaoss,Nv=Nv,k=k)
print("ciao2")
qiniziale=-diag(matricepartenza)
thetainiz=matricepartenza/qiniziale
mpr<-metotalesemplice 
diag(mpr)<-0
Mx<-rowSums(mpr) 
EExsemplice<-as.vector(EEsemplice)
qx<-Mx/EExsemplice
thetaxx<-metotalesemplice/Mx 
nuovamatrice=matrix(0,nrow=nrow(matricepartenza),ncol=nrow(matricepartenza))
diag(nuovamatrice)=qx
for (i in 1:nrow(nuovamatrice)){
  for (j in 1:ncol(nuovamatrice)){
    if (i != j){
      nuovamatrice[i,j]<-thetaxx[i,j]*qx[i]
    }
  }
}
diag(nuovamatrice)=-diag(nuovamatrice)
provaperlogalt=0
for (i in 1:nrow(metotalesemplice)){
  for (j in 1:ncol(metotalesemplice)){
    if (i !=j){
      temp=metotalesemplice[i,j]*log(thetainiz[i,j])
      if (temp !="NaN"){ 
        
        provaperlogalt=provaperlogalt+metotalesemplice[i,j]*log(thetainiz[i,j])
      }
    }}
}
LogVerononNorm=sum(Mx*log(qiniziale)-qiniziale*as.vector(EEsemplice))+provaperlogalt
#[1] [1]--5.751853e-08

TotET=sum(EEsemplice) #5.613631e-08
LogVero=LogVerononNorm/sum(EEsemplice)
# [1] -1.024622
Oggettodaritornare<-list(TotET=TotET,LogVero=LogVero,LogVerononNorm=LogVerononNorm,
                         ET=EEsemplice,Mejk=metotalesemplice,matricefinale=nuovamatrice)
return(Oggettodaritornare)

}
# 
# provone=list(tempo=EEsemplice,matrice=Qsempl)
# provone

start=Sys.time()

(safa=EMstep(matricepartenza =Qsempl ,storiaoss =eustoria2 ,Nv =2 ,k =2 ,N =2)
)
end=Sys.time()
end-start

