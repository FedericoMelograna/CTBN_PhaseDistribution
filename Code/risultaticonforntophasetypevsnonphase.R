#prove risultati migliore impostazione è la mia 
set.seed(123)
Q_risultati=matrix(c(-2,2,0,0,1.5,-2,0.5,0,0,0,-7,7,4,0,3,-7),nrow=4,byrow=T)
Q_partenza=Q_risultati
random=runif(length(Q_partenza[Q_partenza>0]),min=0.5,max=2.5)
Q_partenza[Q_partenza>0]=random
diag(Q_partenza)<-0
diag(Q_partenza)=-rowSums(Q_partenza)
storia_ris=generastorie(matrice =Q_risultati ,tmax =200 ,statix =2 ,Nv =2 ,N =2 )
storia_ris_test=generastorie(matrice =Q_risultati ,tmax =100 ,statix =2 ,Nv =2 ,N =2 )



start=Sys.time()
stima_risultati=EMstep(matricepartenza =Q_partenza ,storiaoss =storia_ris$storiafinale[,c(2,3,4)] ,Nv =2 ,k =2 ,N =2)
end=Sys.time()


start=Sys.time()
stima_risultati2=EMstep(matricepartenza =stima_risultati$matricefinale ,storiaoss =storia_ris$storiafinale[,c(2,3,4)] ,Nv =2 ,k =2 ,N =2)
end=Sys.time()
end-start
# $LogVero
# [1] -0.5247905


EEsemplice=funzioneETtotalesigma(matrice=stima_risultati2$matricefinale,storia=storia_ris_test$storiafinale,Nv=2,k=2,N=2)
# print("ciao")
metotalesemplice=funzMEtotale(matrice=stima_risultati2$matricefinale,storia=storia_ris_test$storiafinale,Nv=2,k=2)
# print("ciao2")

###Modifiche!
sm=sum(EEsemplice)
tempotot=storia_ris_test$storiafinale$fine[nrow(storia_ris_test$storiafinale)]-storia_ris_test$storiafinale$inizio[1]
EEsemplice=(EEsemplice/sm)*tempotot
metotalesemplice=(metotalesemplice/sm)*tempotot
##fine modifiche

qiniziale=-diag(stima_risultati2$matricefinale)
thetainiz=stima_risultati2$matricefinale/qiniziale
mpr<-metotalesemplice 
diag(mpr)<-0
Mx<-rowSums(mpr) 
EExsemplice<-as.vector(EEsemplice)
qx<-Mx/EExsemplice
thetaxx<-metotalesemplice/Mx 
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
#[1] [1]-44.90749

TotET=sum(EEsemplice) #5.613631e-08
LogVero=LogVerononNorm/sum(EEsemplice)
#-0.4490749
calcoloverophase_test=function(matrice,storia_test,Nv,k,N){
  
  EEsemplice=funzioneETtotalesigma(matrice=matrice,storia=storia_test,Nv=Nv,k=k,N=N)
  # print("ciao")
  metotalesemplice=funzMEtotale(matrice=matrice,storia=storia_test,Nv=Nv,k=k)
  # print("ciao2")
  
  ###Modifiche!
  sm=sum(EEsemplice)
  tempotot=storia_test$storiafinale$fine[nrow(storia_test$storiafinale)]-storia_test$storiafinale$inizio[1]
  EEsemplice=(EEsemplice/sm)*tempotot
  metotalesemplice=(metotalesemplice/sm)*tempotot
  ##fine modifiche
  qiniziale=-diag(matrice)
  thetainiz=matrice/qiniziale
  mpr<-metotalesemplice 
  diag(mpr)<-0
  Mx<-rowSums(mpr) 
  EExsemplice<-as.vector(EEsemplice)
  qx<-Mx/EExsemplice
  thetaxx<-metotalesemplice/Mx 
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
  #[1] [1]-44.90749
  
  TotET=sum(EEsemplice) #5.613631e-08
  LogVero=LogVerononNorm/sum(EEsemplice)
  return(LogVerononNorm)
}
lv_phase=calcoloverophase_test(matrice=stima_risultati2$matricefinale,storia_test = storia_ris_test$storiafinale,Nv=2,k=2,N=2)

###
Q_partenza_nonphase=matrix(c(-1,1,2,-2),nrow=2,ncol=2,byrow=T)
EMstep(matricepartenza =stima_risultati$matricefinale ,storiaoss =storia_ris$storiafinale[,c(2,3,4)] ,Nv =1 ,k =2 ,N =1)

primo=mean(storia_ris$storiafinale$fine[storia_ris$storiafinale$stato==1]-storia_ris$storiafinale$inizio[storia_ris$storiafinale$stato==1])
#4.183277
secondo=mean(storia_ris$storiafinale$fine[storia_ris$storiafinale$stato==2]-storia_ris$storiafinale$inizio[storia_ris$storiafinale$stato==2])
#0.4678859
Q_nonphase=matrix(c(-1/primo,1/primo,1/secondo,-1/secondo),nrow=2,ncol=2,byrow=T)

calcolologveronp_binario=function(Q_nonphase,storia_ris_test){#storia_ristest è storia finale
  
  M1=ifelse(nrow(storia_ris_test)%%2!=0,
            nrow(storia_ris_test)/2,
            ifelse(storia_ris_test$stato[1]==1,
                   nrow(storia_ris_test)/2,
                   (nrow(storia_ris_test)/2)-1))
  M2=ifelse(nrow(storia_ris_test)%%2!=0,
            nrow(storia_ris_test)/2,
            ifelse(storia_ris_test$stato[1]==2,
                   nrow(storia_ris_test)/2,
                   (nrow(storia_ris_test)/2)-1))
  lq1np=log(-Q_nonphase[1,1])
  lq2np=log(-Q_nonphase[2,2])
  q1np=-Q_nonphase[1,1]; q2np= -Q_nonphase[2,2]
  t1np=sum(storia_ris_test$fine[storia_ris_test$stato==1]-storia_ris_test$inizio[storia_ris_test$stato==1])
  t2np=sum(storia_ris_test$fine[storia_ris_test$stato==2]-storia_ris_test$inizio[storia_ris_test$stato==2])
  logv=M1*lq1np+M2*lq2np-q1np*t1np-q2np*t2np
  return(logv)
}
calcolologveronp_binario(Q_nonphase=Q_nonphase,storia_ris_test = storia_ris_test$storiafinale)
###-62.9


# confornto con parenty ---------------------------------------------------
Q_risultati2=matrix(c(-7,7,0,0,4,-7,3,0,0,0,-2,2,0.5,0,1,-1.5),nrow=4,byrow=T)


set.seed(1231)
Q_partenza2=Q_risultati2
random=runif(length(Q_partenza2[Q_partenza2>0]),min=0.5,max=2.5)
Q_partenza2[Q_partenza2>0]=random
diag(Q_partenza2)<-0
diag(Q_partenza2)=-rowSums(Q_partenza2)
matrice_y=matrix(c(-0.1,0.1,0.2,-0.2),byrow=T,nrow=2)


#storie
storia_ris_y=generastorie_y(matrice1 = Q_risultati,matrice2=Q_risultati2,matrice_y =matrice_y,
                     tmax =400,statix =2 ,Nv =2 ,N =2 )
storia_ris_senzay=fconvertistoria(storia_ris_y$storiafinale[,c(2,3,4,5)])

storia_ris_test_y=generastorie_y(matrice1 = Q_risultati,matrice2=Q_risultati2,matrice_y =matrice_y,
                                 tmax =300,statix =2 ,Nv =2 ,N =2 )
storia_ris_test_senzay=fconvertistoria(storia_ris_test_y$storiafinale[,c(2,3,4,5)])
###
###Learning parametri senza y
start=Sys.time()
stima_risultati_senza_y=EMstep(matricepartenza =Q_partenza ,storiaoss =storia_ris_senzay ,Nv =2 ,k =2 ,N =2)
end=Sys.time()
end-start

start=Sys.time()
stima_risultati_senza_y2=EMstep(matricepartenza =stima_risultati_senza_y$matricefinale ,storiaoss =storia_ris_senzay ,Nv =2 ,k =2 ,N =2)
end=Sys.time()
end-start
###fine learning senza y

##stima learning y
start=Sys.time()
stima_risultati_y=EMstep_y(matricipartenza = list(matrice1=Q_partenza,
         matrice2=Q_partenza2), storiaoss =storia_ris_y$storiafinale[,c(2,3,4,5)] ,
         Nv =2 ,k =2 ,N = 2)
end=Sys.time()
end-start

start=Sys.time()
stima_risultati_y2=EMstep_y(matricipartenza = list(matrice1=stima_risultati_y[[1]]$matricefinale,
                          matrice2=stima_risultati_y[[2]]$matricefinale), storiaoss =storia_ris_y$storiafinale[,c(2,3,4,5)] ,
                           Nv =2 ,k =2 ,N = 2)
end=Sys.time()
end-start
#########

logv_senzay=calcoloverophase_test(matrice =stima_risultati_senza_y2$matricefinale ,storia_test =storia_ris_test_senzay ,Nv =2 ,k =2 ,N = 2)
# 2279.851
matricipartenza=list(matrice1=stima_risultati_y2[[1]]$matricefinale,
                     matrice2=stima_risultati_y2[[2]]$matricefinale)
storia_test_y=storia_ris_test_y$storiafinale[,c(2,3,4,5)]
calcoloverophase_test_y<-function(matricipartenza,storia_test_y,Nv,k,N){

EEsemplice=funzioneETtotalesigma_y(matrici=matricipartenza,storia=storia_test_y
                                   ,Nv=Nv,k=k,N=N)
# print("ciao")
metotalesemplice=funzMEtotale_y(matrici=matricipartenza,storia=storia_test_y,Nv=Nv,k=k)
# print("ciao2")


###modifiche!!!
for(i in 1:length(matricipartenza)){
  sm[i]=sum(EEsemplice[,,i])
  tempotot[i]=sum(storia_test_y$fine[storia_test_y$y==i]-storia_test_y$inizio[storia_test_y$y==i])
  EEsemplice[,,i]=(EEsemplice[,,i]/sm[i])*tempotot[i]
  metotalesemplice[,,i]=(metotalesemplice[,,i]/sm[i])*tempotot[i]
}

###fine modifiche
i=0
logv=0
for (matricepartenza in matricipartenza){
  i=i+1
  qiniziale=-diag(matricepartenza)
  thetainiz=matricepartenza/qiniziale
  mpr<-metotalesemplice[,,i] 
  diag(mpr)<-0
  Mx<-rowSums(mpr) 
  EExsemplice<-as.vector(EEsemplice[,,i])
  qx<-Mx/EExsemplice
  thetaxx<-metotalesemplice[,,i]/Mx

  provaperlogalt=0
  for (l in 1:nrow(metotalesemplice)){
    for (j in 1:ncol(metotalesemplice)){
      if (l !=j){
        temp=metotalesemplice[l,j,i]*log(thetainiz[l,j])
        if (temp !="NaN"){ 
          
          provaperlogalt=provaperlogalt+metotalesemplice[l,j,i]*log(thetainiz[l,j])
        }
      }}
  }
  LogVerononNorm=sum(Mx*log(qiniziale)-qiniziale*as.vector(EEsemplice[,,i]))+provaperlogalt
  #[1] [1]--5.751853e-08
  
  TotET=sum(EEsemplice[,,i]) #5.613631e-08
  LogVero=LogVerononNorm/sum(EEsemplice[,,i])
  logv=logv+LogVerononNorm
}

return(logv)
}
logv_y=calcoloverophase_test_y(matricipartenza = list(matrice1=stima_risultati_y2[[1]]$matricefinale,
      matrice2=stima_risultati_y2[[2]]$matricefinale), storia_test_y=storia_ris_test_y$storiafinale[,c(2,3,4,5)],2,2,2 )


##[1] 5238.047
calcoloverophase_test_nony= function(matrice,storia_test,Nv,k,N){
  
  
  EEsemplice=funzioneETtotalesigma(matrice=matrice,storia=storia_test,Nv=Nv,k=k,N=N)
  # print("ciao")
  metotalesemplice=funzMEtotale(matrice=matrice,storia=storia_test,Nv=Nv,k=k)
  # print("ciao2")
  
  ###Modifiche!
  sm=sum(EEsemplice)
  tempotot=storia_test$fine[nrow(storia_test)]-storia_test$inizio[1]
  EEsemplice=(EEsemplice/sm)*tempotot
  metotalesemplice=(metotalesemplice/sm)*tempotot
  ##fine modifiche
  qiniziale=-diag(matrice)
  thetainiz=matrice/qiniziale
  mpr<-metotalesemplice 
  diag(mpr)<-0
  Mx<-rowSums(mpr) 
  EExsemplice<-as.vector(EEsemplice)
  qx<-Mx/EExsemplice
  thetaxx<-metotalesemplice/Mx 
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
  #[1] [1]-44.90749
  
  TotET=sum(EEsemplice) #5.613631e-08
  LogVero=LogVerononNorm/sum(EEsemplice)
  return(LogVerononNorm)
}


calcoloverophase_test_y=function(matricipartenza,storia_test_y,Nv,k,N){
  
  EEsemplice=funzioneETtotalesigma_y(matrici=matricipartenza,storia=storia_test_y
                                     ,Nv=Nv,k=k,N=N)
  # print("ciao")
  metotalesemplice=funzMEtotale_y(matrici=matricipartenza,storia=storia_test_y,Nv=Nv,k=k)
  # print("ciao2")
  
  
  ###modifiche!!!
  for(i in 1:length(matricipartenza)){
    sm[i]=sum(EEsemplice[,,i])
    tempotot[i]=sum(storia_test_y$fine[storia_test_y$y==i]-storia_test_y$inizio[storia_test_y$y==i])
    EEsemplice[,,i]=(EEsemplice[,,i]/sm[i])*tempotot[i]
    metotalesemplice[,,i]=(metotalesemplice[,,i]/sm[i])*tempotot[i]
  }
  
  ###fine modifiche
  i=0
  logv=0
  for (matricepartenza in matricipartenza){
    i=i+1
    qiniziale=-diag(matricepartenza)
    thetainiz=matricepartenza/qiniziale
    mpr<-metotalesemplice[,,i] 
    diag(mpr)<-0
    Mx<-rowSums(mpr) 
    EExsemplice<-as.vector(EEsemplice[,,i])
    qx<-Mx/EExsemplice
    thetaxx<-metotalesemplice[,,i]/Mx
    
    provaperlogalt=0
    for (l in 1:nrow(metotalesemplice)){
      for (j in 1:ncol(metotalesemplice)){
        if (l !=j){
          temp=metotalesemplice[l,j,i]*log(thetainiz[l,j])
          if (temp !="NaN"){ 
            
            provaperlogalt=provaperlogalt+metotalesemplice[l,j,i]*log(thetainiz[l,j])
          }
        }}
    }
    LogVerononNorm=sum(Mx*log(qiniziale)-qiniziale*as.vector(EEsemplice[,,i]))+provaperlogalt
    #[1] [1]--5.751853e-08
    
    TotET=sum(EEsemplice[,,i]) #5.613631e-08
    LogVero=LogVerononNorm/sum(EEsemplice[,,i])*tempotot[i]
    logv=logv+LogVero
  }
  
  return(logv)
}
save.image(file="environment2.RData")

