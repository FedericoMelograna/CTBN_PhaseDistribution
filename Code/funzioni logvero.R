calcoloverophase_test=function(matrice,storia_test,Nv,k,N){

    
    EEsemplice=funzioneETtotalesigma(matrice=matrice,storia=storia_test$storiafinale,Nv=Nv,k=k,N=N)
    # print("ciao")
    metotalesemplice=funzMEtotale(matrice=matrice,storia=storia_test$storiafinale,Nv=Nv,k=k)
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
    LogVero=LogVerononNorm/sum(EEsemplice[,,i])*tempotot[i]
    logv=logv+LogVero
  }
  
  return(logv)
}

calcoloverophase_test_nony=calcoloverophase_test=function(matrice,storia_test,Nv,k,N){
  
  
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
