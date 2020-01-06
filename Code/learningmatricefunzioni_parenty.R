
# funzioni ----------------------------------------------------------------

# primabozzalearningparametri ---------------------------------------------
library("expm", lib.loc="~/R/win-library/3.4")
library("dplyr", lib.loc="~/R/win-library/3.4")
#matrice che da Qsi;si
Qsif=function(iesima,matrice,Nv){ #iesima sta a significare quale sottostato siamo
  ###matrice=matrice di partenza 
  Qsi=matrix(NA,nrow=nrow(matrice),ncol=ncol(matrice))
  for (i in 1:nrow(matrice)){
    for(j in 1:ncol(matrice)){
      Qsi[i,j]=ifelse(j>(iesima-1)*Nv && i>(iesima-1)*Nv && i<=iesima*Nv && j<=iesima*Nv,matrice[i,j],0)
    }
  } 
  return(Qsi)
}
##matrice che da Q si0;si1
Qsi0si1f=function(iesima0,iesima1,matrice,Nv){
  #funzione: da dove parte, 0, dove arriva, 1 la matrice di partenza, e Nv. 
  Qsi0si1=matrix(NA,nrow=nrow(matrice),ncol=ncol(matrice))
  for (i in 1:nrow(matrice)){
    for(j in 1:ncol(matrice)){
      Qsi0si1[i,j]=ifelse(j>(iesima1-1)*Nv && i>(iesima0-1)*Nv && i<=iesima0*Nv && j<=iesima1*Nv,matrice[i,j],0)
    }
  } 
  return(Qsi0si1)
  #restituisce la matrice con zero ovunque meno che nella transizione tra stato 0 a 1
}

#singolo passo: da alpha t ---> alpha t+1
alphatpiu1<-function(alphat,Qsi,ti1,ti0,Qsi0_si1){
  #mi computa la formula tramite la quale passo da alphat ad alphat+1
  alphat%*%expm(Qsi*(ti1-ti0))%*%Qsi0_si1
}
#singolo passo: da beta t+1--> beta t
betatif<-function(betat1,Qsi,ti1,ti0,Qsi0_si1){
  ###prenda il beta precedente, la matrice Qsi con zero altrove 
  #la matrice della transizione Qsi0_si1 e i due tempi dove avviene e computa il nuovo beta. 
  betat0<-Qsi0_si1%*%expm(Qsi*(ti1-ti0))%*%(betat1)
  return(betat0)
  #NB Qsi0,si1 sono in realta si-1 , si ; mentre abbiamo ti+1-ti
}

#beta
betacumulatofn_y<-function(matrici,storia,Nv,k ){
  ##### da controllare bene come abbiamo messo i e i-1 ad alpha e beta qui 
  #restituisce un vettore di beta, da beta1 a beta_tau  
  betacumulato<-array(1,c(k*Nv,nrow(storia)+1))
  #vettore di partenza
  #da 9 a 2
  for (i in (nrow(storia)):2){ #ciclo al contrario, per mettere dentro gli elementi ricorsivamente
    iesimainiz=storia$stato[i-1] #si-1
    iesimafin=storia$stato[i] #si
    # print(c(storia$stato[i-1],storia$stato[i]))
    # print("-----")
    matrice=matrici[paste0("matrice",storia$y[i])][[1]]
    #qui è i e non i-1 in quanto sto considerando iesima finale come storia
    ## perche va da sx <---- quindi si fa tutto finale, e poi transiziona in 
    #iniziale
    if (iesimainiz!=iesimafin){
      
    
    Qsi<-Qsif(iesima=iesimafin,matrice=matrice,Nv=Nv) #Qsi
    Qsi1si0<-Qsi0si1f(iesima0 =iesimainiz ,iesima1 =iesimafin, matrice = matrice,Nv=Nv)
    #Qsi-1,Qsi
    ti1<-storia$fine[i] #ti+1?
    # print(storia$fine[i+1])
    ti0<-storia$inizio[i] #ti?
    # print(c(storia$inizio[i],storia$fine[i]))
    betacumulato[,i]<-betatif(betat1 =betacumulato[,i+1],Qsi = Qsi,ti1 = ti1
                              ,ti0= ti0,Qsi0_si1 =Qsi1si0  )
    #computo la funzione
    }
    if (iesimainiz==iesimafin){
      Qsi<-Qsif(iesima=iesimafin,matrice=matrice,Nv=Nv) #da 2.7 a 3
      # Qsi1si0<-Qsi0si1f(iesima0 =iesimainiz ,iesima1 =iesimafin, matrice = matrice,Nv=Nv)
      #ti
      ti1<-storia$fine[i] #3
      # print(storia$fine[i+1])
      #ti-1
      ti0<-storia$inizio[i] #2.7
      # print(c(storia$inizio[i],storia$fine[i]))
      #sarebbe betatau, giusto 
      #beta(t)*Qsi-1*(ti-t(i-1))
      betacumulato[,i]<-betatifpiu(betat1 =betacumulato[,i+1],Qsi = Qsi,ti1 = ti1
                                   ,ti0= ti0 )
      
    }
  }
  return(betacumulato)
}
# bbbsemplice=betacumulatofn(matrice=Qsempl,storia=eustoria2,Nv=2,k=2)
# bbbsemplicey=betacumulatofn_y(matrice=Qsempl,storia=eustoria2,Nv=2,k=2)
# sum(bbbsemplice-bbbsemplicey)

betacumulatof0_y<-function(matrici,storia,betacumulprimo,Nv,k){ 
  #crea beta0 da sovrascrivere poi
  beta0<-array(1,c(k*Nv,1))
  iesima=storia$stato[1] #sia fin che iniz
  #S0
  # # iesimafin=storia$stato[i+1]
  # print(storia$stato[1])
  # print("-----")
  matrice=matrici[paste0("matrice",storia$y[1])][[1]]
  Qsi<-Qsif(iesima=iesima,matrice=matrice,Nv=Nv)
  #Qso
  ti1<-storia$fine[1]
  #t1
  # print(storia$fine[i+1])
  ti0<-storia$inizio[1]
  #t0=0
  # print(c(storia$inizio[1],storia$fine[1]))
  beta0<-betatifpiu(betat1 =betacumulprimo,Qsi = Qsi,ti1 = ti1
                    ,ti0= ti0)
  return(beta0)
}
betacumfinalenonalt_y<-function(matrici,storia,Nv,k){
  #crea vettore beta giusto
  bbb=betacumulatofn_y(matrici=matrici,storia=storia,Nv=Nv,k=k)
  bbb0=betacumulatof0_y(matrici=matrici,storia=storia,betacumulprimo = bbb[,2],Nv=Nv,k=k)
  bbb[,1]<-bbb0
  bbbfin<-bbb
  return(bbbfin)
  
}
bbb_y=betacumfinalenonalt_y(matrici =list(matrice1=Qsempl,matrice2=ssss) ,storia =
                   asdas$storiafinale[,c(2,3,4,5)],Nv =2 ,k =2 )
bbb_y[,c(1:5,130:135)] ##mm strano uno sopra a 1
asdas$storiafinale[c(100:134),]
# 129       4 94.46872  94.54695 2     2
# 130       2 94.54695  95.31737 2     1
# 131       4 95.31737  96.92016 2     2
# 132       2 96.92016  97.25810 2     1
# 133       4 97.25810  98.07106 2     2
# 134       2 98.07106  98.09058 2     1
# 135       4 98.09058 100.00000 2     2
#da beta0 a beta_tau;####versione non alternativa cumulata

#alpha

alphacumulatofn_y<-function(matrici,storia,Nv,k ){
  #computa alpha da 0 (inizializzata da Po che non conosco), finisco a tau-1
  # alpha0<-rep(0.5,k*Nv) #alphao
  alpha0<-rep(0,k*Nv)
  ##probabile che sia meglio inizializzare solo per gli stati non transienti!!
  alpha0[(Nv*(storia$stato[1]-1)+1):(Nv*storia$stato[1])]<-1/Nv
  #alphacumulato<-array(0.5,c(k*Nv,nrow(storia)+1))
  alphacumulato<-cbind(alpha0,array(0.5,c(k*Nv,nrow(storia))))
  alphacumulato<-as.matrix(alphacumulato)
  #da alpha1 a alpha(tau+1)
  for (i in 2:(nrow(storia))){#ciclo che parte dal secondo posto 
    matrice=matrici[paste0("matrice",storia$y[i-1])][[1]]
    iesimainiz=storia$stato[i-1] #Si-2
    iesimafin=storia$stato[i] #Si-1
    # print(c(storia$stato[i-1],storia$stato[i]))
    # print("-----")
    if(iesimafin!=iesimainiz){
    Qsi<-Qsif(iesima=iesimainiz,matrice=matrice,Nv=Nv)
    #Qsi-2
    Qsi1si0<-Qsi0si1f(iesima0 =iesimainiz ,iesima1 =iesimafin, matrice = matrice,Nv=Nv)
    #Qsi-2,si-1
    ti1<-storia$fine[i-1]
    #ti-1
    # print(storia$fine[i+1])
    ti0<-storia$inizio[i-1]
    #ti-2
    
    # print(c(storia$inizio[i-1],storia$fine[i-1]))
    alphacumulato[,i]<-alphatpiu1(alphat =alphacumulato[,i-1],Qsi = Qsi,ti1 = ti1
                                  ,ti0= ti0,Qsi0_si1 =Qsi1si0  )
    #in alpha[i] ci va alpha[i-1];Qsi-2; Qsi-2,si-1; ti-1; ti-2
    #quindi si parta con alpha i-2 
    #e quella che trovo che alpha i-1
    
    }
    if(iesimafin==iesimainiz){
      Qsi<-Qsif(iesima=iesimainiz,matrice=matrice,Nv=Nv)
      # Qsi1si0<-Qsi0si1f(iesima0 =iesimainiz ,iesima1 =iesimafin, matrice = matrice,Nv=Nv)
      ti1<-storia$fine[i-1] #ti-1
      # print(storia$fine[i+1])
      ti0<-storia$inizio[i-1] #ti-2
      # print(c(storia$inizio[i-1],storia$fine[i-1]))
      
      #computa in posizione i l'alpha(ti-1)- perche prende alpha (ti-2)*Qsi-2* ti-1; ti-2
      alphacumulato[,i]<-alphatpiu1meno(alphat =alphacumulato[,i-1],Qsi = Qsi,ti1 = ti1
                                        ,ti0= ti0)
      
      
    }
  } 
  return(alphacumulato)
  #dimensioni(1,10)---> le mappero in 0, tau=9
}
aaa
aaa2=alphacumulatofn_y(matrici =list(matrice1=Qsempl,matrice2=ssss) ,storia =asdas$storiafinale[,c(2,3,4,5)] ,Nv =2 ,k =2 )
aaa[,c(1:5,130:135)]
aaa2[,c(1:5,130:135)]
aaa[,c(9:15)]
aaa2[,c(9:15)]

asdas=generastorie_y(matrice1 = Qsempl,matrice2 = ssss,matrice_y = matrice_y,
                     tmax = 100,statix = 2,Nv=2,N=2)

alphacumulftau_y<-function(matrici,storia,alphacumultaumenouno,Nv,k){ 
  #alphatau; #poi lo sostituirei in alphacumulato[10]
  alphatau<-array(1,c(k*Nv,1)) #vettore iniziale
  iesima=storia$stato[nrow(storia)] #sia fin che iniz
  ###ultimo stato, stato S(tau-1)
  
  # # iesimafin=storia$stato[i+1]
  # print(storia$stato[nrow(storia)])
  # print("-----")
  matrice=matrici[paste0("matrice",storia$y[nrow(storia)])][[1]]
  Qsi<-Qsif(iesima=iesima,matrice=matrice,Nv=Nv)
  #Qs(tau-1)
  ti1<-storia$fine[nrow(storia)]
  #tau
  # print(storia$fine[i+1])
  ti0<-storia$inizio[nrow(storia)]
  #tau-1
  # print(c(storia$inizio[1],storia$fine[1]))
  alphatau<-alphatpiu1meno(alphat =alphacumultaumenouno,Qsi = Qsi,ti1 = ti1
                           ,ti0= ti0)
  return(alphatau)
  #ritorna il vettore alpha(tau)
}
alphacumulftau_y(matrici =list(matrice1=Qsempl,matrice2=ssss) ,storia =
                   asdas$storiafinale[,c(2,3,4,5)],alphacumultaumenouno =aaa2[,ncol(aaa2)-1] ,Nv =2 ,k =2 )
alphacumfinalef_y<-function(matrici,storia,Nv,k){
  aaa=alphacumulatofn_y(matrici=matrici,storia=storia,Nv=Nv,k=k)
  aaa0=alphacumulftau_y(matrici=matrici,storia=storia,alphacumultaumenouno = aaa[,ncol(aaa)-1],Nv=Nv,k=k)
  #sostituisci ad alpha_tau<-- la sua rappresentazione con alphatau
  aaa[,ncol(aaa)]<-aaa0
  colnames(aaa)<-c(0,storia$fine)
  aaafin<-aaa
  return(aaafin)
  
}
alphacumfinalef_y(matrici =list(matrice1=Qsempl,matrice2=ssss) ,storia =
                  asdas$storiafinale[,c(2,3,4,5)],Nv =2,k =2 )
##mette assieme i due alpha e computa alfa globale da 0 a tau


##########a




#due metodi per integrare all'interno di un singolo intervallino:
composite.simpson2 <- function(Qsip,a,b,alpha,beta,Deltaj,n){ #,Qsip,alpha,beta,Deltaj) {#a e b sono low e up
  h <- (b - a) / n
  xj <- seq.int(a, b, length.out = n + 1)
  w=rep(1:length(xj))
  w<-ifelse(w%%2==0,4,2)
  w[1]=1; w[length(w)]=1
  ET=0
  for(i  in 1:length(xj)){
    ET=ET+(alpha%*%expm(Qsip*(xj[i]-a))%*%Deltaj%*%expm(Qsip*(b-xj[i]))%*%beta)*w[i]
  }
  approx <- (h / 3) * ET      
  return(approx)
  
} #questa è prova in provepervederesimposon
integralinosingoloET<-function(Qsip,low,up,alpha,beta,Deltaj){
  Etj=0 #mi computa integrale che ho dentro alla parte di E[t] e anche E[m]
  for (t in seq(low,up,length.out = 10)){ #meglio un numero piu basso?
    #formula
    Etj=Etj+(alpha%*%expm(Qsip*(t-low))%*%Deltaj%*%expm(Qsip*(up-t))%*%beta)/10
    #alpha di low, beta di high, Qs compreso tra low e high. Delta della transizione 
  }
  Etj=Etj*(up-low)
  return(Etj)
}
##useremo composite.simpson2 perche è piu efficiente


matrici
#integra per un certo x selezionato su tutta la storai
funzioneintegraleXselh_y<-function(matrici,storia,xsel,Nv,k){
  
  # xsel=1; 
  Etj=array(0,dim=c(1,Nv,length(matrici)))
  #alpha vettore
  alphatot=alphacumfinalef_y(matrici=matrici,storia=storia,Nv=Nv,k=k)
  #beta vettore
  betatot=betacumfinalenonalt_y(matrici,storia,Nv,k)
  
  #Qsi-1
  #doppio ciclo
  for (i in 1:nrow(storia)){ #per ogni stato #1:134
    if (storia$stato[i]==xsel){ #se xsel== stato i-esimo
      # print(i)
      # print("Lavora!")
      lowf=storia$inizio[i] #ti-1
      # print("lowf=")
      # print(lowf)
      upf=storia$fine[i] #ti
      alphaf=alphatot[,i] #alpha(i-1)
      matrice=matrici[paste0("matrice",storia$y[i])][[1]]
      Qsi=Qsif(iesima = xsel,matrice=matrice,Nv=Nv)
      betaf=betatot[,i+1] #beta(i) perche devo prendere betaw
      for (j in 1:Nv){
        dimensioni=k*(Nv)
        Deltajf=matrix(0,nrow=dimensioni,ncol=dimensioni)
        dim=Nv*(xsel-1)+j
        Deltajf[dim,dim]<-1 
        #Etj è un vettore!!
        Etj[,j,storia$y[i]]<-Etj[,j,storia$y[i]]+composite.simpson2(Qsip=Qsi,a=lowf,b=upf,alpha=alphaf,
                                            beta=betaf, Deltaj=Deltajf,n=10)
      }
    }
  }
  return(Etj) #torna un VETTORE, di lunghezza Nv, quello dell'integrale!
  
}
funzioneintegraleXselh_y(matrici =list(matrice1=Qsempl,matrice2=ssss) 
                         ,storia =asdas$storiafinale[,c(2,3,4,5)] ,xsel =2 ,Nv =2 ,k =2)
## computa la funzione per un singolo xj selezionato, integra su tutte 
## gli intervalli! Per ogni x (osservato), computiamo lo stesso tempo per ogni 
## Xh singolo

#computa il tempo totale passato in ogni stato!
funzioneETtotalesigma_y<-function(matrici,storia,Nv,k,N){
  ETime<-array(NA,dim=c(Nv,k,length(matrici)))
  for (t in 1:k){
    ETime[,t,]=funzioneintegraleXselh_y(matrici=matrici,storia=storia,xsel=t,Nv=Nv,k=k)
  }
  if(N!=(Nv)){ETime[(N+1):Nv,]<-0}
  
  #ETime=ETime*(storia$fine[nrow(storia)-storia$inizio[1]])/sum(ETime)
  return(ETime) #metto in colonna tutti i risultati di prima, ogni colonna un x diverso
}
##NON restituisce il valore normalizzato
ET_Y=funzioneETtotalesigma_y(matrici =list(matrice1=Qsempl,matrice2=ssss) 
                        ,storia =asdas$storiafinale[,c(2,3,4,5)] ,Nv =2 ,k =2 ,N =2 )

##  learning parametriMjk ------------------------------

#computa da alpha1- fino a alpha(tau-1)-
alpha_globale_meno_y<-function(matrici,storia,Nv,k ){
  alpha0<-rep(0,k*Nv)
  alpha0[(Nv*(storia$stato[1]-1)+1):(Nv*storia$stato[1])]<-1/Nv
  #creo array iniziale
  alphacumulato<-array(0.5,c(k*Nv,nrow(storia)))
  
  #Creo array di alphati dal quale origino le alphati-
  alphacumsenzameno<-alphacumfinalef_y(matrici,storia,Nv,k)
  #questo alpha ha in posizione 1 alpha0, posizione 2 alpha1 ecc...
  for (i in 2:(nrow(storia))){ #da alpha1 a alpha(tau=10)
    iesimainiz=storia$stato[i-1] #Si-2
    iesimafin=storia$stato[i] #Si-1
    # print(c(storia$stato[i-1],storia$stato[i]))
    # print("-----")
    #Qsi-2
    matrice=matrici[paste0("matrice",storia$y[i-1])][[1]]
    #matrice è con i-1 secondo me
    Qsi<-Qsif(iesima=iesimainiz,matrice=matrice,Nv=Nv)
    # Qsi1si0<-Qsi0si1f(iesima0 =iesimainiz ,iesima1 =iesimafin, matrice = matrice,Nv=Nv)
    ti1<-storia$fine[i-1] #ti-1
    # print(storia$fine[i+1])
    ti0<-storia$inizio[i-1] #ti-2
    # print(c(storia$inizio[i-1],storia$fine[i-1]))
    
    #computa in posizione i l'alpha(ti-1)- perche prende alpha (ti-2)*Qsi-2* ti-1; ti-2
    alphacumulato[,i]<-alphatpiu1meno(alphat =alphacumsenzameno[,i-1],Qsi = Qsi,ti1 = ti1
                                      ,ti0= ti0)
    
  } 
  return(alphacumulato[,-1]) #nella prima colonna avrei alpha(t=0)<-- Non mi serve, resta iniziale!
  #prima colonna non serve, resterebbe 0.5, il classico alpha0
  
  #### sicuri che la scelta migliore sia shiftare? perche resterebbe l'unica shiftata!!  
}

aaameno_y=alpha_globale_meno_y(matrici =list(matrice1=Qsempl,matrice2=ssss) 
                   ,storia =asdas$storiafinale[,c(2,3,4,5)],Nv =2 ,k =2 )


aaameno_y[,c(1:4,130:133)]

#COMPUTA DA B1+ a B(tau-1)+
beta_globale_piu_y<-function(matrici,storia,Nv,k ){
  ### funzione che computa il vettore beta(t)+
  
  #vettore betat senza +
  betasenzapiu<-betacumfinalenonalt_y(matrici=matrici,storia=storia,Nv=Nv,k=k)
  #beta(t-1)
  #in 1 ce beta0, .. in posizione i ce beta i-1
  betacumulato<-array(1.111,c(k*Nv,nrow(storia)+1))
  #inizializzo il mio vettore beta+
  
  for (i in ((nrow(storia))):2){ #9 a 2
    iesimainiz=storia$stato[i-1] #Si-2 #3
    iesimafin=storia$stato[i] #Si-1 #1 
    # print(c(storia$stato[i-1],storia$stato[i]))
    # print("-----")
    matrice=matrici[paste0("matrice",storia$y[i])][[1]]
    #Qsi-1
    Qsi<-Qsif(iesima=iesimafin,matrice=matrice,Nv=Nv) #da 2.7 a 3
    # Qsi1si0<-Qsi0si1f(iesima0 =iesimainiz ,iesima1 =iesimafin, matrice = matrice,Nv=Nv)
    #ti
    ti1<-storia$fine[i] #3
    # print(storia$fine[i+1])
    #ti-1
    ti0<-storia$inizio[i] #2.7
    # print(c(storia$inizio[i],storia$fine[i]))
    #sarebbe betatau, giusto 
    #beta(t)*Qsi-1*(ti-t(i-1))
    betacumulato[,i]<-betatifpiu(betat1 =betasenzapiu[,i+1],Qsi = Qsi,ti1 = ti1
                                 ,ti0= ti0 )
    ###nella posizioni i-esima sto computando beta(i-1)
    
    ##nb sembra tutto giusto ma bisognera riscalare i beta. in beta9 ci sara quello che è beta
  } 
  colnames(betacumulato)<-c(0,storia$fine)
  
  return(betacumulato[,-c(1,ncol(betacumulato))])
  #sta togliendo la prima e l'ultima colonna, che sarebbe l'11esima, che non ci serve
}

bbpiu_y=beta_globale_piu_y(matrici =list(matrice1=Qsempl,matrice2=ssss) 
                   ,storia =asdas$storiafinale[,c(2,3,4,5)],Nv =2 ,k =2)

bbpiu_y[,c(1:5,130:133)]

#computa MEjk di transizione da un sottostato all'altro!
singoloMEjktrans_y<-function(matrici,storia,Nv,k,j_me,k_me){
  ###SCELTA DRASTICA, SCELGO SOLO DI CONSIDERARE QUANDO 
  ###Y RESTA COSTANTE, ALTRIMENTI INCLUDEREI ANCHE PUNTI CON EVIDENZA COSTANTE
  #FORMULA: alphat- * Deltaj,k * betat+
  alfameno<-alpha_globale_meno_y(matrici=matrici,storia=storia,Nv=Nv,k=k)
  betapiu=beta_globale_piu_y(matrici,storia,Nv,k)
  MEjk=rep(0,length(matrici))
  Deltajk<-matrix(0,nrow=Nv*k,ncol=Nv*k)
  Deltajk[j_me,k_me]<-1
  for (i in 1:(nrow(storia)-1)){
    y=storia$y[i]
    if (y==storia$y[i+1]){
    #da alpha1, beta1 fino a alpha tau-1 beta tau-1
    MEjk[y]<-MEjk[y]+alfameno[,i]%*%Deltajk%*%betapiu[,i]
    # if(alfameno[,i]%*%Deltajk%*%betapiu[,i]!=0){
    #   print("storia")
    #   print(i)
    # }
    }
  }
  elementojk=rep(0,length(matrici))
  for (i in 1:length(matrici)){
    elementojk[i]=matrici[paste0("matrice",i)][[1]][j_me,k_me]
  }
  return(MEjk*elementojk)
}


#debugging
# alfameno<-alpha_globale_meno_y(matrici=list(matrice1=Qsempl,matrice2=ssss),
#                                storia=asdas$storiafinale[,c(2,3,4,5)],Nv=2,2)
# betapiu=beta_globale_piu_y(matrici=list(matrice1=Qsempl,matrice2=ssss),
#                            storia=asdas$storiafinale[,c(2,3,4,5)],Nv=2,2)
# MEjk=rep(0,2)
# Deltajk<-matrix(0,nrow=4,ncol=4)
# Deltajk[1,1]<-1
# for (i in 1:134){
#   # y=storia$y[i]
#   # if (y==storia$y[i+1]){
#     #da alpha1, beta1 fino a alpha tau-1 beta tau-1
#     MEjk[1]<-MEjk[1]+alfameno[,i]%*%Deltajk%*%betapiu[,i]
#     print(i)
#     print(alfameno[,i]%*%Deltajk%*%betapiu[,i])}
# 
# 
# alfameno[,48]%*%Deltajk%*%betapiu[,48]
# betapiu[betapiu<0]
# alfameno[alfameno<0]
# asdas$storiafinale[46:49,]
#da peso ase sei interno a penultimo stato oss!!
singoloMEjktrans_y(matrici =list(matrice1=Qsempl,matrice2=ssss) 
                   ,storia =asdas$storiafinale[,c(2,3,4,5)],Nv =2 ,k =2,1,2)
#computa ME NELLO STESSO SOTTOSISTEMA
funzioneintegraleMEjk_y<-function(matrici,storia,Nv,k,j_me,k_me){
  # FUnzione che computa la funzione MEjk in intervalli di conoscenza costante
  ###
  Etj=rep(0,length(matrici))
  alphatot=alphacumfinalef_y(matrici=matrici,storia=storia,Nv=Nv,k=k)
  betatot=betacumfinalenonalt_y(matrici,storia,Nv=Nv,k=k)
  
  # Qsi=Qsif(iesima = xsel,matrice=matrice,Nv=Nv)
  for (i in 1:(nrow(storia))){
    # if (storia$stato[i]==xsel){
    # print(i)
    # print("Lavora!")
    matrice=matrici[paste0("matrice",storia$y[i])][[1]]
    y=storia$y[i]
    QS<-Qsif(iesima = storia$stato[i],matrice=matrice,Nv=Nv) #So
    lowf=storia$inizio[i] #0
    #print("lowf=")
    # print(lowf)
    upf=storia$fine[i] #0.5
    #print(upf)
    alphaf=alphatot[,i] #alpha_i-1
    betaf=betatot[,i+1] #beta_i
    #perche devo prendere betaw
    #beta1
    dimensioni=k*(Nv)
    Deltajf=matrix(0,nrow=dimensioni,ncol=dimensioni)
    # dim=Nv*(xsel-1)+j
    Deltajf[j_me,k_me]<-1
    # if(integralinosingoloET(Qsip=Qsi,low=lowf,up=upf,alpha=alphaf,
    #                         beta=betaf, Deltaj=Deltajf)>0){ 
    # print(integralinosingoloET(Qsip=Qsi,low=lowf,up=upf,alpha=alphaf,
    #                            beta=betaf, Deltaj=Deltajf))
    # print("")
    # print(i*10)
    # print(j_me)
    # print(k_me)
    # print("")
    # }
    Etj[y]<-Etj[y]+composite.simpson2(Qsip=QS,a=lowf,b=upf,alpha=alphaf,
                                beta=betaf, Deltaj=Deltajf,n=10)
  }
  #ritorna qjk*E[ME]
  elementojk=rep(0,length(matrici))
  for (i in 1:length(matrici)){
    elementojk[i]=matrici[paste0("matrice",i)][[1]][j_me,k_me]
  }
  return(elementojk*Etj)
  
}
funzioneintegraleMEjk_y(matrici =list(matrice1=Qsempl,matrice2=ssss) 
                      ,storia =asdas$storiafinale[,c(2,3,4,5)],
                      Nv =2 ,k =2 ,j_me =3 ,k_me = 1)

#dovrebbe funzionare

#potrei farlo piu veloce!
funzMEtotale_y<-function(matrici,storia,Nv,k){ 
  MEtotale2=array(0,dim=c(k*Nv,k*Nv,length(matrici)))
  for( el in 1:(k*Nv)){
    for (el2 in 1:(k*Nv)){
      k2=el%/%Nv
      #print("")
      #if(!el2 %in% seq(k2*Nv+1,(k2+1)*Nv)){
      if(matrici$matrice1[el,el2]!=0){
        
        MEtotale2[el,el2,]=MEtotale2[el,el2,]+singoloMEjktrans_y(matrici=matrici,storia=storia
                                                             ,Nv=Nv,k=k,j_me=el,k_me=el2)
        MEtotale2[el,el2,]=MEtotale2[el,el2,]+funzioneintegraleMEjk_y(matrici=matrici,storia=storia,Nv=Nv,k=k,
                                                                  j_me=el,k_me=el2)
      }
      #se il primo numero è divisibile per k
      ## e il secondo numero è diverso da el-1 e el-2
      
      #}
      
      print(paste0("iterazione numero: ",el," ", el2 ))
    }
  }
  return(MEtotale2)
}

funzMEtotale_y(matrici =list(matrice1=Qsempl,matrice2=ssss) 
               ,storia =asdas$storiafinale[,c(2,3,4,5)],Nv = 2,k = 2)
#dovrebbe funzionare

# DEVO DARE UNA CONDIZIONE PER FARLO PIU VELOCE!!! ------------------------
# attenzione alla condizione per ME ---------------------------------------
##il !


####   funzionegeneratricestoria------
#genera una storia da una certa matrice
generastorie_y<-function(matrice1,matrice2,matrice_y,tmax,statix,Nv,N){
  #### un problemino da dire, quando scatta h o x y resta sempre uguale
  #### quando invece scatta y anche h o x si azzera!
  
  statixh=seq(from=1,to=nrow(matrice1),1)
  statiy=seq(from=1,to=nrow(matrice_y),1)
  x <- sample(x =c(1:nrow(matrice1)) ,size =1 ,replace = T,prob = rep(1/nrow(matrice1),nrow(matrice1)))
  y <- sample(x =c(1:nrow(matrice_y)) ,size =1 ,replace = T,prob = rep(1/nrow(matrice_y),nrow(matrice_y)))
  
  collezionestorie<-data.frame("statoxh"=NA,"inizio"=NA,"fine"=NA,"y"=NA)
  i=1
  tempo=0
  
  
  Qsecondo_y=matrice_y
  diag(Qsecondo_y)=0
  
  repeat {
    if (y==1) matrice=matrice1 else{matrice=matrice2}
    Qsecondo=matrice
    diag(Qsecondo)=0
    c1=x
    c2=tempo
    c_y=y
    tempo_i=rexp(1,rate=abs(matrice[x,x]))
    tempo_y=rexp(1,rate=abs(matrice_y[y,y]))
    if (tempo_i< tempo_y){
      repeat{
        tempo_y<-tempo_y-tempo_i
        if (y==1) matrice=matrice1 else{matrice=matrice2}
        tempo=tempo+tempo_i
        c3=tempo
        collezionestorie<-rbind(collezionestorie,c(c1,c2,c3,c_y))
        x<-sample(x =c(1:nrow(matrice)) ,size =1 ,replace = T, prob = c(Qsecondo[x,]))
        c1=x
        c2=tempo
        i=i+1
        tempo_i=rexp(1,rate=abs(matrice[x,x]))
        if(tempo_i> tempo_y | tempo >= tmax){
          break}
        
      }      
      
    }
    
    if (tempo_y<=tempo_i & tempo<tmax){
      tempo=tempo+tempo_y
      c3=tempo
      collezionestorie<-rbind(collezionestorie,c(c1,c2,c3,c_y))
      y<-sample(x =c(1:nrow(matrice_y)) ,size =1 ,replace = T, prob = c(Qsecondo_y[y,]))
      c_y=y
      c2=tempo
      
    }
    
    if (tempo >= tmax){
      break
    }
  }
  
  collezionestorie$fine[nrow(collezionestorie)]<-tmax
  stat=((collezionestorie$statoxh-1)%/%Nv)+1
  collezionestorie$stato<-stat
  collezionestorie<-collezionestorie[-1,]
  i2=2
  storiaqsempl<-collezionestorie
  
  for (i in 2:nrow(storiaqsempl)){
    
    if(collezionestorie$stato[i]!=collezionestorie$stato[i-1] | 
       collezionestorie$y[i]!=collezionestorie$y[i-1]){
      storiaqsempl[i2,]<-collezionestorie[i,]
      i2=i2+1
    }}
  storiaqsempl
  storiaqsempl<-storiaqsempl[1:(i2-1),]
  storiaqsempl$fine<-c(storiaqsempl$inizio[-1],tmax)
  storiaqsempl
  collezionestorie
  return(list(storiafinale=storiaqsempl,storiavera=collezionestorie))
}

matrice1=Qsempl
matrice2=ssss
matrice_y=matrix(c(-10.1,10.1,10.2,-10.2),nrow=2,byrow = T)
tQsempl=generastorie(matrice =Qsempl ,tmax =100 ,statix =2 ,Nv =2 ,N = 2)
(ssss=matrix(c(-4.5,1.5,3,0,1,-6,0,5,2,0,-5.5,3.5,0,1,0.5,-1.5
),nrow=4,byrow = T))
matrice_y=
asdas=generastorie_y(matrice1 = Qsempl,matrice2 = ssss,matrice_y = matrice_y,
                     tmax = 100,statix = 2,Nv=2,N=2)

collezionestorie=asdas$storiavera
round(asdas$storiafinale[c(1:30),],2) 
round(asdas$storiafinale[c(1150:1157),],5)
round(asdas$storiavera[c(1:30),],2)

matrice_y=matrix(c(-0.1,0.1,0.2,-0.2),nrow=2,byrow = T)
statix =2 ;Nv =2 ;N = 2
#statix= quanti stati ha la x?

##funzione che fa uno step di expmax con pero 1! sigmaoss

######EMSTEP #Da cambiare la normalizzazione #ci ho provato!
EMstep_y<-function(matricipartenza,storiaoss,Nv,k,N){
  # print("aa")
  EEsemplice=funzioneETtotalesigma_y(matrici=matricipartenza,storia=storiaoss,Nv=Nv,k=k,N=N)
  # print("ciao")
  metotalesemplice=funzMEtotale_y(matrici=matricipartenza,storia=storiaoss,Nv=Nv,k=k)
  # print("ciao2")
  
  
  ###modifiche!!!
  for(i in 1:length(matricipartenza)){
    sm[i]=sum(EEsemplice[,,i])
    tempotot[i]=sum(storiaoss$fine[storiaoss$y==i]-storiaoss$inizio[storiaoss$y==i])
    EEsemplice[,,i]=(EEsemplice[,,i]/sm[i])*tempotot[i]
    metotalesemplice[,,i]=(metotalesemplice[,,i]/sm[i])*tempotot[i]
  }

  ###fine modifiche
  
  
  Oggettodaritornare= vector(mode = "list", length = length(matricipartenza))
  
  i=0
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
  nuovamatrice=matrix(0,nrow=nrow(matricepartenza),ncol=nrow(matricepartenza))
  diag(nuovamatrice)=qx
  for (k in 1:nrow(nuovamatrice)){
    for (j in 1:ncol(nuovamatrice)){
      if (k != j){
        nuovamatrice[k,j]<-thetaxx[k,j]*qx[k]
      }
    }
  }
  diag(nuovamatrice)=-diag(nuovamatrice)
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
  # [1] -1.024622
  Oggettodaritornare[[i]]<-list(TotET=TotET,LogVero=LogVero,LogVerononNorm=LogVerononNorm,
                           ET=EEsemplice[,,i],Mejk=metotalesemplice[,,i],matricefinale=nuovamatrice)
}
  return(Oggettodaritornare)
  
}

start=Sys.time()
eee_y2=EMstep_y(matricipartenza =list(matrice1=Qsempl,matrice2=ssss) 
         ,storiaoss =asdas$storiafinale[,c(2,3,4,5)],Nv =2 ,k =2 ,N =2)
end=Sys.time()
end-start
eee_y[[1]]$LogVero
eee_y2[[2]]$LogVero
eeesecondogiro[[1]]$LogVero
eeesecondogiro[[2]]$LogVero
eeeterzogiro[[1]]$LogVero
eeeterzogiro[[2]]$LogVero
eeesecondogiro=EMstep_y(matricipartenza =list(matrice1=eee_y[[1]]$matricefinale,matrice2=
                                 eee_y[[2]]$matricefinale) 
         ,storiaoss =asdas$storiafinale[,c(2,3,4,5)],Nv =2 ,k =2 ,N =2)


eeeterzogiro=EMstep_y(matricipartenza =list(matrice1=eeesecondogiro[[1]]$matricefinale,matrice2=
                                              eeesecondogiro[[2]]$matricefinale) 
                      ,storiaoss =asdas$storiafinale[,c(2,3,4,5)],Nv =2 ,k =2 ,N =2)
##come sopra ma prende in input una LISTA di storie!
#possibile aggiunta, per ora considera le storie tutte stessa lunghezza
#dovrei cambiarla come cosa!
EMstepdoppio<-function(matricepartenza,storie,Nv,k,N){
  
  #inizializzazione
  EEsemplicetot=matrix(0,nrow=Nv,ncol=k)
  metotale=matrix(0,nrow=k*Nv,ncol=k*Nv)
  
  
  for (storiaoss in storie){
    #per ogni storia, ogni traiettoria
    
    # T[j]
    EEsemplice=funzioneETtotalesigma(matrice=matricepartenza,storia=storiaoss,Nv=Nv,k=k,N=N)
    
    #le due somme del tempo
    sm=sum(EEsemplice)
    tempotot=storiaoss$fine[nrow(storiaoss)]-storiaoss$inizio[1]
    #normalizziamo
    EEsemplice=(EEsemplice/sm)*tempotot
    
    #ME[j,k] e norm
    metotalesemplice=funzMEtotale(matrice=matricepartenza,storia=storiaoss,Nv=Nv,k=k)
    metotalesemplice=(metotalesemplice/sm)*tempotot
    
    #somma
    EEsemplicetot=EEsemplicetot+EEsemplice
    metotale=metotale+metotalesemplice
  }
  qiniziale=-diag(matricepartenza)
  thetainiz=matricepartenza/qiniziale
  mpr<-metotale
  diag(mpr)<-0
  Mx<-rowSums(mpr) 
  EExsemplice<-as.vector(EEsemplicetot)
  qx<-Mx/EExsemplice
  thetaxx<-metotale/Mx 
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
  for (i in 1:nrow(metotale)){
    for (j in 1:ncol(metotale)){
      if (i !=j){
        temp=metotale[i,j]*log(thetainiz[i,j])
        if (temp !="NaN"){ 
          
          provaperlogalt=provaperlogalt+metotale[i,j]*log(thetainiz[i,j])
        }
      }}
  }
  LogVerononNorm=sum(Mx*log(qiniziale)-qiniziale*as.vector(EEsemplicetot))+provaperlogalt
  #[1] [1]--5.751853e-08
  
  TotET=sum(EEsemplicetot) #5.613631e-08
  LogVero=LogVerononNorm/sum(EEsemplicetot)
  # [1] -1.024622
  Oggettodaritornare<-list(TotET=TotET,LogVero=LogVero,LogVerononNorm=LogVerononNorm,
                           ET=EEsemplicetot,Mejk=metotale,matricefinale=nuovamatrice)
  return(Oggettodaritornare)
  
}
#versione migliroe con timetot


EMstepdoppio_y<-function(matricipartenza,storie,Nv,k,N){
  ##errore errore!!! usa solo la prima storia!!
  #fa uno step di exp max con un parent y
  EEtot=array(0,dim=c(Nv,k,length(matricipartenza))) #rray tridimensionale vuoto
  metotale=array(0,dim=c(Nv*k,Nv*k,length(matricipartenza)))
  #terza dimensioni, quante matrici di partenza, tante quanti gli y
  for (storiaoss in storie){ #per ogni raiettoria osservata
    # storiaoss=storie$storia1 #che coglionata
    
    EEsemplice=funzioneETtotalesigma_y(matrici=matricipartenza,storia=storiaoss,Nv=Nv,k=k,N=N)
    #dim(EEsemplice) #[1] 3 4 2
    #Calcola il tempo medio trascorso in ogni stato in quella storia
    
    mesemplice=funzMEtotale_y(matrici=matricipartenza,storia=storiaoss,Nv=Nv,k=k)
    #12*12*2
    #calcola transizione medie da ogni stato ad ogni altro. 
    
    #inizializzazioni vettori per normalizzare ME[j,k] e T[j]
    sm=vector()
    tempotot=vector()
    for (i in 1:length(matricipartenza)){ #per ogni storia, per ogni matrice (ogni y)
      
      #tempo trascorso in y==1, y==2,...
      tempotot[i]=sum(storiaoss$fine[storiaoss$y==i]-storiaoss$inizio[storiaoss$y==i])
      #somma per normalizzare a 1 T[j| y==1]
      sm[i]=sum(EEsemplice[,,i]) 
      
      #normalizzazioni
      EEsemplice[,,i]=(EEsemplice[,,i]/sm[i])*tempotot[i]
      mesemplice[,,i]=(mesemplice[,,i]/sm[i])*tempotot[i]
      
      #sommo l'iesimo contributo al totale. 
      EEtot[,,i]=EEtot[,,i]+EEsemplice[,,i]
      metotale[,,i]=metotale[,,i]+mesemplice[,,i]
    }
    
  }

  #inizializzazione di i, per ogni matrice (ogni istanza di y==1,..)
  i=0
  Oggettodaritornare= vector(mode = "list", length = length(matricipartenza))
  
  for (matricepartenza in matricipartenza){
    #matrice finale dato un certo y==istanza
    i=i+1
    qiniziale=-diag(matricepartenza) #iniziale
    thetainiz=matricepartenza/qiniziale #thetainiz
    
    
    mpr<-metotale[,,i] 
    diag(mpr)<-0
    
    #somma delle righe, il me totale
    Mx<-rowSums(mpr) 
    EExsemplice<-as.vector(EEtot[,,i])
    
    #definizione di qx e thetax
    qx<-Mx/EExsemplice
    thetaxx<-metotale[,,i]/Mx 
    
    #inizializzazione matrice
    nuovamatrice=matrix(0,nrow=nrow(matricepartenza),ncol=nrow(matricepartenza))
    
    diag(nuovamatrice)=qx #tutte le intensità USCENTI dallo stato
    
    for (k in 1:nrow(nuovamatrice)){
      for (j in 1:ncol(nuovamatrice)){
        if (k != j){
          nuovamatrice[k,j]<-thetaxx[k,j]*qx[k] #do i valori ad ogni intensita k->j
        }
      }
    }
    diag(nuovamatrice)=-diag(nuovamatrice) #da qjj--> -qjj
    
  #computazione loglikelihood
    
    #inizializzazion  
    provaperlogalt=0
    for (l in 1:nrow(metotale)){
      for (j in 1:ncol(metotale)){
        if (l !=j){
          temp=metotale[l,j,i]*log(thetainiz[l,j])
          if (temp !="NaN"){ 
            
            provaperlogalt=provaperlogalt+metotale[l,j,i]*log(thetainiz[l,j])
          }
        }}
    }
    #computazione logvero
    
    LogVerononNorm=sum(Mx*log(qiniziale)-qiniziale*as.vector(EEtot[,,i]))+provaperlogalt
    #[1] [1]--5.751853e-08
    
    TotET=sum(EEtot[,,i]) #5.613631e-08
    LogVero=LogVerononNorm/sum(EEtot[,,i])
    # [1] -1.024622
    Oggettodaritornare[[i]]<-list(TotET=TotET,LogVero=LogVero,LogVerononNorm=LogVerononNorm,
                                  ET=EEtot[,,i],Mejk=metotale[,,i],matricefinale=nuovamatrice)
  }
  
  return(Oggettodaritornare)
}
##funzione che converti la storia da y a normale 
fconvertistoria=function(storia){
  ss=storia
  i2=0
  for(i in 2:nrow(storia)){
    
    if (storia$y[i]==storia$y[i-1]){
      i2=i2+1
      ss[i2,]=storia[i,]
    }
    
  }
  ss2=rbind(storia[1,],ss)
  for(i in 1:(nrow(ss2)-1)){
    ss2$fine[i]=ss2$inizio[i+1]
  }
  ss2$fine[i2+1]=storia$fine[nrow(storia)]
  return(ss2[(1:(i2+1)),c(1,2,4)])
  
}


# problema: se storia troppo lunga alpha e beta diventano zero ------------

##tuttp a puttane


# problema: MEjk trans DA PESO ANCHE INTRA TRANSIZIONE!! ------------------
##daun problema in MEjk[j==penultimo_stato_oss,k==appartenti a quello stato]
#Es: j==3; k==1 
#oppure j==1; k==3
##i beta+ sono visibilmente errati! pero non danno problemi, tolti quelli!


##anche integraleMEjk stesso problema, con ULTIMO STATO OSSERVATO
#PROBLEMA IN STATO FINALE DI ULTIMO STATO OSSERVATO
#ES J==6; K=1 o K=7; o K=10
#dovrebbero essere zero ma in realta sono diversi! 
#qui pero i beta e gli alfa sono corretti!! cio che è sbaglaito è concettualmente!

