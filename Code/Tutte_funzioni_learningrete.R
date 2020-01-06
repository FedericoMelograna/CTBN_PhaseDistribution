
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

alphatpiu1meno<-function(alphat,Qsi,ti1,ti0){
  alphat%*%expm(Qsi*(ti1-ti0))
}

betatifpiu<-function(betat1,Qsi,ti1,ti0){ #P[sigmat:tau | Xt+=i]
  betat0<-expm(Qsi*(ti1-ti0))%*%betat1
  return(betat0)
}

betacumulatofn<-function(matrice,storia,Nv,k ){
  #restituisce un vettore di beta, da beta1 a beta_tau  
  betacumulato<-array(1,c(k*Nv,nrow(storia)+1))
  #vettore di partenza
  #da 9 a 2
  for (i in (nrow(storia)):2){ #ciclo al contrario, per mettere dentro gli elementi ricorsivamente
    iesimainiz=storia$stato[i-1] #si-1
    iesimafin=storia$stato[i] #si
    # print(c(storia$stato[i-1],storia$stato[i]))
    # print("-----")
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
  return(betacumulato)
}
betacumulatof0<-function(matrice,storia,betacumulprimo,Nv,k){ 
  #crea beta0 da sovrascrivere poi
  beta0<-array(1,c(k*Nv,1))
  iesima=storia$stato[1] #sia fin che iniz
  #S0
  # # iesimafin=storia$stato[i+1]
  # print(storia$stato[1])
  # print("-----")
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
betacumfinalenonalt<-function(matrice,storia,Nv,k){
  #crea vettore beta giusto
  bbb=betacumulatofn(matrice=matrice,storia=storia,Nv=Nv,k=k)
  bbb0=betacumulatof0(matrice=matrice,storia=storia,betacumulprimo = bbb[,2],Nv=Nv,k=k)
  bbb[,1]<-bbb0
  bbbfin<-bbb
  return(bbbfin)
  
}
#da beta0 a beta_tau;####versione non alternativa cumulata

#alpha
alphacumulatofn<-function(matrice,storia,Nv,k ){
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
    iesimainiz=storia$stato[i-1] #Si-2
    iesimafin=storia$stato[i] #Si-1
    # print(c(storia$stato[i-1],storia$stato[i]))
    # print("-----")
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
  return(alphacumulato)
  #dimensioni(1,10)---> le mappero in 0, tau=9
}
alphacumulftau<-function(matrice,storia,alphacumultaumenouno,Nv,k){ 
  #alphatau; #poi lo sostituirei in alphacumulato[10]
  alphatau<-array(1,c(k*Nv,1)) #vettore iniziale
  iesima=storia$stato[nrow(storia)] #sia fin che iniz
  ###ultimo stato, stato S(tau-1)
  
  # # iesimafin=storia$stato[i+1]
  # print(storia$stato[nrow(storia)])
  # print("-----")
  
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
alphacumfinalef<-function(matrice,storia,Nv,k){
  aaa=alphacumulatofn(matrice=matrice,storia=storia,Nv=Nv,k=k)
  aaa0=alphacumulftau(matrice=matrice,storia=storia,alphacumultaumenouno = aaa[,ncol(aaa)-1],Nv=Nv,k=k)
  #sostituisci ad alpha_tau<-- la sua rappresentazione con alphatau
  aaa[,ncol(aaa)]<-aaa0
  colnames(aaa)<-c(0,storia$fine)
  aaafin<-aaa
  return(aaafin)
  
}
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


#integra per un certo x selezionato su tutta la storai
funzioneintegraleXselh<-function(matrice,storia,xsel,Nv,k){
  # xsel=1; 
  Etj=matrix(0,nrow=1,ncol=Nv)
  #alpha vettore
  alphatot=alphacumfinalef(matrice=matrice,storia=storia,Nv=Nv,k=k)
  #beta vettore
  betatot=betacumulatofn(matrice,storia,Nv,k)
  
  #Qsi-1
  Qsi=Qsif(iesima = xsel,matrice=matrice,Nv=Nv)
  #doppio ciclo
  for (i in 1:nrow(storia)){ #per ogni stato
    if (storia$stato[i]==xsel){ #se xsel== stato i-esimo
      # print(i)
      # print("Lavora!")
      lowf=storia$inizio[i] #ti-1
      # print("lowf=")
      # print(lowf)
      upf=storia$fine[i] #ti
      alphaf=alphatot[,i] #alpha(i-1)
      betaf=betatot[,i+1] #beta(i) perche devo prendere betaw
      for (j in 1:Nv){
        dimensioni=k*(Nv)
        Deltajf=matrix(0,nrow=dimensioni,ncol=dimensioni)
        dim=Nv*(xsel-1)+j
        Deltajf[dim,dim]<-1 
        #Etj è un vettore!!
        Etj[,j]<-Etj[,j]+composite.simpson2(Qsip=Qsi,a=lowf,b=upf,alpha=alphaf,
                                            beta=betaf, Deltaj=Deltajf,n=10)
      }
    }
  }
  return(Etj) #torna un VETTORE, di lunghezza Nv, quello dell'integrale!
  
}
## computa la funzione per un singolo xj selezionato, integra su tutte 
## gli intervalli! Per ogni x (osservato), computiamo lo stesso tempo per ogni 
## Xh singolo

#computa il tempo totale passato in ogni stato!
funzioneETtotalesigma<-function(matrice,storia,Nv,k,N){
  ETime<-matrix(NA,nrow=Nv,ncol=k)
  for (t in 1:k){
    ETime[,t]=funzioneintegraleXselh(matrice=matrice,storia=storia,xsel=t,Nv=Nv,k=k)
  }
  if(N!=(Nv)){ETime[(N+1):Nv,]<-0}
  
  #ETime=ETime*(storia$fine[nrow(storia)-storia$inizio[1]])/sum(ETime)
  return(ETime) #metto in colonna tutti i risultati di prima, ogni colonna un x diverso
}
##NON restituisce il valore normalizzato


##  learning parametriMjk ------------------------------

#computa da alpha1- fino a alpha(tau-1)-
alpha_globale_meno<-function(matrice,storia,Nv,k ){
  alpha0<-rep(0,k*Nv)
  alpha0[(Nv*(storia$stato[1]-1)+1):(Nv*storia$stato[1])]<-1/Nv
  #creo array iniziale
  alphacumulato<-array(0.5,c(k*Nv,nrow(storia)))
  
  #Creo array di alphati dal quale origino le alphati-
  alphacumsenzameno<-alphacumfinalef(matrice,storia,Nv,k)
  #questo alpha ha in posizione 1 alpha0, posizione 2 alpha1 ecc...
  for (i in 2:(nrow(storia))){ #da alpha1 a alpha(tau=10)
    iesimainiz=storia$stato[i-1] #Si-2
    iesimafin=storia$stato[i] #Si-1
    # print(c(storia$stato[i-1],storia$stato[i]))
    # print("-----")
    #Qsi-2
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


#COMPUTA DA B1+ a B(tau-1)+
beta_globale_piu<-function(matrice,storia,Nv,k ){
  ### funzione che computa il vettore beta(t)+
  
  #vettore betat senza +
  betasenzapiu<-betacumfinalenonalt(matrice=matrice,storia=storia,Nv=Nv,k=k)
  #beta(t-1)
  #in 1 ce beta0, .. in posizione i ce beta i-1
  betacumulato<-array(1.111,c(k*Nv,nrow(storia)+1))
  #inizializzo il mio vettore beta+
  
  for (i in ((nrow(storia))):2){ #9 a 2
    iesimainiz=storia$stato[i-1] #Si-2 #3
    iesimafin=storia$stato[i] #Si-1 #1 
    # print(c(storia$stato[i-1],storia$stato[i]))
    # print("-----")
    
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

#computa MEjk di transizione da un sottostato all'altro!
singoloMEjktrans<-function(matrice,storia,Nv,k,j_me,k_me){
  #FORMULA: alphat- * Deltaj,k * betat+
  alfameno<-alpha_globale_meno(matrice=matrice,storia=storia,Nv=Nv,k=k)
  betapiu=beta_globale_piu(matrice,storia,Nv,k)
  MEjk=0
  Deltajk<-matrix(0,nrow=Nv*k,ncol=Nv*k)
  Deltajk[j_me,k_me]<-1
  for (i in 1:(nrow(storia)-1)){
    MEjk<-MEjk+alfameno[,i]%*%Deltajk%*%betapiu[,i]
    # if(alfameno[,i]%*%Deltajk%*%betapiu[,i]!=0){
    #   print("storia")
    #   print(i)
    # }
  }
  
  return(MEjk*matrice[j_me,k_me])
}

#computa ME NELLO STESSO SOTTOSISTEMA
funzioneintegraleMEjk<-function(matrice,storia,Nv,k,j_me,k_me){
  # FUnzione che computa la funzione MEjk in intervalli di conoscenza costante
  ###
  Etj=0
  alphatot=alphacumfinalef(matrice=matrice,storia=storia,Nv=Nv,k=k)
  betatot=betacumfinalenonalt(matrice,storia,Nv=Nv,k=k)
  
  # Qsi=Qsif(iesima = xsel,matrice=matrice,Nv=Nv)
  for (i in 1:(nrow(storia))){
    # if (storia$stato[i]==xsel){
    # print(i)
    # print("Lavora!")
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
    Etj<-Etj+composite.simpson2(Qsip=QS,a=lowf,b=upf,alpha=alphaf,
                                beta=betaf, Deltaj=Deltajf,n=10)
  }
  #ritorna qjk*E[ME]
  return(matrice[j_me,k_me]*Etj)
  
}



#potrei farlo piu veloce!
funzMEtotale<-function(matrice,storia,Nv,k){ 
  MEtotale2=matrix(0,nrow=k*Nv,ncol=k*Nv)
  for( el in 1:(k*Nv)){
    for (el2 in 1:(k*Nv)){
      k2=el%/%Nv
      #print("")
      #if(!el2 %in% seq(k2*Nv+1,(k2+1)*Nv)){
      MEtotale2[el,el2]=funzioneintegraleMEjk(matrice=matrice,storia=storia,Nv=Nv,k=k,
                                              j_me=el,k_me=el2)
      #}
      MEtotale2[el,el2]=MEtotale2[el,el2]+singoloMEjktrans(matrice=matrice,storia=storia
                                                           ,Nv=Nv,k=k,j_me=el,k_me=el2)
      
      print(paste0("iterazione numero: ",el," ", el2 ))
    }
  }
  return(MEtotale2)
}

# DEVO DARE UNA CONDIZIONE PER FARLO PIU VELOCE!!! ------------------------
funzMEtotale<-function(matrice,storia,Nv,k){ 
  MEtotale2=matrix(0,nrow=k*Nv,ncol=k*Nv)
  for( el in 1:(k*Nv)){
    for (el2 in 1:(k*Nv)){
      k2=el%/%Nv
      #print("")
      #if(!el2 %in% seq(k2*Nv+1,(k2+1)*Nv)){
      if(matrice[el,el2]!=0){
        
        MEtotale2[el,el2]=MEtotale2[el,el2]+singoloMEjktrans(matrice=matrice,storia=storia
                                                             ,Nv=Nv,k=k,j_me=el,k_me=el2)
        MEtotale2[el,el2]=MEtotale2[el,el2]+funzioneintegraleMEjk(matrice=matrice,storia=storia,Nv=Nv,k=k,
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
#cosi è piu veloce
# attenzione alla condizione per ME ---------------------------------------
##il !


####   funzionegeneratricestoria------
#genera una storia da una certa matrice
generastorie<-function(matrice,tmax,statix,Nv,N){
  statixh=seq(from=1,to=nrow(matrice),1)
  x <- sample(x =c(1:nrow(matrice)) ,size =1 ,replace = T,prob = rep(1/nrow(matrice),nrow(matrice)))
  collezionestorie<-data.frame("statoxh"=NA,"inizio"=NA,"fine"=NA)
  i=1
  tempo=0
  Qsecondo=matrice
  diag(Qsecondo)=0
  repeat {
    c1=x
    c2=tempo
    tempo_i=rexp(1,rate=abs(matrice[x,x]))
    x<-sample(x =c(1:nrow(matrice)) ,size =1 ,replace = T, prob = c(Qsecondo[x,]))
    tempo=tempo+tempo_i
    c3=tempo
    collezionestorie<-rbind(collezionestorie,c(c1,c2,c3))
    if (tempo >= tmax){
      break
    }
    i=i+1
  }
  
  collezionestorie$fine[nrow(collezionestorie)]<-tmax
  stat=((collezionestorie$statoxh-1)%/%Nv)+1
  collezionestorie$stato<-stat
  collezionestorie<-collezionestorie[-1,]
  i2=2
  storiaqsempl<-collezionestorie
  for (i in 2:nrow(storiaqsempl)){
    if(collezionestorie$stato[i]!=collezionestorie$stato[i-1]){
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
#statix= quanti stati ha la x?

##funzione che fa uno step di expmax con pero 1! sigmaoss

######EMSTEP
EMstep<-function(matricepartenza,storiaoss,Nv,k,N){
  # print("aa")
  EEsemplice=funzioneETtotalesigma(matrice=matricepartenza,storia=storiaoss,Nv=Nv,k=k,N=N)
  # print("ciao")
  metotalesemplice=funzMEtotale(matrice=matricepartenza,storia=storiaoss,Nv=Nv,k=k)
  # print("ciao2")
  
  ###Modifiche!
  sm=sum(EEsemplice)
  tempotot=storiaoss$fine[nrow(storiaoss)]-storiaoss$inizio[1]
  EEsemplice=(EEsemplice/sm)*tempotot
  metotalesemplice=(metotalesemplice/sm)*tempotot
  ##fine modifiche
  
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


##come sopra ma prende in input una LISTA di storie!
#possibile aggiunta, per ora considera le storie tutte stessa lunghezza
#dovrei cambiarla come cosa!
EMstepdoppio<-function(matricepartenza,storie,Nv,k,N){
  EEsemplicetot=matrix(0,nrow=Nv,ncol=k)
  metotale=matrix(0,nrow=k*Nv,ncol=k*Nv)
  for (storiaoss in storie){
    
    EEsemplice=funzioneETtotalesigma(matrice=matricepartenza,storia=storiaoss,Nv=Nv,k=k,N=N)
    sm=sum(EEsemplice)
    EEsemplice=EEsemplice/sm
    metotalesemplice=funzMEtotale(matrice=matricepartenza,storia=storiaoss,Nv=Nv,k=k)
    metotalesemplice=metotalesemplice/sm
    EEsemplicetot=EEsemplicetot+EEsemplice
    metotale=metotale+metotalesemplice
  }
  metotale=metotale/length(storie)
  EEsemplicetot=EEsemplicetot/length(storie)
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

#versione migliorata con anche il tempotot moltiplicato
EMstepdoppio<-function(matricepartenza,storie,Nv,k,N){
  EEsemplicetot=matrix(0,nrow=Nv,ncol=k)
  metotale=matrix(0,nrow=k*Nv,ncol=k*Nv)
  for (storiaoss in storie){
    
    EEsemplice=funzioneETtotalesigma(matrice=matricepartenza,storia=storiaoss,Nv=Nv,k=k,N=N)
    sm=sum(EEsemplice)
    tempotot=storiaoss$fine[nrow(storiaoss)]-storiaoss$inizio[1]
    EEsemplice=(EEsemplice/sm)*tempotot
    metotalesemplice=funzMEtotale(matrice=matricepartenza,storia=storiaoss,Nv=Nv,k=k)
    metotalesemplice=(metotalesemplice/sm)*tempotot
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

