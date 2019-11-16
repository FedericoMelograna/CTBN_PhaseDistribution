#########proviamo a imparareuna rete bayesiana
k=4
Nv=7
df=df+0.01
#storie di prova
# df2 è df con max 100
sigmaosservato1=data.frame(stato=c(1,3,2,1,3,4,1,3,1),inizio=c(0,0.5,0.9,1.3,1.5,1.9,2.2,2.5,2.7),
                           fine=c(0.5,0.9,1.3,1.5,1.9,2.2,2.5,2.7,3))

#f(sigmaosservato1)
sigmaosservato2=data.frame(stato=c(1,2,1,2,1,2,1,2,1),inizio=c(0,0.4,0.9,1.1,1.5,2,2.2,2.5,2.7),
                           fine=c(0.4,0.9,1.1,1.5,2,2.2,2.5,2.7,3))

storia=sigmaosservato1
#matrice di prova
m=matrix(c(-1,1,2,-2),nrow=2,byrow = T)
#del tutto inutile
library("expm", lib.loc="~/R/win-library/3.4")
#libreria per esponenzializzazione matrice

# computare alphax e betat ------------------------------------------------

#inizializzo a caso beta con la probabilità di trovarsi negli stadi. 
##con #di stadi pari al numero totale di stadi 
# (K*(n+v))
#supponiamo passa da stato 1 a 2 a 0.7 
df #matrice del prof
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
#restituisce una matrice QS dove ha mandato a zero tutte le entrate diverse dallo stato i-esimo
Qsiprova=Qsif(3,matrice=df,Nv=7)#[1:15,1:15]

Qsip=Qsif(1,matrice=df,Nv=7)
##all'interno del sistema. 
Qsiprova[15:21,15:21]
Qsip[1:7,1:7]
df[1:7,1:7]

# PROBLEMA;: METTO DENTRO INTENSITA NEGATIVA O TUTTE POSITIVE? ------------


#cambiamento di shtato. 
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

Qsi0si1prova=Qsi0si1f(1,3,matrice=df,Nv=7)
Qsi0si1prova[1:7,15:21]


Qsi0_3si1_1prova=Qsi0si1f(3,1,matrice=df,Nv=7)
Qsi0_3si1_1prova[15:21,1:7]
##viste e corrette 20/07

#inizializzazioni per debugging
K=4;k=4;n=3;v=4
#inizializzazione di alpha: non sappiamo niente!
#pero dovremmo inizializzarla magari a caso ma solo sugli stati di 1, dove 
#parte la nostra storia 

alpha0<-rep(0,K*(n+v))
alpha0[1:Nv]<-1/Nv
alpha0
#strutturazione canonica
(alpha0<-rep(1/(K*(n+v)),K*(n+v)) )
#1) DOmanda per il prof: come si inizializza alpha? il Po dove lo tr --------


#funzione che forwarda alphat
alphatpiu1<-function(alphat,Qsi,ti1,ti0,Qsi0_si1){
  #mi computa la formula tramite la quale passo da alphat ad alphat+1
  alphat%*%expm(Qsi*(ti1-ti0))%*%Qsi0_si1
}

#chiamata di prova alla funzione che va da 0 a 0.7 da stato 1 a 2
round(expm(Qsip*0.7),2) #fuori dal primo quadrante mette tutti zero e 1 sulla diagionale
(alpha1=alphatpiu1(alpha0,Qsi=Qsip,ti1=0.7,ti0=0,Qsi0_si1 = Qsi0si1prova))

#alpha solo per fare la prova: anche cosi esce un problema con prob>1
alphaprovazero=rep(0,k*Nv); alphaprovazero[15:21]<-1/k*Nv
(alpha1=alphatpiu1(alphaprovazero,Qsi=Qsiprova,ti1=0.7,ti0=0,Qsi0_si1 = Qsi0_3si1_1prova))

alphaprovazero=rep(0,k*Nv); alphaprovazero[8:15]<-1/k*Nv
(alpha1=alphatpiu1(alphaprovazero,Qsi=Qsif(2,matrice=df,Nv=7),ti1=0.7,ti0=0,Qsi0_si1 = Qsi0si1f(2,1,matrice=df,Nv=7)))
#esce problema anche qui!!
# anche con df2 con max 100 mi esce un valore cosi alto 
#esce alpha1 con valorei diversi da 0 solamente in 15:21! giusto, siamo 
#P(xt=i E osservare quella storia)

#gormulazione alpha0 peso uguale ovunque
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]       [,15]
# [1,]    0    0    0    0    0    0    0    0    0     0     0     0     0     0 0.005675817
#             [,16]       [,17]     [,18]       [,19]       [,20]       [,21] [,22] [,23] [,24]
# [1,] 0.005675817 0.005675817 0.1929331 0.005675817 0.005675817 0.005675817     0     0     0
#       [,25] [,26] [,27] [,28]
# [1,]     0     0     0     0

# #formulazione alpha0 peso uguale sono in stato1
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]      [,15]
# [1,]    0    0    0    0    0    0    0    0    0     0     0     0     0     0 0.02270327
#            [,16]      [,17]     [,18]      [,19]      [,20]      [,21] [,22] [,23] [,24] [,25]
# [1,] 0.02270327 0.02270327 0.7717324 0.02270327 0.02270327 0.02270327     0     0     0     0
#       [,26] [,27] [,28]
# [1,]     0     0     0

#alph è naturalmente 1*28
#varie prove di dimensionalità
alphatpiu1(alpha0,Qsi=Qsip,ti1=0.7,ti0=0,Qsi0_si1 = Qsi0si1prova)
df;str(df);str(expm(df*0.7));str(alpha0)
df<-as.matrix(df)
# alpha0[1]=100 #nb qui cambia akpga, rompe le palle poi??
dim(alpha0%*%expm(df))
dim(expm(df)%*%Qsi0si1prova)

##NB MANCA LA FUNZIONE CHE CALCOLA ALPHATAU FINALE
###adesso calcoliamo la funzione che invece fa i beta

#inizializzazione di betatau=1
betatau<-rep(1,K*(n+v))
dim(betatau)
betatif<-function(betat1,Qsi,ti1,ti0,Qsi0_si1){
  ###prenda il beta precedente, la matrice Qsi con zero altrove 
  #la matrice della transizione Qsi0_si1 e i due tempi dove avviene e computa il nuovo beta. 
  betat0<-Qsi0_si1%*%expm(Qsi*(ti1-ti0))%*%(betat1)
  return(betat0)
  #NB Qsi0,si1 sono in realta si-1 , si ; mentre abbiamo ti+1-ti
}

#prove
#come fosse beta2.7--> stato nuovo è1, quindi Qsip
(beta_taumeno1=betatif(betatau,Qsi=Qsip,ti1=3,ti0=2.7,Qsi0_si1 = Qsi0_3si1_1prova))
(beta_taumeno2=betatif(beta_taumeno1,Qsi=Qsiprova,ti1=2.7,ti0=2.5,Qsi0_si1=Qsi0si1prova))
round(t(beta_taumeno1),2)
#è assurdo avere una probabilità amggiore di 1!
# > round(t(beta_taumeno1),2)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17]
# [1,]    0    0    0    0    0    0    0    0    0     0     0     0  0 0  0.15  0.15    0.15
#      [,18]  [,19]   [,20] [,21] [,22] [,23] [,24] [,25] [,26] [,27] [,28]
# [1,] 0.15 104010.7  0.15  0.15     0     0     0     0     0     0     0
#vediamo che escono con proababilita maggiore di zero solo gli appartenentei al terzo sottosistema 
##con PROB altissima se è 19, per via della TRANSIZIONE IMMEDIATA in 19



#2) PROBLEMA: PROB MAGGIORE DI 1!!! -----------------------------------------


dim(beta_taumeno1)
dim(expm(Qsiprova*(0.7))%*%t(betatau))
dim(betatau)
betatau
#beta è 28*1 va bene???? anche se NON so bene perche!

#computa alpha(t+1)-
alphatpiu1meno<-function(alphat,Qsi,ti1,ti0){
  alphat%*%expm(Qsi*(ti1-ti0))
}
 
#prove
alpha1=alphatpiu1meno(alpha0,Qsi=Qsip,ti1=0.7,ti0=0)
round(alpha1,2)
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16]
# [1,] 1.21 0.57 0.17    0    0 0.16 0.16    0    0     0     0     0     0     0     0     0
# [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26] [,27] [,28]
# [1,]     0     0     0     0     0     0     0     0     0     0     0     0

#con alpha diffuso ovunque
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16]
# [1,]  0.3 0.14 0.04    0    0 0.04 0.04 0.04 0.04  0.04  0.04  0.04  0.04  0.04  0.04  0.04
# [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26] [,27] [,28]
# [1,]  0.04  0.04  0.04  0.04  0.04  0.04  0.04  0.04  0.04  0.04  0.04  0.04

#vediamo che anche da qua si vede come senza la transizione, la probabilita è positiva SOLO
#prima della transizione. alphat-=P[xt-=i E sigma o:t]
##quindi vuole dire che è GIUSTO alpha non diffuso
dim(alpha1)




#funzione che computa beta(t)+
betatifpiu<-function(betat1,Qsi,ti1,ti0){ #P[sigmat:tau | Xt+=i]
  betat0<-expm(Qsi*(ti1-ti0))%*%betat1
  return(betat0)
}

#prove #Qsip è ancora una matrice in 1 come stato
(betataumeno1meno=betatifpiu(betatau,Qsi=Qsip,ti1=3,ti0=2.7))
(betataumeno2piu=betatifpiu(beta_taumeno1,Qsi=Qsiprova,ti1=2.7,ti0=2.4))
#vediamo pero che toglie tutti i numeri= ad 1 con il passaggio rpima invece
#è come se al posto che farmi una evidenza diffusa facesse un point evidence. 

# Problema: qunado mi esponenzializza una matrice con degli zeri ----------
##### mi tornano 1!!!

dim(betataumeno1meno)
round(betataumeno1meno, 2) ##assurdo, le probabilita non possono essere cosi alte!!!
dim(as.vector(betataumeno1meno))
dim(betataumeno1meno)
# [1] 28  1
###???!??!?!?!!?

# beta cumulato -----------------------------------------------------------


# PROBLEM:: QSI AL POSTO CHE QSI-1 ----------------------------------------

# possibile errore, sto prendendo lo stato sbagliato probabilmente --------
### pero teoricamente puo funzionare perche tutte le funzioni partono "indietro di un giro"
##sia sti che Sti-1 che ti+1, ti
#partono tutte indietro

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
#mi ritorna betacumulato shiftato di una posizione!


sigmaosservato1

bbb=betacumulatofn(matrice=df,storia=sigmaosservato1,Nv=7,k=4)
bbb

##NB giusto betacumulato fn alternativo qua 

# betacumulatofnalt<-function(matrice,storia,Nv,k ){
#   #versione alternativa che usa Qsi,si+1 al posto che Qsi-1,si
#   betacumulato<-array(1,c(k*Nv,nrow(storia)))
#   for (i in (nrow(storia)-1):1){
#     iesimainiz=storia$stato[i]
#     iesimafin=storia$stato[i+1]
#     print(c(storia$stato[i],storia$stato[i+1]))
#     print("-----")
#     Qsi<-Qsif(iesima=iesimainiz,matrice=matrice,Nv=Nv)
#     Qsi1si0<-Qsi0si1f(iesima0 =iesimainiz ,iesima1 =iesimafin, matrice = matrice,Nv=Nv)
#     ti1<-storia$fine[i+1]
#     # print(storia$fine[i+1])
#     ti0<-storia$inizio[i+1]
#     print(c(storia$inizio[i+1],storia$fine[i+1]))
#     betacumulato[,i]<-betatif(betat1 =betacumulato[,i+1],Qsi = Qsi,ti1 = ti1
#                               ,ti0= ti0,Qsi0_si1 =Qsi1si0  )
#     
#   } 
#   return(betacumulato)
# }

# #prove
# bbb2=betacumulatofnalt(matrice=df,storia=sigmaosservato1,Nv=7,k=4)
# bbb2
# #secondo me ha piu senso bbb2, al quale aggiungerei  in questo caso betacumulatof0
# sigmaosservato1


#funzione che computa beta0
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
#prove
#bbb2[,1]


#versione alternativa non funziona, meglio l'altra
# bbb0=betacumulatof0(matrice=df,storia=sigmaosservato1,betacumulprimo = bbb2[,1],Nv=7,k=4)
# betacumulatof0()
#versione base 
bbb0=betacumulatof0(matrice=df,storia=sigmaosservato1,betacumulprimo = bbb[,2],Nv=7,k=4)


#dovrebbe funzionare!
##secondo me funziona con bbb2[,1]
bbb0 #diverso solo in 1 da 0, certo!!



# mettiamo assieme i due  -------------------------------------------------

# betacumfinalealt<-function(matrice,storia,Nv,k){
#   bbb2=betacumulatofnalt(matrice=matrice,storia=storia,Nv=7,k=4)
#   bbb0=betacumulatof0(matrice=matrice,storia=storia,betacumulprimo = bbb2[,1],Nv=7,k=4)
#   
#   bbbfin<-cbind(bbb0,bbb2)
#   return(bbbfin)
#   
# }
# betacumfinalealt(matrice=df,storia=sigmaosservato1,Nv=7,k=4)


####versione non alternativa cumulata
betacumfinalenonalt<-function(matrice,storia,Nv,k){
  #crea vettore beta giusto
  bbb=betacumulatofn(matrice=matrice,storia=sigmaosservato1,Nv=7,k=4)
  bbb0=betacumulatof0(matrice=matrice,storia=sigmaosservato1,betacumulprimo = bbb[,2],Nv=7,k=4)
  bbb[,1]<-bbb0
  bbbfin<-bbb
  return(bbbfin)
  
}

####dovrebbe funzionare!
#prove: 
betacmfgn=betacumfinalenonalt(matrice=df,storia=sigmaosservato1,Nv=7,k=4)

round(betatot[1:8,])
round(betacmfgn[,9],5)
beta_taumeno1
round(betacmfgn[,8],5)
beta_taumeno2
#la funzione è uguale alle singole... il problema è altrove!!

beta;sigmaosservato1


# alpha cumulata ----------------------------------------------------------
##NNB per noi ti+1 è il fine, altrimenti 
##NB alpha1 è in corrispondenza di quando ho osservato la prima transizione
alphacumulatofn<-function(matrice,storia,Nv,k ){
  #computa alpha da 0 (inizializzata da Po che non conosco), finisco a tau-1
  alpha0<-rep(0.5,k*Nv) #alphao
  alphacumulato<-array(0.5,c(k*Nv,nrow(storia)+1))
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

#prove
aaa=alphacumulatofn(matrice=df,storia=sigmaosservato1,Nv=7,k=4)
sigmaosservato1
colnames(aaa)<-c(0,sigmaosservato1$fine)
aaa


# adesso fare la stesa cosa ma per beta+ e alpha- -------------------------

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


#prove
aaa
alphacumulftau(matrice=df,storia=sigmaosservato1,alphacumultaumenouno = aaa[,ncol(aaa)-1],7,4)


#Adesso manca la funzione per computare alphaTau

alphacumfinalef<-function(matrice,storia,Nv,k){
  aaa=alphacumulatofn(matrice=matrice,storia=storia,Nv=7,k=4)
  aaa0=betacumulatof0(matrice=matrice,storia=storia,betacumulprimo = aaa[,ncol(aaa)-1],Nv=7,k=4)
  #sostituisci ad alpha_tau<-- la sua rappresentazione con alphatau
  aaa[,ncol(aaa)]<-aaa0
  aaafin<-aaa
  return(aaafin)
  
}

# question per il prof: ha senso che alpha e beta di dimensioni di --------

###provo a fare un esempio di transizione osservata

sigmaosservato1=data.frame(stato=c(1,3,2,1,3,4,1,3,1),inizio=c(0,0.5,0.9,1.3,1.5,1.9,2.2,2.5,2.7),
           fine=c(0.5,0.9,1.3,1.5,1.9,2.2,2.5,2.7,3))

#f(sigmaosservato1)
sigmaosservato2=data.frame(stato=c(1,2,1,2,1,2,1,2,1),inizio=c(0,0.4,0.9,1.1,1.5,2,2.2,2.5,2.7),
                          fine=c(0.4,0.9,1.1,1.5,2,2.2,2.5,2.7,3))
##una transizoone osservat



# DOMANDA PER IL PROF: QUANDO PRENDO LE MATRICI QS1,S2 BASTA CHE P --------

##IN QUESTO CASO SE VOGLIO SAPERE TEMPO MEDIO atteso PASSATO IN Qx1h1
#basta fare 
# sum() di in intervalli di evidenza costante quando stiamo nell'uno dell'integrale
##singolo integralino 
#es: stato x1h1
#primo integrale da 0 a 0.5 
Qsip=Qsif(1,matrice=df,Nv=7)#[1:15,1:15]
#crea una matrice di prova in campo1, x1, matrice=df, amalgamata del prof. 
##Nv= somma di N+v

dimensioni=k*(N+v) #le dimensioni totali della matrice
Deltaj=matrix(0,nrow=dimensioni,ncol=dimensioni)
Deltaj[1,1]<-1 


# Questa sara la funzione che ci dara il singolo contributino, il singolo integralino

integralinosingoloET<-function(Qsip,low,up,alpha,beta,Deltaj){
  Etj=0 #mi computa integrale che ho dentro alla parte di E[t]
  for (t in seq(low,up,length.out = 10)){ #meglio un numero piu basso?
    #formula
    Etj=Etj+(alpha%*%expm(Qsip*(t-low))%*%Deltaj%*%expm(Qsip*(up-t))%*%beta)/10 
  }
  Etj=Etj*(up-low)
  return(Etj)
}

#prove
betaf
integralinosingoloET(Qsip=Qsi,low=0,up=0.5,alpha=alphaf,beta=betaf,Deltaj = Deltaj)
##dovrebbe funzionare.:!
integralinosingoloET(Qsip=Qsip,low=0,up=0.5,alpha=alpha0,beta=betazerocinque,Deltaj = Deltaj)
##Qsip, matrice Qs con resto zero. 


# Adesso creo la funzioncina che somma su tutti gli intervalli di  --------


# per qunato riguarda tutti gli h di AMXH di un certo x -------------------

## computa la funzione per un singolo xj selezionato, integra su tutte 
## gli intervalli! Per ogni x (osservato), computiamo lo stesso tempo per ogni 
## Xh singolo
funzioneintegraleXselh<-function(matrice,storia,xsel,Nv,k){
  # xsel=1; 
  Etj=matrix(0,nrow=1,ncol=Nv)
  #alpha vettore
  alphatot=alphacumfinalef(matrice=matrice,storia=storia,Nv=Nv,k=k)
  #.beta vettore
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
        dimensioni=k*(N+v)
        Deltajf=matrix(0,nrow=dimensioni,ncol=dimensioni)
        dim=Nv*(xsel-1)+j
        Deltajf[dim,dim]<-1 
       Etj[,j]<-Etj[,j]+integralinosingoloET(Qsip=Qsi,low=lowf,up=upf,alpha=alphaf,
                                              beta=betaf, Deltaj=Deltajf)
      }
    }
  }
  return(Etj) #torna un VETTORE, di lunghezza Nv, quello dell'integrale!
  
}

#prove
fxsel=funzioneintegraleXselh(matrice=df,storia=sigmaosservato1,xsel=3,Nv=7,k=4)
dim(fxsel)
# [1] 1 7
fxsel
#debugging
Etj=matrix(0,nrow=1,ncol=7)
alphatot=alphacumfinalef(matrice=df,storia=storia,Nv=7,k=4)
alphatot
betatot=betacumulatofn(df,storia,7,4)
betatot
xsel=3
Qsi=Qsif(iesima = xsel,matrice=df,Nv=7)
Qsi[15:21,15:21]
# for (i in 1:nrow(storia)){
i=1
#false
i=2
# i=5
# i=8
# temp=matrix(NA,nrow=7,ncol=4)
# temp[,1]<-fxsel
# temp
# 
# [,1]     [,2]     [,3]      [,4]      [,5]      [,6]      [,7]
# [1,] 13.65099 3.721552 6.877058 0.1955395 0.1543382 0.4267164 0.4267164




# Adesso facciamo funzione che cicla funzione sopra per ogni stato --------

funzioneETtotalesigma<-function(matrice,storia,Nv,k){
  ETime<-matrix(NA,nrow=Nv,ncol=k)
  for (t in 1:k){
    ETime[,t]=funzioneintegraleXselh(matrice=matrice,storia=storia,xsel=t,Nv=Nv,k=k)
  }
  return(ETime) #metto in colonna tutti i risultati di prima, ogni colonna un x diverso
  ##NB potrbbe esserci un problema perche funzioneXselh ritorna una matrice, fa broadcasting?
}


#prove
EE=funzioneETtotalesigma(matrice=df,storia=sigmaosservato1,Nv=7,k=4)
# dd=matrix(1,nrow=28,ncol=28)
# dim(df)
funzioneintegraleXselh(matrice=df,storia=sigmaosservato1,xsel=3,Nv=7,k=4)
round(EE,2)
# Adesso creiamo funzione E[Tj] che somma, a meno di p(sigma0 tau) tutti --------
storia
##gli integralini
sigmaosservato1



##struttura: prende tutta il sigma osservato e dove osserva la x corrispondente va a calcolare.
#prima di far quello faccio una funzione che calcola all'indietro e avanti beta e alpha. 