####Functions for Learning Parameters


# FINITA PARTE DI ET, con una singola sigma -------------------------------
###manca fare la media su tutte le transizioni, e poi la costante p(sigma0:t)



# parte Expected number of transition -------------------------------------


###SI divide in due parti, prima è dove non ci sono transizioni, seconda le transizioni!


# per farlo devo computare la lista di alpha- e beta+ da 1 a 8  -----------

#1--> transizione x1 x3, 2, transizione 3-->2,....#8, tranziosne 3-->1

#     stato inizio fine
# 1     1    0.0  0.5
# 2     3    0.5  0.9
# 3     2    0.9  1.3
# 4     1    1.3  1.5
# 5     3    1.5  1.9
# 6     4    1.9  2.2
# 7     1    2.2  2.5
# 8     3    2.5  2.7
# 9     1    2.7  3.0
#il tempo sara tempo inizio 

# alphatpiu1meno<-function(alphat,Qsi,ti1,ti0){
#   alphat%*%expm(Qsi*(ti1-ti0))
# }

#prove
(alpha1=alphatpiu1meno(alpha0,Qsi=Qsiprova,ti1=0.7,ti0=0))
dim(alpha1) # [1]  1 28


#crea vettore delle alphati-; ci serve per computare M
alpha_globale_meno<-function(matrice,storia,Nv,k ){
  alpha0<-rep(0.5,k*Nv) #da modificare poi 
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
#### sicuri che la scelta migliore sia shiftare? --------------------------
#perche resterebbe l'unica shiftata!!  al contrario di betat e alpha t




#prove
alfatotmeno=alpha_globale_meno(matrice=df,storia=sigmaosservato1,Nv=7,k=4)
#in ogni posizione di alpha globale_meno ci sono gli alpha delle transizioni effettuati in quella 
#posizione nel vettore di stato. 
### NB DIVERSAMENTE DA ALPHA NORMALE E BETA NROMALE DOVE SONO SHIFTATE DI UNO PER PERMETTERE DI ESSERCI 
#QUELLI CON 0. 



# riassunto  --------------------------------------------------------------

#mmm mi sermbra tutto sbagliato questo
alphatot[1:9,] #sono 10 colonne, rispettivamente da alpha0 fino ad alphatau 
#in learning parametri si usano da 1 a 9, rispettivamente da alpha0 a alphatau-1
betatot[1:9,] #sono 10 colonne, rispettivamente da beta0 fino ad betatau 
#in learning si usano da 2 a 10, rispettivamente da beta1 a betatau. 
alfatotmeno[1:8,] #ci sono solamente 8 colonne, da alpha 1 a alpha 8, come 8 sono le transizioni
#che vediamo. manca alpha0 e alphatau. SI USERANNO TUTTE!!


# STESSA COSA CON BETAPIU -------------------------------------------------


##nb per i nostri scopi, usare betacumulatononalt, oppure betacumulatofn( che non computa beta0)
#è la stesa cosa, in quanto non siamo interessati veramente a uare beta0
# betatifpiu<-function(betat1,Qsi,ti1,ti0){
#   
#   
#   betat0<-expm(Qsi*(ti1-ti0))%*%betat1
#   return(betat0)
# }

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

  return(betacumulato[,-c(1,ncol(betacumulato))])
  #sta togliendo la prima e l'ultima colonna, che sarebbe l'11esima, che non ci serve
}


#prove
sigmaosservato1
bbb=beta_globale_piu(matrice=df,storia=sigmaosservato1,Nv=7,k=4)
bbb
##teoricamente adesso è a posto bene!


storia


# COmputiamo adesso intervallo con transizione! ---------------------------
Deltajk<-matrix(0,nrow=28,ncol=28)
j=1; k=1
Deltajk[j,k]<-1
Deltajk[7,8]<-1
Deltajk[15,11]<-1

#alpha_globale_meno()
#MEjk con transizioni.
singoloMEjktrans<-function(matrice,storia,Nv,k,j_me,k_me){
  alfameno<-alpha_globale_meno(matrice=matrice,storia=storia,Nv=Nv,k=k)
  betapiu=beta_globale_piu(matrice,storia,Nv,k)
  MEjk=0
  Deltajk<-matrix(0,nrow=Nv*k,ncol=Nv*k)
  Deltajk[j_me,k_me]<-1
  for (i in 1:(nrow(storia)-1)){
    MEjk<-MEjk+alfameno[,i]%*%Deltajk%*%betapiu[,i]
  }
  
  return(MEjk)
}

alfameno<-alpha_globale_meno(matrice=df,storia=sigmaosservato1,Nv=7,k=4)

betapiu=beta_globale_piu(matrice=df,storia,Nv=7,4)
MEjk=0
for (i in 1:(nrow(storia)-1)){
  MEjk<-MEjk+alfameno[,i]%*%Deltajk%*%betapiu[,i]
  }

MEjk
dim(alfameno)
dim(Deltajk)
dim(betapiu)


##dovrebbe essere giusta questa sommatoria che è il primo membro. per un singolo j,k


##NB MANCA QJK



# Passo successivo: l'integrale di quando non ci sono transizioni ---------




funzioneintegraleMEjk<-function(matrice,storia,Nv,k,j_me,k_me){
  # xsel=1
  Etj=0
  alphatot=alphacumfinalef(matrice=matrice,storia=storia,Nv=Nv,k=k)
  betatot=betacumulatofn(matrice,storia,Nv,k=k)
  
  # Qsi=Qsif(iesima = xsel,matrice=matrice,Nv=Nv)
  for (i in 1:nrow(storia)){
    # if (storia$stato[i]==xsel){
      # print(i)
      # print("Lavora!")
      QS<-Qsif(iesima = storia$stato[i],matrice=matrice,Nv=Nv)
      lowf=storia$inizio[i] #0
      #print("lowf=")
      # print(lowf)
      upf=storia$fine[i] #0.5
      alphaf=alphatot[,i] #alpha0
      betaf=betatot[,i+1] #perche devo prendere betaw
      #beta1
        dimensioni=k*(Nv)
        Deltajf=matrix(0,nrow=dimensioni,ncol=dimensioni)
        # dim=Nv*(xsel-1)+j
        Deltajf[j_me,k_me]<-1 
        Etj<-Etj+integralinosingoloET(Qsip=Qsi,low=lowf,up=upf,alpha=alphaf,
                                              beta=betaf, Deltaj=Deltajf)
  }
  return(Etj)
  
}
#df=df+0.01
df
funzioneintegraleMEjk(matrice=df,storia=sigmaosservato1,Nv=7,k=4,j_me=1,k_me=8)
#
#            [,1]
# [1,] 0.003221869
#transizione tra 1 e 8 sommata su tutto l'oroscopo. 




##adesso dobbiamo fare la somma tra i due e fare il ciclo su tutti i valori

df=df-0.01
df
MEtotale=matrix(0,nrow=k*Nv,ncol=k*Nv)
for( el in 1:(k*Nv)){
  for (el2 in 1:(k*Nv)){
    MEtotale[el,el2]=funzioneintegraleMEjk(matrice=df,storia=sigmaosservato1,Nv=Nv,k=k,
                                           j_me=el,k_me=el2)
    MEtotale[el,el2]=MEtotale[el,el2]+singoloMEjktrans(matrice=df,storia=sigmaosservato1
                                                       ,Nv=Nv,k=k,j_me=el,k_me=el2)

    cat("terminata iterazione numero ", el, el2)
  }
}


1MEtotale

# round(MEtotale,2)
# round(MEtotale,2)
# [,1]  [,2]   [,3]      [,4]   [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12] [,13]  [,14] [,15]
# [1,]  0.27  0.27   0.27 193693.31   0.27  0.27  0.27  0.00  0.00  0.00  0.00  0.00  0.00   0.00 14.09
# [2,]  0.05  0.05   0.05  32839.51   0.05  0.05  0.05  0.00  0.00  0.00  0.00  0.00  0.00   0.00  3.92
# [3,]  0.05  0.05   0.05  32839.51   0.05  0.05  0.05  0.00  0.00  0.00  0.00  0.00  0.00   0.00  1.16
# [4,]  0.05  0.05   0.05  32839.51   0.05  0.05  0.05  0.00  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [5,]  0.61  0.61   0.61  32840.07   0.61  0.61  0.61  0.56  0.56  0.56  0.56  0.56  0.56   0.56  0.55
# [6,]  0.05  0.05   0.05  32839.51   0.05  0.05  0.05  0.00  0.00  0.00  0.00  0.00  0.00   0.00  1.00
# [7,]  0.05  0.05   0.05  32839.51   0.05  0.05  0.05  0.00  0.00  0.00  0.00  0.00  0.00   0.00  1.00
# [8,] 11.53 25.22 107.60  19256.14 114.85  7.37  7.37  0.57  0.16  0.16  0.16  0.16  0.16   0.16  0.38
# [9,] 19.69 43.09 183.85  19256.14 196.23 12.57 12.57  0.57  0.16  0.16  0.16  0.16  0.16   0.16  0.38
# [10,] 18.30 40.05 170.90  19256.14 182.41 11.69 11.69 10.98  3.12  3.12  3.12  3.12  3.12   3.12  0.38
# [11,]  1.94  4.21  17.88  19256.14  19.08  1.25  1.25  0.57  0.16  0.16  0.16  0.16  0.16   0.16  0.38
# [12,]  1.94  4.21  17.88  19256.14  19.08  1.25  1.25  0.57  0.16  0.16  0.16  0.16  0.16   0.16  0.38
# [13,]  0.03  0.04   0.07  19256.14   0.08  0.03  0.03  0.57  0.16  0.16  0.16  0.16  0.16   0.16  0.38
# [14,]  0.03  0.03   0.03  19256.14   0.03  0.03  0.03  0.57  0.16  0.16  0.16  0.16  0.16   0.16  0.38
# [15,]  0.60  0.59   0.51 102607.45   4.70  0.60  0.60 19.28 22.91 22.75 16.62 16.62  0.44 223.58  9.33
# [16,]  0.59  0.58   0.50  88845.38   4.78  0.59  0.59 22.04 26.21 26.03 19.00 19.00  0.45 256.26  6.59
# [17,]  0.24  0.23   0.20  38526.34   1.89  0.24  0.24  9.48 11.27 11.20  8.17  8.17  0.18 110.37  2.77
# [18,]  0.00  0.00   0.00    196.49   0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.01  0.41
# [19,]  0.00  0.00   0.00    234.78   0.00  0.00  0.00  0.01  0.01  0.01  0.01  0.01  0.00   0.12  0.39
# [20,]  0.03  0.03   0.03  19950.76   0.09  0.03  0.03  0.40  0.47  0.47  0.34  0.34  0.01   4.62  0.48
# [21,]  0.03  0.03   0.03  19950.76   0.09  0.03  0.03  0.40  0.47  0.47  0.34  0.34  0.01   4.62  0.48
# [22,]  6.96 14.54  41.89  19256.14  69.22  3.82  3.82  0.00  0.00  0.00  0.00  0.00  0.00   0.00  0.38
# [23,] 32.36 67.73 195.33  19256.15 322.83 17.70 17.70  0.00  0.00  0.00  0.00  0.00  0.00   0.00  0.38
# [24,] 19.19 40.15 115.76  19256.14 191.32 10.50 10.50  0.00  0.00  0.00  0.00  0.00  0.00   0.00  0.38
# [25,]  2.10  4.37  12.56  19256.14  20.75  1.16  1.16  0.00  0.00  0.00  0.00  0.00  0.00   0.00  0.38
# [26,]  2.10  4.37  12.56  19256.14  20.75  1.16  1.16  0.00  0.00  0.00  0.00  0.00  0.00   0.00  0.38
# [27,]  0.03  0.03   0.03  19256.14   0.03  0.03  0.03  0.00  0.00  0.00  0.00  0.00  0.00   0.00  0.38
# [28,]  0.03  0.04   0.05  19256.14   0.07  0.03  0.03  0.00  0.00  0.00  0.00  0.00  0.00   0.00  0.38
# [,16] [,17]  [,18]  [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]  [,27] [,28]
# [1,] 14.09 14.24 140.89   0.00  5.43  5.43  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [2,]  3.92  3.95  39.21   0.00  1.52  1.52  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [3,]  1.16  1.15  11.60   0.00  0.46  0.46  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [4,]  0.00  0.00   0.04   0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [5,]  0.55  0.44   5.45   0.01  0.58  0.58  0.56  0.56  0.56  0.56  0.56   0.56  0.56
# [6,]  1.00  0.97  10.00   0.00  0.41  0.41  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [7,]  1.00  0.97  10.00   0.00  0.41  0.41  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [8,]  0.38  0.30   3.85   0.00  0.21  0.21  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [9,]  0.38  0.30   3.85   0.00  0.21  0.21  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [10,]  0.38  0.30   3.85   0.00  0.21  0.21  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [11,]  0.38  0.30   3.85   0.00  0.21  0.21  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [12,]  0.38  0.30   3.85   0.00  0.21  0.21  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [13,]  0.38  0.30   3.85   0.00  0.21  0.21  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [14,]  0.38  0.30   3.85   0.00  0.21  0.21  0.00  0.00  0.00  0.00  0.00   0.00  0.00
# [15,]  9.33 16.38  89.02 105.51  4.55  4.55 24.24 22.64 22.64 15.55 15.55 222.44  0.44
# [16,]  6.59 13.37  61.45 107.75  3.53  3.53 27.73 25.90 25.90 17.77 17.77 254.93  0.45
# [17,]  2.77  5.44  25.96  42.41  1.51  1.51 11.93 11.14 11.14  7.64  7.64 109.81  0.18
# [18,]  0.41  0.32   4.05   0.00  0.22  0.22  0.00  0.00  0.00  0.00  0.00   0.01  0.00
# [19,]  0.39  0.30   3.88   0.05  0.21  0.21  0.01  0.01  0.01  0.01  0.01   0.12  0.00
# [20,]  0.48  0.50   4.76   1.51  0.26  0.26  0.57  0.53  0.53  0.36  0.36   5.24  0.01
# [21,]  0.48  0.50   4.76   1.51  0.26  0.26  0.57  0.53  0.53  0.36  0.36   5.24  0.01
# [22,]  0.38  0.30   3.85   0.00  0.21  0.21  0.68  0.11  0.11  0.11  0.11   0.11  0.11
# [23,]  0.38  0.30   3.85   0.00  0.21  0.21  0.68  0.11  0.11  0.11  0.11   0.11  0.11
# [24,]  0.38  0.30   3.85   0.00  0.21  0.21 13.18  2.22  2.22  2.22  2.22   2.22  2.22
# [25,]  0.38  0.30   3.85   0.00  0.21  0.21  0.68  0.11  0.11  0.11  0.11   0.11  0.11
# [26,]  0.38  0.30   3.85   0.00  0.21  0.21  0.68  0.11  0.11  0.11  0.11   0.11  0.11
# [27,]  0.38  0.30   3.85   0.00  0.21  0.21  0.68  0.11  0.11  0.11  0.11   0.11  0.11
# [28,]  0.38  0.30   3.85   0.00  0.21  0.21  0.68  0.11  0.11  0.11  0.11   0.11  0.11



MEtotale2=matrix(0,nrow=k*Nv,ncol=k*Nv)
for( el in 1:(k*Nv)){
  for (el2 in 1:(k*Nv)){
    MEtotale2[el,el2]=funzioneintegraleMEjk(matrice=df,storia=sigmaosservato2,Nv=Nv,k=k,
                                           j_me=el,k_me=el2)
    MEtotale2[el,el2]=MEtotale[el,el2]+singoloMEjktrans(matrice=df,storia=sigmaosservato2
                                                       ,Nv=Nv,k=k,j_me=el,k_me=el2)
    
  }
}

round(MEtotale,2)
round(MEtotale2,2)
#ce un possibile problema, che non da peso alle transizioni intra sottostato. 
# write.csv(MEtotale,"primosigma.csv")
# write.csv(MEtotale2,"secondosigma.csv")

df



# Funzione che computa la somma di M --------------------------------------


funzMEtotale<-function(matrice,storia,Nv,k){ 
  MEtotale2=matrix(0,nrow=k*Nv,ncol=k*Nv)
  for( el in 1:(k*Nv)){
    for (el2 in 1:(k*Nv)){
      MEtotale2[el,el2]=funzioneintegraleMEjk(matrice=matrice,storia=storia,Nv=Nv,k=k,
                                            j_me=el,k_me=el2)
      MEtotale2[el,el2]=MEtotale[el,el2]+singoloMEjktrans(matrice=matrice,storia=storia
                                                        ,Nv=Nv,k=k,j_me=el,k_me=el2)
    
      print(paste0("iterazione numero: ",el," ", el2 ))
    }
  }
  return(MEtotale2)
}


funzMEtotale(matrice=df,storia=storia,Nv=7,k=4)


# domanda per il prof, ma questa va fatta anche per la diagonale?? --------

##manca solo qjk nella funzione 





# Funzione adesso che computa E[t] totale per tutti i 28 membri  ----------






funzioneintegraleTEjj<-function(matrice,storia,Nv,k,j_me){
  # xsel=1
  Etj=0
  alphatot=alphacumfinalef(matrice=matrice,storia=storia,Nv=Nv,k=k)
  betatot=betacumulatofn(matrice,storia,Nv,k=k)
  
  # Qsi=Qsif(iesima = xsel,matrice=matrice,Nv=Nv)
  for (i in 1:nrow(storia)){
    # if (storia$stato[i]==xsel){
    # print(i)
    # print("Lavora!")
    QS<-Qsif(iesima = storia$stato[i],matrice=matrice,Nv=Nv)
    lowf=storia$inizio[i] #0
    #print("lowf=")
    # print(lowf)
    upf=storia$fine[i] #0.5
    alphaf=alphatot[,i] #alpha0
    betaf=betatot[,i+1] #perche devo prendere betaw
    #beta1
    dimensioni=k*(Nv)
    Deltajf=matrix(0,nrow=dimensioni,ncol=dimensioni)
    # dim=Nv*(xsel-1)+j
    Deltajf[j_me,j_me]<-1 
    Etj<-Etj+integralinosingoloET(Qsip=Qsi,low=lowf,up=upf,alpha=alphaf,
                                  beta=betaf, Deltaj=Deltajf)
  }
  return(Etj)
  
}



funzioneintegraleTEjj(matrice=df,storia=sigmaosservato1,Nv=7,k=4,j_me=4)


#adesso questa funzione la scrivo complessivamente per ogni j 


ETtotale<-function(matrice,storia,Nv,k){
vettoreET<-vector()
for (i in 1:(k*Nv)){
  vettoreET[i]<-funzioneintegraleTEjj(matrice=matrice,storia=storia,Nv=Nv,k=k,j_me=i)

}
return(vettoreET)
}

ET<-ETtotale(matrice=df,storia=sigmaosservato1,Nv=7,k=4)
# 
# [1] 0.000000e+00 0.000000e+00 0.000000e+00 9.902333e+01 1.841265e-03 0.000000e+00 0.000000e+00
# [8] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
# [15] 2.093136e-02 1.495073e-02 2.072629e-02 6.442006e-03 1.324248e-03 1.945292e-06 1.945292e-06
# [22] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00


##ce sicuramente qualcosa che non va

