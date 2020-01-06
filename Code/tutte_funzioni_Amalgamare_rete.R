
# funzioni amalgamazione e indietro ---------------------------------------
library(stringr)

# creazionefullRapprdatoDirect --------------------------------------------
directnames=function(k_int,N_int,matriced){
  vector=c()
  for (i in 1:k_int){
    for (j in 1:N_int){
      vector=c(vector,paste0(j,"_",i))
      
    }
  }
  # print(vector)
  rownames(matriced)=vector;colnames(matriced)=vector  
  return(matriced)
} 
##importantissima anche per dopo!!


# RAPPRESENTAZIONE FULL: DA DIRECT A FULL

#funzione che crea singola Qx|h1
##funzione che crea Qx|h1 in rappresentazione FULL
singolax=function(k_int,j_int){
  #prende come input due numeri, k e j, e crea la matrice corrisponendte
  mm=matrix(ncol=k_int,nrow=k_int)
  for (i in 1:k_int){
    for (m in 1:k_int){
      mm[i,m]=ifelse(i==m && i!=j_int,-Inf,ifelse(j_int==m & i!=j_int, Inf, 0))
    }
    
  }
  return(mm)
}
vuota=function(k_int){
  m_int=matrix(0,ncol=k_int,nrow=k_int)
  return(m_int)
  
} #crea matrice vuota X!Hi
Funz_tuttex=function(k,N){
  m_full_x=array(NA,c(k,k,N*k))
  
  for (l in 1:(N*k)){ #N*k)
    j=(l+N-1)/N; #print(j)
    ii=ifelse(j!= floor(j), T,F)# vuota(k_int=k),singolax(k_int=k,j_int=j))
    if (ii){
      m_full_x[,,l]=vuota(k_int=k)
    }
    else 
      m_full_x[,,l]=singolax(k_int=k,j_int=j)
  }         
  
  return(m_full_x)
  
  
} #creo tutte le matrici Qx|H

singolah=function(k_int,N_int,j_int,matriced){
  mm=matrix(ncol=k_int*N_int,nrow=k_int*N_int)
  for (l in 1:(k_int*N_int)){
    for (m in 1:(k_int*N_int)){
      # print(l)
      # print(m)
      # print("-----")
      i=l-N_int*(j_int-1)
      u=ceiling(m/N_int)
      r=m-N_int*(u-1)
      mm[l,m]=ifelse(l>N_int*(j_int-1) && l<=N_int*j_int, matriced[paste0(i,"_",j_int),paste0(r,"_",u)],0)
      #mm[i,m]=ifelse(i==m && i!=j_int,-Inf,ifelse(j_int==m & i!=j_int, Inf, 0))
    }
    
  }
  return(mm)
} #singola H|xi
Funz_tutteh_rappresentazioneFULL=function(k,N,matrice){
  m_full_x=array(NA,c(N*k,N*k,k))
  
  for (j in 1:(k)){ #N*k)
    m_full_x[,,j]=singolah(k_int=k,N_int=N,j_int=j,matriced=matrice)
  }
  return(m_full_x)
  
} #tutte 
#gli cambio nome!! senno entra in conflitto


# va,vb,v, ordinamento ----------------------------------------------------

#transizioni non allowate per va(pari) e vb (dispari)
funz_va<-function(k_int){
  ###OVVERO LA FORMULA è: k/2*(k/2 -1) con approx per eccesso
  return(ceiling(k_int/2)*(ceiling(k_int/2)-1))
}
funz_vb<-function(k_int){
  ###OVVERO LA FORMULA è: k/2*(k/2 -1) con approx per DIFETTO. 
  return(floor(k_int/2)*(floor(k_int/2)-1))
}
funz_v<-function(k_int){
  ### combina semplicemente le due sopra, restituisce il numero totale v
  ### di transizioni non allowate
  return(funz_vb(k_int)+funz_va(k_int))
}



#funzioni di ordinamento, per gli eta (servono alle N)
ordinamentoX_H<-function(N,K){ 
  #############a
  #ordinamento è una funzione che si occupa di creare un ordine 
  #congiunto eta_X_H. come parametri nel nostro caso specifico gli diamo
  #numero di stati di X,N, e numero di stati di ogNiX_H H,K.
  ##ritorna l'ordinamento. 
  #############a 
  mm=array(NA,dim=c(N*K,K,2))#ncol=K,nrow=(N*K))
  succ_K=paste0("K",seq(1,K,by=1))
  succ_N=paste0("N",seq(1,N,by=1))
  succe_KN=array(NA,dim=c(N,K,2))
  for(j in 1:length(succ_K)){
    
    
    for(i in 1:length(succ_N)){
      succe_KN[i,j,1]=succ_K[j]
      succe_KN[i,j,2]=succ_N[i]
      
    }
  }
  m_f=matrix(nrow=length(c(succe_KN[,,1])),ncol=2)
  for (i in 1:length(c(succe_KN[,,1]))){
    m_f[i,1]=c(succe_KN[,,1])[i] #ordinamento opposto: c(t(succe....))
    m_f[i,2]=c(succe_KN[,,2])[i] #ordinamento opposto: c(t(succe....))
  }
  return(m_f)
} #ordine congiunto eta_X_H
#ordinamenti parziali invece di X e H 
ordinamentoX<-function(K) return(paste0("K",seq(1,K,by=1)))
ordinamentoH<-function(N) return(paste0("N",seq(1,N,by=1)))

#funzioni che computano le matrici NX|H e NH|X che servono per amalgamare
singolaNX_hi<-function(K,N,i,ord,ordH,ordX){
  ### calcola la matrice NX|hi, avendo in input tutti i vari ordinamenti, a quali 
  ### indice i ci riferiamo, numero di stati di X e #di stati di N 
  NiX_H=matrix(NA,nrow=K*N,ncol=K)
  
  for (j in 1:nrow(NiX_H)){
    for (k in 1:ncol(NiX_H)){
      NiX_H[j,k]=ifelse(ord[j,2]==ordH[i] & ord[j,1]==ordX[k],1,0 )
    }
  }
  
  return(NiX_H)
} #singola NX|H, per h specifico
singolaNH_xi<-function(K,N,i,ord,ordH,ordX){
  ### calcola la matrice NH_xi, avendo in input tutti i vari ordinamenti, a quali 
  ### indice i ci riferiamo, numero di stati di X e #di stati di N 
  NiH_X=matrix(NA,nrow=K*N,ncol=N)
  for (j in 1:nrow(NiH_X)){
    for (k in 1:ncol(NiH_X)){
      NiH_X[j,k]=ifelse(ord[j,1]==ordX[i] & ord[j,2]==ordH[k],1,0 )
    }
    
  }
  return(NiH_X)
} #singola NH|X
completaNX_H<-function(K,N){
  ###mette assieme tutte le matricini NX|hi di prima per formare
  #un array tridimensionale di dimensioni K*N, K, N 
  ##dove in ogni elemento bidimensionale K*N, K CE una singola matrice NX|hi
  # K=2;N=3
  NX_H=array(NA,dim=c(K*N,K,N))
  ordn=ordinamentoX_H(N,K)
  ordnH=ordinamentoH(N)
  ordnX=ordinamentoX(K)
  
  for (i in 1:N){
    temp=singolaNX_hi(K,N,i,ordn,ordnH,ordnX)
    NX_H[,,i]=temp
    
  }
  return(NX_H)
}
completaNH_X<-function(K,N){
  
  ###mette assieme tutte le matricini NH_xi di prima per formare
  #un array tridimensionale di dimensioni K*N, N, K
  ##dove in ogni elemento bidimensionale K*N, N CE una singola matrice NH_xi
  # K=2;N=3
  # K=2;N=3
  NH_X=array(NA,dim=c(K*N,N,K))
  ordn=ordinamentoX_H(N,K)
  ordnH=ordinamentoH(N)
  ordnX=ordinamentoX(K)
  
  for (i in 1:K){
    temp=singolaNH_xi(K,N,i,ordn,ordnH,ordnX)
    NH_X[,,i]=temp
    
  }
  return(NH_X)
} #tutte le matrici NH|X



# bipartitica -------------------------------------------------------------

#funzioni di supporto e ordinameto per QX_l
funzsupporto<-function(k) {
  #prende in input k #stati di X e restituisce le transazioni non allowate 
  A<-seq(1,k,2)
  B<-seq(2,k,2)
  nl=vector(); nn=vector();  cont=1
  for (i in 1:length(A)){
    for (j in 1:length(A)){
      if (i!=j) {
        nl[cont]<-paste0(A[i],"_",A[j])
        nn[cont]<-paste0(B[i],"_",B[j])
        cont=cont+1
      }
    }
  }
  nn<-nn[!str_detect(nn,"NA")]
  return(c(nn,nl)  )
}
funzioneordinamento<-function(k){
  ff=funzsupporto(k)
  z=data.frame(prim=as.numeric(substr(ff,1,1)),sec=as.numeric(substr(ff,3,3)))
  z$terz=ifelse(z$prim>z$sec,1,0)
  z1=z[z$terz==0,];(z1<-z1[order(z1$prim,z1$sec),]);z1$ris=2*(1:nrow(z1)-1)+1
  z2=z[z$terz!=0,]; z2<-z2[order(z2$sec,z2$prim),];z2$ris=2*(1:nrow(z2))#;z2
  zfin=rbind(z1,z2);zfin=zfin[order(zfin$ris),c(1,2,4)]
  return(zfin)
  
}


#funzioni per computare QX|H nei vari casi
QX_h1<-function(K,N,mm){
  ##prende in input #stati di X, #stati di N, e la matrice mm in forma diretta
  ## restituisce una matrice QX|h1: la prima matrice bipartitica
  matrice=matrix(NA,ncol=K,nrow=K)
  for (i in 1:K){
    for (j in 1:K){
      matrice[i,j]=ifelse(i==j & i%%2==0,NA,ifelse(i!=j & i%%2==0 & j%%2!=0, mm[paste0(N,"_",i),paste0(1,"_",j)],0 ))
      # print(matrice)                    
      
    }
    temp=matrice[i,]
    temp[is.na(temp)]<--sum(temp,na.rm=T)
    matrice[i,]<-temp
  }
  return(matrice)
} ## Matrice bipartitica dato h1
QX_hN<-function(K,N,mm){
  ###altro caso particolare è la matrice QX|Hn, anche in questo caso input sono K,N e 
  ###la matrice in forma diretta. Restituisce una matrice K*K. 
  matrice=matrix(NA,ncol=K,nrow=K)
  for (i in 1:K){
    for (j in 1:K){
      matrice[i,j]=ifelse(i==j & i%%2!=0,NA,ifelse(i!=j & i%%2!=0 & j%%2==0, mm[paste0(N,"_",i),paste0(1,"_",j)],0 ))
      # print(matrice)                    
      
    }
    temp=matrice[i,]
    temp[is.na(temp)]<--sum(temp,na.rm=T)
    matrice[i,]<-temp
  }
  
  return(matrice)
} # matrice bipartitica dato hN
QX_h0<-function(K,N=1,mm=matrix(c(-1,1,0,0,1,-2,1,0,0,0,-0.5,0.5,0.6,0,0.4,-1)
,ncol=4,byrow=T)){##matrice vuota di zeri quando l sta tra 1 e N 
  #dovrei cambiare input: non mi serve a niente avere una matrice. Anche N è tralasciabile 
  return(matrice=matrix(0,ncol=K,nrow=K))
  
} # matrice quando 1<l<N ---> matrice vuota di zeri
QX_l<-function(K,N,l){
  ###FUNZONE che computa la matrice QX|l quando l>N ; l<=N+v
  ###ovvero quelli stati di h ausiliari che servono per far avvenire le 
  ### transazioni non allowate. 
  ### Input: K,N, e l'indice l. Restituisce una matrice K*K
  ff=funzioneordinamento(k=K)
  a=ff[ff$ris==(l-N),1]
  b=ff[ff$ris==(l-N),2]
  # l=4
  matrice=matrix(NA,ncol=K,nrow=K)
  for (i in 1:K){
    for (j in 1:K){
      matrice[i,j]=ifelse(i==j & i==a,-Inf,ifelse(i==a & j==b, Inf,0))
      # print(matrice)                    
    }
  }
  return(matrice)
} ###FUNZONE che computa la matrice QX|l quando l>N ; l<=N+v
QX_Htotale<-function(k,N,matrice){
  ### funzione che mette assieme i pezzi costruiti prima per costruire tutte le 
  ### matrici di QX|H. 
  ### Input: k,N e la matrice in forma diretta. 
  ### Restituisce un array K, K, N+v. 
  ### dove in ogni singolo elemento K,K ce una matrice QX|hi
  v=funz_v(k)
  m_full_x=array(NA,c(k,k,N+v))
  
  for (j in 1:(N+v)){ #N*k)
    tt=data.frame()
    if (j==1){
      m_full_x[,,j]=QX_h1(K=k,N=N,mm=matrice)
      
    } else if(j<N){
      m_full_x[,,j]=QX_h0(K=k)
      
      
    } else if(j==N){
      m_full_x[,,j]=QX_hN(K=k,N=N,mm=matrice)
      
      
      
    } else {
      m_full_x[,,j]=QX_l(K=k,N=N,l=j)
      
      
    }
    # print(tt)
    # m_full_x[,,j]=tt
  }
  return(m_full_x)
  
}
#crea tutte le matrici QX|H

#funzione di supporto per QH|X
funz_quartopuntoqhlhm_x<-function(K,N,j){
  ##funzione che mi serve per quarto punto, sia apri che dispari
  #della matrice QH|X
  v=funz_v(K) #ordinamento
  c<-0
  s<-vector()
  for (l in (N+1):(N+v)){
    a=QX_l(K=K,N=N,l=l) ##rifaccio al contrario passo QX|Hl con l>N
    #per ogni l>N 
    # print(a)
    # print(l)
    if (Inf %in% a[,j]){
      #se in tale matrice,alla j-esima colonna, è presente un Inf, 
      ##allora scrivo tale l nel mio vettore s 
      c<-c+1
      s[c]<-l
    }
    
  }
  ifelse(length(s)>0,return(s),return(FALSE))
  ##controllo per restituire sempre qualcosa di diverso dall'insieme vuoto. 
} ##si trova in esempioprimarete3_3
#computo tutte le matrici QH|X
QH_Xdispari<-function(N,K,v,va,j,mm){
  ##Funzione che computa QH|xj, quando l'indice della x, j, è dispari. 
  ###prende in input, N,K #di stati rispettivamente di H e X, 
  ### v, va transizioni non allowate totali, e transizioni non allowate 
  ### all'interno del gruppo A (quello dispari.)
  ### j, indice della X, e mm matrice di rappresentazione diretta. 
  
  ### Restituisce una matrice N+v * N+v
  s<-funz_quartopuntoqhlhm_x(K=K,N=N,j=j)
  matrice=matrix(NA,ncol=N+v,nrow=N+v)
  ff=funzioneordinamento(k=K)
  for (l in 1:(N+v)){
    for (m in 1:(N+v)){
      if (l==m & m<N){
        matrice[l,m]=mm[paste0(l,"_",j),paste0(l,"_",j)]
        
      } else if(l==m & m==N){
        tempo=0
        for (u in 1:K){
          tempo<-tempo+ifelse(u%%2==0,mm[paste0(N,"_",j),paste0(1,"_",u)],0)
          ##ciclo che mi porta dentro la somma quando u è pari
          # print(tempo)
          
        }
        matrice[l,m]=mm[paste0(N,"_",j),paste0(N,"_",j)]+tempo
        # print("bbbb")
        # print(matrice[l,m])
        
      } else if(l!=m & l<=N & m<=N){
        matrice[l,m]=mm[paste0(l,"_",j),paste0(m,"_",j)]
        
        
        
        
      } else if(l==N & m>N & m<=N+va ){
        if( ff[ff$ris==abs((m-N)),1]==j){ #ho messo un if dentro poichè altrimenti 
          #dava errore che ci sono casi dove la chiamata dell'if è indefinita (out of bounds)
          a=ff[ff$ris==(m-N),1]
          b=ff[ff$ris==abs((m-N)),2]
          # print(a)
          # print("...")
          # print(b)
          # print("...")
          
          matrice[l,m]=mm[paste0(N,"_",j),paste0(1,"_",b)]
        } else matrice[l,m]=0
        # print("ddddd")
        # print(matrice[l,m])
        
        
      }
      else if(l>N & l<=N+va & l %in% s & m==l){
        ###per inserire i +-Infinito NON FUNZIONA la formula del Prof, ho usato questa formula 
        ###con s, con s funzione AD-Hoc per questo punto. 
        matrice[l,m]=-Inf
        #problema: SBAGLIATO!!!
        # l>N & l<=N+va & l==ff[ff$ris==K+1-j,2]+N & m==l
        
        
        
      }#j 3--->1  j 1 --->2
      ##possibile sol K-j per dispari e K-j+1 per pari 
      else if(l>N & l<=N+va & l %in% s & m==1){
        matrice[l,m]=Inf
        #problema: SBAGLIATO!!!
        
        
        
        
      }
      else {
        matrice[l,m]=0
        
        
      }
      # print(tt)
      # m_full_x[,,j]=tt
    }
    
    
  }
  
  return(matrice)  
} #versione sovrascritta!!!!
QH_Xpari<-function(N,K,v,va,j,mm){
  
  ##Funzione che computa QH|xj, quando l'indice della x, j, è PARI. 
  ###prende in input, N,K #di stati rispettivamente di H e X, 
  ### v, va transizioni non allowate totali, e transizioni non allowate 
  ### all'interno del gruppo A (quello dispari.) (IN QUANTO andra ad agire da N+va fino a N+v)
  ### j, indice della X, e mm matrice di rappresentazione diretta. 
  
  ### Restituisce una matrice N+v * N+v
  matrice=matrix(NA,ncol=N+v,nrow=N+v)
  s<-funz_quartopuntoqhlhm_x(K=K,N=N,j=j)
  ff=funzioneordinamento(k=K)
  for (l in 1:(N+v)){
    for (m in 1:(N+v)){
      if (l==m & m<=N & l>1){
        matrice[l,m]=mm[paste0(N-l+1,"_",j),paste0(N-l+1,"_",j)]
        
      } else if(l==m & m==1){
        tempo=0
        for (u in 1:K){
          tempo<-tempo+ifelse(u%%2==1,mm[paste0(N,"_",j),paste0(1,"_",u)],0)
          ##DIFFERENZA CON IL PROF: 1_u alposto che N_u! 
          # print(tempo)
          
        }
        # print(tempo)
        # print(mm[paste0(N,"_",j),paste0(N,"_",j)])
        # print(matrice[l,m])
        matrice[l,m]=mm[paste0(N,"_",j),paste0(N,"_",j)]+tempo
        # print("bbbb")
        # print(matrice[l,m])
        
      } else if(l!=m & l<=N & m<=N){ #3
        matrice[l,m]=mm[paste0(N-l+1,"_",j),paste0(N-m+1,"_",j)]
        
        
        
        
      } else if(l==1 & m>N+va & m<=N+v ){#4
        if( ff[ff$ris==abs((m-N)),1]==j){##come sopra: if interno altrimenti va out of bound
          ##si basa dul dire che la a deve essere a==j
          a=ff[ff$ris==(m-N),1]
          b=ff[ff$ris==(m-N),2]
          matrice[l,m]=mm[paste0(N,"_",j),paste0(1,"_",b)]
          ## DIVERSO DAL PROF!!!!
        } else matrice[l,m]=0
        # print("ddddd")
        # print(matrice[l,m])
        
        # DIVERSO DAL PROF #4 PARI!!
        
        
        
      }
      else if(l>N+va & l<=N+v & l %in% s & m==l){
        ####come in dispari anche questo 5,6 diversi dal PROF
        matrice[l,m]=-Inf
        #sbaglaito
        # l>N+va & l<=N+v & l==ff[ff$ris==K+1-j,2]+N & m==l
        
        
        
      }
      else if(l>N+va & l<=N+v & l %in% s & m==N){
        matrice[l,m]=Inf
        
        
        
        
      }
      else {
        matrice[l,m]=0
        
        
      }
      # print(tt)
      # m_full_x[,,j]=tt
    }
    
    
  }
  
  return(matrice)  
} ##versione sovrascritta!!!!
Funz_tutteh=function(k,N,matrice){
  ###combina il pari e dispari di prima per creare un array, 
  # di dimensioni N+v, N+v, k, dove ogni elemento n+v*n+v è una matrice
  ## QH|xj. 
  ## input: k,N soliti e la matrice in forma diretta
  va<-funz_va(k_int=k)
  v=funz_v(k_int=k)
  m_full_h=array(NA,c(N+v,N+v,k))
  for (j in 1:k){
    if (j%%2==0 ){
      m_full_h[,,j]=QH_Xpari(K=k,N=N,va=va,v=v,j=j,mm=matrice)
      
    } else if(j%%2!=0  ){
      m_full_h[,,j]=QH_Xdispari(K=k,N=N,va=va,v=v,j=j,mm=matrice)
      
      
    } 
    # print(tt)
    # m_full_x[,,j]=tt
  }
  
  return(m_full_h)
}
#funzione che mi crea tutte le matrici QH|X


#matrice amalgamata H|X
amalgamataH_X<-function(k,N,matrice){
  ###funzione che amalgama tutte le matrice QH|X in una matriciona di 
  ### dimensioni k*(N+v), k*(N+v) amalgamata. 
  matrice=directnames(matriced=matrice,k_int=k,N_int = N)
  aus=Funz_tutteh(k=k,N=N,matrice=matrice)
  aus[aus==Inf]<-10000; aus[aus==-Inf]<--1000
  v=funz_v(k_int=k)
  tot=N+v
  NN=completaNH_X(K=k,N=tot) 
  Matr_amalgH_X<-matrix(0,nrow=k*(N+v),ncol=k*(N+v))
  for(i in 1:k){
    Matr_amalgH_X<-Matr_amalgH_X+NN[,,i]%*%aus[,,i]%*%t(NN[,,i])
    
  }
  return(Matr_amalgH_X)
  
  
}
#matrice amalgamata X|H
amalgamataX_H<-function(k,N,matrice){
  ###funzione che amalgama tutte le matrice QX|H in una matriciona di 
  ### dimensioni k*(N+v), k*(N+v) amalgamata. 
  
  # matrice=directnames(matriced=matrice,k_int=k,N_int = N)
  aus2=QX_Htotale(k=k,N=N,matrice=matrice)
  aus2[aus2==Inf]<-10000
  aus2[aus2==-Inf]<--1000
  v=funz_v(k_int=k)
  tot=N+v
  NN=completaNX_H(K=k,N=tot)
  Matr_amalgX_H<-matrix(0,nrow=k*(N+v),ncol=k*(N+v))
  
  for(i in 1:(N+v)){
    Matr_amalgX_H<-Matr_amalgX_H+NN[,,i]%*%aus2[,,i]%*%t(NN[,,i])
    
  }
  return(Matr_amalgX_H)
}

# NuovaversionediQH_Xpariedispari QH|X -----------------------------------------
#strutturazione alternativa di QH|X e QX|H, corretta!! 
#e consistente per ogni dimensione
#creano matrice QH|X_pari e se x è dispari 
QH_Xpari<-function(N,K,v,va,j,mm){
  va=1
  
  ##Funzione che computa QH|xj, quando l'indice della x, j, è PARI. 
  ###prende in input, N,K #di stati rispettivamente di H e X, 
  ### v, va transizioni non allowate totali, e transizioni non allowate 
  ### all'interno del gruppo A (quello dispari.) (IN QUANTO andra ad agire da N+va fino a N+v)
  ### j, indice della X, e mm matrice di rappresentazione diretta. 
  
  ### Restituisce una matrice N+v * N+v
  matrice=matrix(NA,ncol=N+v,nrow=N+v)
  s<-funz_quartopuntoqhlhm_x(K=K,N=N,j=j)
  ff=funzioneordinamento(k=K)
  for (l in 1:(N+v)){
    for (m in 1:(N+v)){
      if (l==m & m<=N & l>1){
        matrice[l,m]=mm[paste0(N-l+1,"_",j),paste0(N-l+1,"_",j)]
        
      } else if(l==m & m==1){
        tempo=0
        for (u in 1:K){
          tempo<-tempo+ifelse(u%%2==1,mm[paste0(N,"_",j),paste0(1,"_",u)],0)
          ##DIFFERENZA CON IL PROF: 1_u alposto che N_u! 
          # print(tempo)
          
        }
        # print(tempo)
        # print(mm[paste0(N,"_",j),paste0(N,"_",j)])
        # print(matrice[l,m])
        matrice[l,m]=mm[paste0(N,"_",j),paste0(N,"_",j)]+tempo
        # print("bbbb")
        # print(matrice[l,m])
        
      } else if(l!=m & l<=N & m<=N){ #3
        matrice[l,m]=mm[paste0(N-l+1,"_",j),paste0(N-m+1,"_",j)]
        
        
        
        
      } else if(l==1 & m>N+va & m<=N+v ){#4
        if( ff[ff$ris==abs((m-N)),1]==j){##come sopra: if interno altrimenti va out of bound
          ##si basa dul dire che la a deve essere a==j
          a=ff[ff$ris==(m-N),1]
          b=ff[ff$ris==(m-N),2]
          matrice[l,m]=mm[paste0(N,"_",j),paste0(1,"_",b)]
          ## DIVERSO DAL PROF!!!!
        } else matrice[l,m]=0
        # print("ddddd")
        # print(matrice[l,m])
        
        # DIVERSO DAL PROF #4 PARI
        
        
      }
      else if(l>N+va & l<=N+v & l %in% s & m==l){
        ####come in dispari anche questo 5,6 diversi dal PROF
        matrice[l,m]=-Inf
        #sbaglaito
        # l>N+va & l<=N+v & l==ff[ff$ris==K+1-j,2]+N & m==l
        
        
        
      }
      else if(l>N+va & l<=N+v & l %in% s & m==N){
        matrice[l,m]=Inf
        
        
        
        
      }
      else {
        matrice[l,m]=0
        
        
      }
      # print(tt)
      # m_full_x[,,j]=tt
    }
    
    
  }
  
  return(matrice)  
}
QH_Xdispari<-function(N,K,v,va,j,mm){
  va=v
  ##Funzione che computa QH|xj, quando l'indice della x, j, è dispari. 
  ###prende in input, N,K #di stati rispettivamente di H e X, 
  ### v, va transizioni non allowate totali, e transizioni non allowate 
  ### all'interno del gruppo A (quello dispari.)
  ### j, indice della X, e mm matrice di rappresentazione diretta. 
  
  ### Restituisce una matrice N+v * N+v
  s<-funz_quartopuntoqhlhm_x(K=K,N=N,j=j)
  matrice=matrix(NA,ncol=N+v,nrow=N+v)
  ff=funzioneordinamento(k=K)
  for (l in 1:(N+v)){
    for (m in 1:(N+v)){
      if (l==m & m<N){
        matrice[l,m]=mm[paste0(l,"_",j),paste0(l,"_",j)]
        
      } else if(l==m & m==N){
        tempo=0
        for (u in 1:K){
          tempo<-tempo+ifelse(u%%2==0,mm[paste0(N,"_",j),paste0(1,"_",u)],0)
          ##ciclo che mi porta dentro la somma quando u è pari
          # print(tempo)
          
        }
        matrice[l,m]=mm[paste0(N,"_",j),paste0(N,"_",j)]+tempo
        # print("bbbb")
        # print(matrice[l,m])
        
      } else if(l!=m & l<=N & m<=N){
        matrice[l,m]=mm[paste0(l,"_",j),paste0(m,"_",j)]
        
        
        
        
      } else if(l==N & m>N & m<=N+va ){
        if( ff[ff$ris==abs((m-N)),1]==j){ #ho messo un if dentro poichè altrimenti 
          #dava errore che ci sono casi dove la chiamata dell'if è indefinita (out of bounds)
          a=ff[ff$ris==(m-N),1]
          b=ff[ff$ris==abs((m-N)),2]
          # print(a)
          # print("...")
          # print(b)
          # print("...")
          
          matrice[l,m]=mm[paste0(N,"_",j),paste0(1,"_",b)]
        } else matrice[l,m]=0
        # print("ddddd")
        # print(matrice[l,m])
        
        
      }
      else if(l>N & l<=N+va & l %in% s & m==l){
        ###per inserire i +-Infinito NON FUNZIONA la formula del Prof, ho usato questa formula 
        ###con s, con s funzione AD-Hoc per questo punto. 
        matrice[l,m]=-Inf
        #problema: SBAGLIATO!!!
        # l>N & l<=N+va & l==ff[ff$ris==K+1-j,2]+N & m==l
        
        
        
      }#j 3--->1  j 1 --->2
      ##possibile sol K-j per dispari e K-j+1 per pari 
      else if(l>N & l<=N+va & l %in% s & m==1){
        matrice[l,m]=Inf
        #problema: SBAGLIATO!!!
        
        
        
        
      }
      else {
        matrice[l,m]=0
        
        
      }
      # print(tt)
      # m_full_x[,,j]=tt
    }
    
    
  }
  
  return(matrice)  
}



# fileamalgamato indietro -------------------------------------------------

#creano rispettivamente le CIM QH|x1,...xk e QX|h1,...hn
funzioneAM_QH_X<-function(k,N,matriceam){
  v=funz_v(k)
  matr_full_H=array(NA,c(N+v,N+v,k))
  for (l in 1:k){
    for (i in 1:(N+v)){
      for (j in 1:(N+v)){
        if (i!=j){
          print(c(i,j,l))
          matr_full_H[i,j,l]=matriceam[(l-1)*(N+v)+i,(l-1)*(N+v)+j]
          
        } else if(i==j){
          tempo=0
          for (u in 1:(N+v)){
            tempo<-tempo-ifelse(u!=i,matriceam[(l-1)*(N+v)+i,(l-1)*(N+v)+u],0)
            ##ciclo che mi porta dentro la somma quando u è pari
            # print(tempo)
            print(c(i,j,l))
          }
          matr_full_H[i,j,l]=tempo
        }
        
      }
    }
  }
  
  return(matr_full_H)
}
funzioneAM_QX_H<-function(k,N,matriceam){
  v=funz_v(k)
  matr_full_X=array(NA,c(k,k,N+v))
  for (l in 1:(N+v)){
    for (i in 1:(k)){
      for (j in 1:(k)){
        if (i!=j){
          print(c(i,j,l,N,v))
          matr_full_X[i,j,l]=matriceam[(i-1)*(N+v)+l,(j-1)*(N+v)+l]
          
        } else if(i==j){
          tempo=0
          for (u in 1:(k)){
            print(c(i,j,l,u))
            
            tempo<-tempo-ifelse(u!=i,matriceam[(i-1)*(N+v)+l,(u-1)*(N+v)+l],0)
            ##ciclo che mi porta dentro la somma quando u è pari
            print(tempo)
            
          }
          matr_full_X[i,j,l]=tempo
        }
        
      }
    }
  }
  
  return(matr_full_X)
}
