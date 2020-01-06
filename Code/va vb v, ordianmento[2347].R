
# Funzione che computa il numero V di transizioNiX_H non legittime ----------



# Funzione di transizione non allowate ------------------------------------
#va--> non allowate nel gruppo A
funz_va<-function(k_int){
  ###OVVERO LA FORMULA è: k/2*(k/2 -1) con approx per eccesso
  return(ceiling(k_int/2)*(ceiling(k_int/2)-1))
}

# funz_va(3)
# funz_va(1)
# funz_va(2)

funz_vb<-function(k_int){
  ###OVVERO LA FORMULA è: k/2*(k/2 -1) con approx per DIFETTO. 
  return(floor(k_int/2)*(floor(k_int/2)-1))
}

funz_v<-function(k_int){
  ### combina semplicemente le due sopra, restituisce il numero totale v
  ### di transizioni non allowate
  return(funz_vb(k_int)+funz_va(k_int))
}

# funz_v(3)
# funz_v(4)






# Amalgamazione, matrice Ns|ci --------------------------------------------
# N=3
# K=5
#FUNZIONE ORDINAMENTO CONGIUNTO 

# DA CAMBIARE!!!!!!!!!!!!!! -----------------------------------------------
#cambiata

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
}
#ordinamenti parziali invece di X e H 
ordinamentoX<-function(K) return(paste0("K",seq(1,K,by=1)))
ordinamentoH<-function(N) return(paste0("N",seq(1,N,by=1)))

# ESEMPI
ordinamentoX_H(3,10)
ord=ordinamentoX_H(N=3,K=5)
ordX=ordinamentoX(5)
ordH=ordinamentoH(3)
# FINE EX

###################k


# Creazione matrice N di x|h1 e H|x1 --------------------------------------

###esempio con N=3,K=2,i=2

# K=2;N=3

# funzione singola matricina NX|hi -----------------------------------------

####k=#stati di X
####N=#stati di H
####i=quale N stiato calcolando
####ord=ordinamento congiunto di X e H con precedenza a X
####ordX=ordinamento solamente di X
####ordH=ordinamento di H

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
}


# funzione singola matricina N H|xi -----------------------------------------


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
}


# esempi
singolaNH_xi(K=2,N=3,i=2,ord= ordinamentoX_H(N=3,K=2),
             ordH=ordinamentoH(N=3)
             ,ordX=ordinamentoX(K=2))

singolaNX_hi(2,3,1,ord=ord,ordH = ordH,ordX=ordX)
# funzia!




# Adesso ci creiamo tutte le funzioni NX|hi per ogni i --------------------

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
}


#esempio
completaNH_X(2,3)




# Inizio di come dovrebbe essere amalgamazione per bipartitica ------------


m1=matrix(c(-5,2,3,2,-6,4,2,5,-7),byrow=T,ncol=3)
m1
(CC=completaNH_X(2,3) )
CC[,,1]%*%m1%*%t(CC[,,1])


##verra fatto alla fine della bipartitica
