library(stringr)
###matrice di prova
mm=matrix(c(-1,0.9,0.1,0,0,0,0,0,0,0,0,0,1.2,-2,0.8,0,0,0,0,0,0,0,0,0,
            0.5,1.5,-3,0.5,0,0,0.3,0,0,0.2,0,0,
            0,0,0,-4,3,1,0,0,0,0,0,0,0,0,0,1,-5,4,0,0,0,0,0,0,0.6,0,0,2,4,-8,1,0,0,0.4,0,0,
            0,0,0,0,0,0,-9,7,2,0,0,0,0,0,0,0,0,0,5,-7,2,0,0,0,1,0,0,1,0,0,5,1,-9,1,0,0,
            0,0,0,0,0,0,0,0,0,-5,4,1,0,0,0,0,0,0,0,0,0,1,-2,1,1,0,0,0.5,0,0,0.5,0,0,2,3,-7),ncol=12,byrow=T)
mm


mm=directnames(matriced=mm,k_int=4,N_int = 3)
mm
K=4
N=3
##valori di K e N della matrice di prova

# Matrice bipartitica dato h1 ---------------------------------------------

##Qxixj| H1
  
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
}

QX_h1(K=4,N=3,mm=mm)
#       [,1] [,2] [,3] [,4]
# [1,]  0.0  0.0  0.0  0.0
# [2,]  0.6 -1.6  1.0  0.0
# [3,]  0.0  0.0  0.0  0.0
# [4,]  1.0  0.0  0.5 -1.5


# matrice bipartitica dato hN ---------------------------------------------

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
}
#esempio

QX_hN(K=4,N=3,mm=mm)
#       [,1] [,2] [,3] [,4]
# [1,] -0.7  0.5    0  0.2
# [2,]  0.0  0.0    0  0.0
# [3,]  0.0  1.0   -2  1.0
# [4,]  0.0  0.0    0  0.0
###


# matrice quando 1<l<N ---> matrice vuota di zeri -------------------------

QX_h0<-function(K,N=1,mm=matrix(c(-1,1,0,0,1,-2,1,0,0,0,-0.5,0.5,0.6,0,0.4,-1),ncol=4,byrow=T)
){##matrice vuota di zeri quando l sta tra 1 e N 
  #dovrei cambiare input: non mi serve a niente avere una matrice. Anche N è tralasciabile 
  return(matrice=matrix(0,ncol=K,nrow=K))
  
}
#esempio
QX_h0(K=4,mm=mm)
#       [,1] [,2] [,3] [,4]
# [1,]    0    0    0    0
# [2,]    0    0    0    0
# [3,]    0    0    0    0
# [4,]    0    0    0    0
##
# K=4--> ci sono

# numero trans non legittime ----------------------------------------------

#numero di transizioni non legittime
(v=funz_v(4))
#4 transazioni non legittime
# -->1,3 in A, 2,4 in B 

##transazioni non allowate: 
# 1--->3
# 3--->1
# 2--->4
# 4--->2

# Matrice X dato h quando l>N ---------------------------------------------

#funzione ausiliaria 1

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

funzsupporto(k=5)
# [1] "2_4" "4_2" "1_3" "3_1"
funzsupporto(k=5)
# [1] "2_4" "4_2" "1_3" "1_5" "3_1" "3_5" "5_1" "5_3"

## funzione ausiliaria 2  
funzioneordinamento<-function(k){
  ff=funzsupporto(k)
  z=data.frame(prim=as.numeric(substr(ff,1,1)),sec=as.numeric(substr(ff,3,3)))
  z$terz=ifelse(z$prim>z$sec,1,0)
  z1=z[z$terz==0,];(z1<-z1[order(z1$prim,z1$sec),]);z1$ris=2*(1:nrow(z1)-1)+1
  z2=z[z$terz!=0,]; z2<-z2[order(z2$sec,z2$prim),];z2$ris=2*(1:nrow(z2))#;z2
  zfin=rbind(z1,z2);zfin=zfin[order(zfin$ris),c(1,2,4)]
  return(zfin)

}


# PUNTO IMPORTANTE DA CONTROLLARE CON IL PROFESSORE -----------------------

#QUESTO ORDINAMENTO è CORRETTO?
funzioneordinamento(k=5)

#   prim sec ris
# 3    1   3   1
# 5    3   1   2
# 4    1   5   3
# 7    5   1   4
# 1    2   4   5
# 2    4   2   6
# 6    3   5   7
# 8    5   3   8


###funzione vera e propria
# N=3;l=6;K=4


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
}

#esempio
QX_l(K=4,N=3,l=7)
#       [,1] [,2] [,3] [,4]
# [1,]    0    0    0    0
# [2,]    0    0    0    0
# [3,]    0    0    0    0
# [4,]    0  Inf    0 -Inf
###



# FUNZIONE COMPLETA QX|H --------------------------------------------------


# Funz_tutteh=function(k,N,matrice){
  
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


#esempio
QX_Htotale(k=4,N=3,matrice=mm)


##





# passiamo ora alla fase successiva QH|X ----------------------------------

###due casi diversi, j dispari e j pari


# QH|X j pari -------------------------------------------------------------

##prove iniziali
for (l in 1:(N+v)){
  for (m in 1:N+v){
    print(m);print(l)
    if (l==m & m<N){
      print(l); print(m)
    }}
}

l=1
m=1
N=3

j=3
v=4
N=3
va=2

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
}




QH_Xdispari(N=3,K=4,v=4,va=2,j=3,mm=mm)
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7]
# [1,]   -9    7    2    0    0    0    0
# [2,]    5   -7    2    0    0    0    0
# [3,]    5    1   -7    0    1    0    0
# [4,]  Inf    0    0 -Inf    0    0    0
# [5,]    0    0    0    0    0    0    0
# [6,]    0    0    0    0    0    0    0
# [7,]    0    0    0    0    0    0    0
QH_Xdispari(N=3,K=5,v=8,va=6,j=1,mm=matrice55)
v=funz_va(5)


# QH|X j DISPARI ----------------------------------------------------------
##prime prove
v=4
# l=2
# m=1
N=3
matrice=matrix(NA,ncol=N+v,nrow=N+v)
j=4
v=4
N=3
va=2


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

# DIVERSO DAL PROF #4 PARI ------------------------------------------------

        
        
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





QH_Xpari(N=3,K=4,v=4,va=2,j=4,mm=mm)

# Adesso completa QH|X_j --------------------------------------------------

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


Funz_tutteh(k=4,N=3,matrice=mm)




# PASSO AMALGAMAZIONE BIPARTITICA -----------------------------------------

# aus=Funz_tutteh(k=4,N=3,matrice=mm);aus[aus==Inf]<-10000;aus[aus==-Inf]<--1000;aus
# k=4; m1=matrix(c(-5,2,3,2,-6,4,2,5,-7),byrow=T,ncol=3); m1
# tot=N+v; CC=completaNH_X(K=4,N=tot) ;CC[,,1]%*%aus[,,1]%*%t(CC[,,1])


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

# esempio
bb=amalgamataH_X(k=4,N=3,matrice=mm)
bb[1:5,1:5]
#       [,1] [,2] [,3] [,4] [,5]
# [1,] -1.0  0.0    0  0.0  0.9
# [2,]  0.0 -6.4    0  0.0  0.0
# [3,]  0.0  0.0   -9  0.0  0.0
# [4,]  0.0  0.0    0 -5.5  0.0
# [5,]  1.2  0.0    0  0.0 -2.0

# amalgamazione X|H
# aus2=Funz_tuttex(k=4,N=3); aus2[aus2==Inf]<-10000 ;aus2[aus2==-Inf]<--1000; aus2
# CC2=completaNX_H(K=4,N=tot); CC2[,,1]%*%aus2[,,1]%*%t(CC2[,,1])

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
# esempio
amalgamataX_H(k=4,N=3,matrice=mm)

## voglio vedere se si sovrappongono tra di loro le due, e si lo fanno!

aa=amalgamataX_H(k=4,N=3,matrice=mm)

bb=amalgamataH_X(k=4,N=3,matrice=mm)
a<-matrix(0,ncol=ncol(aa),nrow=nrow(aa))
a[aa!=0]<-1
b<-matrix(0,ncol=ncol(bb),nrow=nrow(bb))
b[bb!=0]<-1
a+b# 
##si sovrappongono!!!


