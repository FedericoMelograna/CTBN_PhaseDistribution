
# Mappatura da amalgamata a singole matricini -----------------------------

##funzione che si occupa di rimappare l'amalgaata in QH|X

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



matrice33

N=3
v=2
k=3

salvataggio=funzioneAM_QH_X(k=k,N=N,matriceam=d33)
salvataggio
###come utilizzarla concretamente


funzioneAM_QX_H(k,N,d33)

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

###funzia per adesso!!
i=1;j=1;l=4;u=5;
DF55[(i-1)*(N+v)+l,(u-1)*(N+v)+l]
# Proviamo con altre matrici ----------------------------------------------

k=5;N=3

salvataggio53=funzioneAM_QX_H(k=k,N=N,matriceam=DF53)
t1[t1==Inf]<-10000
t1[t1==-Inf]<--10000
sum((salvataggio53-t1)^2)
##taaac 

#vediamo se funziona anche l'inversa
salvataggio53H_X=funzioneAM_QH_X(k=k,N=N,matriceam=DF53)
t2[t2==Inf]<-10000
t2
t2[t2==-Inf]<--10000
i=3
sum(salvataggio53H_X[]-t2[])
salvataggio53H_X[,,3]
t2[,,3]

salvataggioprof=funzioneAM_QX_H(k=4,N=3,matriceam=df)
salvataggioprof

salvataggioprofinversa=funzioneAM_QH_X(k=4,N=3,matriceam=df)
salvataggioprofinversa

##con la 4*3 funziona
