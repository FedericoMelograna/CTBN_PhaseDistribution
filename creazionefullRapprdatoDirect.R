
###prima matrice di esempio, una matrice N=2,k=2
#IN RAPPRESENTAZIONE DIRETTA. 
 md=matrix(c(-1,1,0,0,1,-2,1,0,0,0,-0.5,0.5,0.6,0,0.4,-1),ncol=4,byrow=T)
 md
#primo 2*2= interna a S1, 2*2 a dx, va da S1 a S2, 
#e al contrario seconda riga


# RAPPRESENTAZIONE FULL: DA DIRECT A FULL ---------------------------------

 
 
# prima parte, funzione che crea singola Qx|h1 ----------------------------
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


##funzia

singolax(k_int = 2,j_int=1)


# seconda parte: funzione che crea matrice vuota X|hi ---------------------


vuota=function(k_int){
  m_int=matrix(0,ncol=k_int,nrow=k_int)
  return(m_int)
  
}
vuota(6)


# Terza parte: creo tutte le matrici Qx|H ---------------------------------

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


}


###funzia
##parametri di funz_tuttex, K=#di stati s1,..sk e N=#di stati nascosti per ogni si 

Funz_tuttex(k=4,N=3)
#funz tuttex che ci serve per le matrici X della rappresentazione full dell'esempio del prog

Funz_tuttex(k=2,N=2)
#per il mio esempio della 2*2

##NB le matrici X|Hi non sono influenzati dai dati, dalle intenstià interne alle matrici, 
#ma sono bensì fisse, dati numero di hidden state e di shallow state. 




###adesso proviamo a fare h dato x

#mia matrice di esempio
md=matrix(c(-1,1,0,0,1,-2,1,0,0,0,-0.5,0.5,0.6,0,0.4,-1),ncol=4,byrow=T)
md
#quella del prof
mm=matrix(c(-1,0.9,0.1,0,0,0,0,0,0,0,0,0,1.2,-2,0.8,0,0,0,0,0,0,0,0,0,
            0.5,1.5,-3,0.5,0,0,0.3,0,0,0.2,0,0,
            0,0,0,-4,3,1,0,0,0,0,0,0,0,0,0,1,-5,4,0,0,0,0,0,0,0.6,0,0,2,4,-8,1,0,0,0.4,0,0,
            0,0,0,0,0,0,-9,7,2,0,0,0,0,0,0,0,0,0,5,-7,2,0,0,0,1,0,0,1,0,0,5,1,-9,1,0,0,
            0,0,0,0,0,0,0,0,0,-5,4,1,0,0,0,0,0,0,0,0,0,1,-2,1,1,0,0,0.5,0,0,0.5,0,0,2,3,-7),ncol=12,byrow=T)
mm
# quinta parte, funzione che rinomina colonne e righe di direct ----------------------------
directnames=function(k_int,N_int,matriced){
  vector=c()
  for (i in 1:k_int){
    for (j in 1:N_int){
      vector=c(vector,paste0(j,"_",i))
      
    }
  }
  print(vector)
  rownames(matriced)=vector;colnames(matriced)=vector  
  return(matriced)
}

md=directnames(matriced = md,k_int=2,2)
md
mm=directnames(matriced=mm,k_int=4,N_int = 3)
mm
# singolah, singola QH|x1 -------------------------------------------------


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
}

##esempi di singola H|xi per  md
singolah(2,2,1,matriced = md)
singolah(2,2,2,matriced = md)
##esempi di singola H|xi per  mm, matrice del prof
singolah(k_int=4,N_int=3,j_int=1,matriced=mm)
singolah(k_int=4,N_int=3,j_int=4,matriced=mm)
#j va da 1 a k_int

# insieme delle h  --------------------------------------------------------



Funz_tutteh=function(k,N,matrice){
  m_full_x=array(NA,c(N*k,N*k,k))
  
  for (j in 1:(k)){ #N*k)
      m_full_x[,,j]=singolah(k_int=k,N_int=N,j_int=j,matriced=matrice)
  }
  return(m_full_x)
      
}



Funz_tutteh(k=4,N=3,matrice=mm)





#############fine full rappresentation


##inizio rappresentazione bipartita. 
