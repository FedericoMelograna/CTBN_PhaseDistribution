
# Strutturazione alternativa QH|X -----------------------------------------


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



