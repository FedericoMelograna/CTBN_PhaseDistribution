j=1
v=4
N=3
va=2

l=3
m=1
N=3

matrice=matrix(NA,ncol=N+v,nrow=N+v)

for (m in 1:(N+v)){
  if (l==m & m<N){
    matrice[l,m]=mm[paste0(l,"_",j),paste0(l,"_",j)]
    print("aaaa")
  } else if(l==m & m==N){
    tempo=0
    for (u in 1:K){
      tempo<-tempo+ifelse(u%%2==0,mm[paste0(N,"_",j),paste0(1,"_",u)],0)
      print(tempo)
      
    }
    matrice[l,m]=mm[paste0(N,"_",j),paste0(N,"_",j)]+tempo
    print("bbbb")
    print(matrice[l,m])
    
  } else if(l!=m & l<=N & m<=N){
    matrice[l,m]=mm[paste0(l,"_",j),paste0(m,"_",j)]
    print("ccccc")
    print(matrice[l,m])
    
    
  } else if(l==N & m>N & m<=N+va & ff[ff$ris==(m-N),1]==j){
    a=ff[ff$ris==(m-N),1]
    b=ff[ff$prim==a,2]
    matrice[l,m]=mm[paste0(N,"_",j),paste0(1,"_",b)]
    print("ddddd")
    print(matrice[l,m])
    
    
  }
  else if(l>N & l<=N+va & l%%2!=0 & m==l){
    matrice[l,m]=-Inf
    print("eeeee")
    
    
    
    
  }
  else if(l>N & l<=N+va & l%%2!=0 & m==1){
    matrice[l,m]=-Inf
    print("ffff")
  }
  else {
    matrice[l,m]=0
    print("ggggg")
    print(matrice[l,m])
  }
  # print(tt)
  # m_full_x[,,j]=tt
}
matrice

