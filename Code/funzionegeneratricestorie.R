
# generazione di una storia da una matrice --------------------------------

#esempio con Qsempl

# vediamo quanto ci mette a generare --------------------------------------

Qsempl
start=Sys.time()
stato_iniz<-rmultinom(10000,1,c(0.25,0.1,0.5,0.15))
end=Sys.time()
end-start
stato_iniz
rowSums(stato_iniz)
rmultinom(n=1,size=)
start=Sys.time()
ss=sample(x =c(1,2,3,4) ,size =10000 ,replace = T,prob = c(0.25,0.1,0.6,0.15))
end=Sys.time()
end-start
table(ss)

rexp(1,rate=1000)

Qsecondo=Qsempl
diag(Qsecondo)=0
stato_iniz=sample(x =c(1,2,3,4) ,size =1 ,replace = T,prob = c(0.25,0.25,0.25,0.25))
tempo_iniz=rexp(1,rate=abs(Qsempl[stato_iniz,stato_iniz]))
stato
table(sample(x =c(1,2,3,4) ,size =40000 ,replace = T, prob = c(Qsecondo[stato_iniz,]))
)

Qsempl

x <- sample(x =c(1,2,3,4) ,size =1 ,replace = T,prob = c(0.25,0.25,0.25,0.25))
collezionestorie<-data.frame("statoxh"=NA,"inizio"=NA,"fine"=NA)
i=1
tempo=0
Qsecondo=Qsempl
diag(Qsecondo)=0
repeat {
  c1=x
  c2=tempo
  tempo_i=rexp(1,rate=abs(Qsempl[x,x]))
  x<-sample(x =c(1,2,3,4) ,size =1 ,replace = T, prob = c(Qsecondo[x,]))
  tempo=tempo+tempo_i
  c3=tempo
  collezionestorie<-rbind(collezionestorie,c(c1,c2,c3))
  if (tempo >= 100){
    break
  }
  i=i+1
}

collezionestorie$fine[nrow(collezionestorie)]<-100
collezionestorie$stato<-ifelse(collezionestorie$statoxh==1 |collezionestorie$statoxh==2,1,2 )
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
storiaqsempl$fine<-c(storiaqsempl$inizio[-1],100)
storiaqsempl
collezionestorie
# 218       1 94.21474779  95.07565786     1
# 219       3 95.07565786  95.16463730     2
# 220       1 95.16463730  96.04145715     1
# 221       2 96.04145715  98.41094948     1
# 222       4 98.41094948  99.58567047     2
# 223       3 99.58567047  99.62577399     2
# 224       4 99.62577399 100.00000000     2


statixh=seq(from=1,to=28,1)
x <- sample(x =statixh ,size =1000 ,replace = T,prob = rep(1/28,28))
table(x)



# qui la funzione importante! ---------------------------------------------



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


ggg=generastorie(matrice =Qsempl ,tmax =100 ,statix =2 ,Nv =2 ,N =2 )
ggg$storiafinale
ggg$storiavera



# da fare domani: mettere delle transienze per matrice df -----------------

#probabilmente va gia bene cosi
