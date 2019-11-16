### Matrici di prova e di lavoro 

m1=matrix(c(-5,2,3,2,-6,4,2,5,-7),byrow=T,ncol=3)
m1
#inutile

md=matrix(c(-1,1,0,0,1,-2,1,0,0,0,-0.5,0.5,0.6,0,0.4,-1),ncol=4,byrow=T)
md
#primo 2*2= interna a S1, 2*2 a dx, va da S1 a S2, 
#e al contrario seconda riga

md=directnames(matriced = md,k_int=2,2)
md
##matrice PROF
mm=matrix(c(-1,0.9,0.1,0,0,0,0,0,0,0,0,0,1.2,-2,0.8,0,0,0,0,0,0,0,0,0,
            0.5,1.5,-3,0.5,0,0,0.3,0,0,0.2,0,0,
            0,0,0,-4,3,1,0,0,0,0,0,0,0,0,0,1,-5,4,0,0,0,0,0,0,0.6,0,0,2,4,-8,1,0,0,0.4,0,0,
            0,0,0,0,0,0,-9,7,2,0,0,0,0,0,0,0,0,0,5,-7,2,0,0,0,1,0,0,1,0,0,5,1,-9,1,0,0,
            0,0,0,0,0,0,0,0,0,-5,4,1,0,0,0,0,0,0,0,0,0,1,-2,1,1,0,0,0.5,0,0,0.5,0,0,2,3,-7),ncol=12,byrow=T)
mm

mm=directnames(matriced=mm,k_int=4,N_int = 3)
mm

aa=amalgamataH_X(k=4,N=3,matrice=mm)
bb=amalgamataX_H(k=4,N=3,matrice=mm)
df=as.data.frame(aa+bb)
print(df)
#write.csv2(df,"amalgamataprofcorretta.csv")
#qx|h
QX_Htotale(k=4,N=3,matrice = mm)
#QH|X
Funz_tutteh(k=4,N=3,matrice=mm)



###matrice33
matrice33<-matrix(c(-1,0.7,0.3,0,0,0,0,0,0,
                    1.8,-3,1.2,0,0,0,0,0,0,
                    3,6,-10,0.5,0,0,0.5,0,0,
                    0,0,0,-2,1.4,0.6,0,0,0,
                    0,0,0,2.7,-5,2.3,0,0,0,
                    2.1,0,0,4,7,-15,1.9,0,0,
                    0,0,0,0,0,0,-0.5,0.26,0.24,
                    0,0,0,0,0,0,0.1,-0.3,0.2,
                    1.5,0,0,3.5,0,0,13,12,-30),ncol=9,byrow=T)
matrice33<-directnames(matriced=matrice33,k_int=3,N_int = 3)
matrice33

a1=amalgamataH_X(k=3,N=3,matrice=matrice33)
b1=amalgamataX_H(k=3,N=3,matrice=matrice33)
d33=a1+b1
#write.csv2(d33,"amalgamata33.csv")

QX_Htotale(k=3,N=3,matrice = matrice33)
Funz_tutteh(k=3,N=3,matrice=matrice33)

#write.csv2(QX_Htotale(k=3,N=3,matrice = matrice33),"QX_H33")
#write.csv2(Funz_tutteh(k=3,N=3,matrice=matrice33),"QH_X33")



##matriceBIG
matrice53<-matrix(c(-1,0.3,0.7,0,0,0,0,0,0,0,0,0,0,0,0,
                    1.4,-3,1.6,0,0,0,0,0,0,0,0,0,0,0,0,
                    3,6,-10,0.5,0,0,0.5,0,0,0,0,0,0,0,0,
                    0,0,0,-2,1.9,0.1,0,0,0,0,0,0,0,0,0,
                    0,0,0,2.7,-5.1,2.4,0,0,0,0,0,0,0,0,0,
                    2.1,0,0,4,6.5,-15.2,1.3,0,0,0.6,0,0,0.7,0,0,
                    0,0,0,0,0,0,-0.5,0.26,0.24,0,0,0,0,0,0,
                    0,0,0,0,0,0,0.1,-0.3,0.2,0,0,0,0,0,0,
                    1.5,0,0,3.5,0,0,13,12,-31,0,0,0,1,0,0,
                    0,0,0,0,0,0,0,0,0,-8,4,4,0,0,0,
                    0,0,0,0,0,0,0,0,0,9,-16,7,0,0,0,
                    2,0,0,1,0,0,0,0,0,3,0.5,-7,0.5,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,-3,1.6,1.4,
                    0,0,0,0,0,0,0,0,0,0,0,0,3.5,-7,3.5,
                    0,0,0,5,0,0,1,0,0,0,0,0,11,1,-18),ncol=15,byrow=T)
matrice53<-directnames(matriced=matrice53,k_int=5,N_int = 3)
matrice53
sum(matrice53)
round(rowSums(matrice53),2)
matrice53
tt53=amalgamataH_X(k=5,N=3,matrice=matrice53)
tt53<-as.data.frame(tt53)
tt53
# write.csv(tt55,"matriceamalgamataH_X55")
tt53bis=amalgamataX_H(k=5,N=3,matrice=matrice53)
tt53bis<-as.data.frame(tt53bis)
tt53bis
# write.csv(tt56,"matriceamalgamataX_H56")
DF53=tt53+tt53bis
#write.csv2(DF55,"Amalgamata55.csv")
(t1=QX_Htotale(k=5,N=3,matrice = matrice53))
(t2=Funz_tutteh(k=5,N=3,matrice=matrice53))
# 1_2  0.0  0.0   0.0 |-2.0  1.9   0.1 | 0.0  0.00   0.00  0.0   0.0   0  0.0  0.0   0.0
# 2_2  0.0  0.0   0.0 | 2.7 -5.1   2.4 | 0.0  0.00   0.00  0.0   0.0   0  0.0  0.0   0.0
# 3_2  2.1  0.0   0.0 | 4.0  6.5 -15.2 |1 .3  0.00   0.00  0.6   0.0   0  0.7  0.0   0.0
