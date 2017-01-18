#FA<-c("16_0","18_0","18_1","20_0","20_4","20_5","22_0","22_6","24_0","24_1")
#length<-c(16,18,18,20,20,20,22,22,24,24)
#db<-c(0,0,1,0,4,5,0,6,0,1)
#fatty.acids<-data.frame(FA=FA,length=length,db=db)
#BB<-c("glycerol","sphingoid")
#HG<-c("OH","PC","PE","PS","PI")

lipids<-read.csv("lipids.csv")
temp<-lipids

lipids1<-as.matrix(lipids[,2:10])
lipids2<-as.matrix(lipids[,12:16])
linear<-lm(lipids2[,5]~lipids1)
summary(linear)

C<-lipids$PC*5+lipids$PE*2+lipids$PS*3+lipids$PI*6+lipids$C+3
H<-lipids$PC*12+lipids$PE*6+lipids$PS*6+lipids$PI*11+lipids$C*2-lipids$DB*2-lipids$sn2*2-lipids$sn3*2+6
N<-lipids$PC*1+lipids$PE*1+lipids$PS*1
O<-lipids$PC*3+lipids$PE*3+lipids$PS*5+lipids$PI*8+lipids$sn2*1+lipids$sn3*1+4
P<-lipids$PC*1+lipids$PE*1+lipids$PS*1+lipids$PI*1
temp$C.1<-temp$C.1-C
temp$H<-temp$H-H
temp$N<-temp$N-N
temp$O<-temp$O-O
temp$P<-temp$P-P
