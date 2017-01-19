FA<-c("16_0","18_0","18_1","20_0","20_4","20_5","22_0","22_6","24_0","24_1")
length<-c(16,18,18,20,20,20,22,22,24,24)
db<-c(0,0,1,0,4,5,0,6,0,1)
fattyAcids<-data.frame(FA=FA,length=length,db=db)

##  Regression based on real data
lipids<-read.csv("lipids.csv")
temp<-lipids
lipids1<-as.matrix(lipids[,2:24])
lipids2<-as.matrix(lipids[,26:31])
linear<-lm(lipids2[,4]~lipids1)
summary(linear)
###

##  Coefficients
C<-c(3,27,0,5,2,3,6,0,0,6,12,6,1,0,1,0,1,0,0,0,0,0,0,0)
H<-c(3,41,0,12,6,6,11,0,1,10,20,10,2,-2,2,-2,2,-2,0,2,-2,-2,-2,5)
N<-c(-1,-1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
O<-c(2,0,0,3,3,5,8,0,3,5,10,8,0,0,0,0,0,0,-1,-1,1,1,1,1)
P<-c(0,0,0,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
S<-c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
elementalCoef<-matrix(c(C,H,N,O,P,S),ncol=6,dimnames = list(c("Glycerol","Cholesterol","Sphingoid","PC","PE","PS","PI","OH","Pi","Glu","Lac","Sulfatide","SN1-chain","SN1-DB","SN2-chain","SN2-DB","SN3-Chain","SN3-DB","(P)","(O)","sn1","sn2","sn3","offset"),NULL))
lipidsOffset<-cbind(lipids1,rep(1,nrow(lipids1)))
lipids2[,1:6]-lipidsOffset%*%elementalCoef


c("Glycerol","Cholesterol","Sphingoid","PC","PE","PS","PI","OH","Pi","Glu","Lac","Sulfatide","Chain","DB","Sph-Chain","Sph-DB","(P)","(O)","sn1","sn2","sn3","offset")

generateLipid<-function(backBone="Glycerol",headGroup="PC",acylSN1="18_0",acylSN2="18_1",acylSN3=NA,modification=NA) {
  composition<-numeric(24)
  if(backBone=="Glycerol") {
    composition[1:3]<-c(1,0,0)
  } else if(backBone=="Cholesterol") {
    composition[1:3]<-c(0,1,0)
    if(!is.na(acylSN1)) {
      headGroup=NA
    }
    acylSN2<-NA
    acylSN3<-NA
  } else if(backBone=="Sphingoid") {
    composition[1:3]<-c(0,0,1)
    acylSN3<-NA
  }
  if(!is.na(headGroup)) {
    if(headGroup=="PC") {
      composition[4:12]<-c(1,0,0,0,0,0,0,0,0)
      acylSN3<-NA
    } else if(headGroup=="PE") {
      composition[4:12]<-c(0,1,0,0,0,0,0,0,0)
      acylSN3<-NA
    } else if(headGroup=="PS") {
      composition[4:12]<-c(0,0,1,0,0,0,0,0,0)
      acylSN3<-NA
    } else if(headGroup=="PI") {
      composition[4:12]<-c(0,0,0,1,0,0,0,0,0)
      acylSN3<-NA
    } else if(headGroup=="OH") {
      composition[4:12]<-c(0,0,0,0,1,0,0,0,0)
      acylSN3<-NA
    } else if(headGroup=="Pi") {
      composition[4:12]<-c(0,0,0,0,0,1,0,0,0)
    } else if(headGroup=="Glu") {
      composition[4:12]<-c(0,0,0,0,0,0,1,0,0)
      acylSN3<-NA
    } else if(headGroup=="Lac") {
      composition[4:12]<-c(0,0,0,0,0,0,0,1,0)
      acylSN3<-NA
    } else if(headGroup=="Sulfatide") {
      composition[4:12]<-c(0,0,0,0,0,0,0,0,1)
      acylSN3<-NA
    } 
  } else {
    composition[4:12]<-c(0,0,0,0,0,0,0,0,0)
  }
  if(backBone=="Glycerol"||(backBone=="Cholesterol"&!is.na(acylSN1))) {
    composition[13]<-fattyAcids$length[which(fattyAcids$FA==acylSN1)]
    composition[14]<-fattyAcids$db[which(fattyAcids$FA==acylSN1)]
    if(!is.na(acylSN2)) {
      composition[15]<-fattyAcids$length[which(fattyAcids$FA==acylSN2)]
      composition[16]<-fattyAcids$db[which(fattyAcids$FA==acylSN2)]
    }
    if(!is.na(acylSN3)) {
      composition[17]<-fattyAcids$length[which(fattyAcids$FA==acylSN3)]
      composition[18]<-fattyAcids$db[which(fattyAcids$FA==acylSN3)]
    }
  } else if(backBone=="Sphingoid") {
    composition[13]<-fattyAcids$length[which(fattyAcids$FA==acylSN1)]
    composition[14]<-fattyAcids$db[which(fattyAcids$FA==acylSN1)]
    if(!is.na(acylSN2)) {
      composition[15]<-fattyAcids$length[which(fattyAcids$FA==acylSN2)]
      composition[16]<-fattyAcids$db[which(fattyAcids$FA==acylSN2)]
    }
  }
  if(!is.na(modification)) {
    if(modification=="(P)") {
      composition[19:20]<-c(1,0)
    } else if(modification=="(O)") {
      composition[19:20]<-c(0,1)
    }
  }
  if(!is.na(acylSN1)) {
    composition[21]<-1
  }
  if(!is.na(acylSN2)) {
    composition[22]<-1
  }
  if(!is.na(acylSN3)) {
    composition[23]<-1
  }
  composition[24]<-1
  return(composition)
}

generateLipid("Sphingoid",headGroup="OH",acylSN1="18_1",acylSN2 = NA)

generateLipid("Sphingoid",headGroup = "OH",acylSN1 = "18_1",acylSN2 = "18_0")%*%elementalCoef
