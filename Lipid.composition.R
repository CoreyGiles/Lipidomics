library(reshape2)

FA<-    c("16_0","16_1","18_0","18_1","18_2","18_3","20_0","20_4","20_5","22_0","22_1","22_5","22_6","24_0","24_1","26_0","26_1")
length<-c( 16,    16,    18,    18,    18,    18,    20,    20,    20,    22,    22,    22,    22,    24,    24,    26,    26)
db<-    c( 0,     1,     0,     1,     2,     3,     0,     4,     5,     0,     1,     5,     6,     0,     1,     0,     1)
fattyAcids<-data.frame(FA=FA,length=length,db=db)

headGroup<-c("PC","PE","PS","PI","OH","TG","Pi","Glu","Lac","Sulfatide")

##  Regression based on lipid data
#lipids<-read.csv("lipids.csv")
#temp<-lipids
#lipids1<-as.matrix(lipids[,2:24])
#lipids2<-as.matrix(lipids[,26:31])
#linear<-lm(lipids2[,4]~lipids1)
#summary(linear)
#lipidsOffset<-cbind(lipids1,rep(1,nrow(lipids1)))
#lipids2[,1:6]-lipidsOffset%*%elementalCoef
###

##  Coefficients
#   Gly Chol Sph PC  PE  PS  PI  OH  Pi Glu Lac Sul SN1  DB SN2  DB SN3  DB (P) (O) sn1 sn2 sn3 offset
C<-c( 3, 27,  0,  5,  2,  3,  6,  0,  0,  6, 12,  6,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0)
H<-c( 3, 41,  0, 12,  6,  6, 11,  0,  1, 10, 20, 10,  2, -2,  2, -2,  2, -2,  0,  2, -2, -2, -2,  5)
N<-c(-1, -1,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1)
O<-c( 2,  0,  0,  3,  3,  5,  8,  0,  3,  5, 10,  8,  0,  0,  0,  0,  0,  0, -1, -1,  1,  1,  1,  1)
P<-c( 0,  0,  0,  1,  1,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
S<-c( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
elementalCoef<-matrix(c(C,H,N,O,P,S),ncol=6,dimnames = list(c("Glycerol","Cholesterol","Sphingoid","PC","PE","PS","PI","OH","Pi","Glu","Lac","Sulfatide","SN1-chain","SN1-DB","SN2-chain","SN2-DB","SN3-Chain","SN3-DB","(P)","(O)","sn1","sn2","sn3","offset"),NULL))

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

generateLipid("Sphingoid",headGroup ="OH",acylSN1="18_1",acylSN2="18_0")%*%elementalCoef


lipidSet<-matrix(0,nrow=4,ncol=nrow(fattyAcids),dimnames=list(c("Glycerol","SphingoidBase","SphingoidFA","Cholesterol"),fattyAcids$FA))
lipidSet[c(1,4),]<-1
lipidSet[2,]<-fattyAcids$FA%in%c("18_0","18_1")
lipidSet[3,]<-fattyAcids$FA%in%c("16_0","18_0","18_1","20_0","22_0","22_1","24_0","24_1","26_0","26_1")

headGroupSet<-matrix(0,nrow=3,ncol=length(headGroup),dimnames = list(c("Glycerol","Sphingoid","Cholesterol"),headGroup))
headGroupSet[1,c(1,2,3,4,5,6,7)]<-1
headGroupSet[2,c(1,5,7,8,9,10)]<-1
headGroupSet[3,5]<-1

lysoGroupSet<-matrix(0,nrow=3,ncol=length(headGroup),dimnames = list(c("Glycerol","Sphingoid","Cholesterol"),headGroup))
lysoGroupSet[1,c(1,2,4)]<-1
lysoGroupSet[2,c(5,7)]<-1

vinylEtherGroupSet<-matrix(0,nrow=3,ncol=length(headGroup),dimnames = list(c("Glycerol","Sphingoid","Cholesterol"),headGroup))
vinylEtherGroupSet[1,c(1,2)]<-1

etherGroupSet<-matrix(0,nrow=3,ncol=length(headGroup),dimnames = list(c("Glycerol","Sphingoid","Cholesterol"),headGroup))
etherGroupSet[1,c(1,2)]<-1

lipidMatrix<-list()
lipidMatrix$lipidSet<-lipidSet
lipidMatrix$headGroupSet<-headGroupSet
lipidMatrix$lysoGroupSet<-lysoGroupSet
lipidMatrix$vinylEtherGroupSet<-vinylEtherGroupSet
lipidMatrix$etherGroupSet<-etherGroupSet

glyceroFattyAcidMatrix<-t(lipidSet[1,,drop=F])%*%lipidSet[1,,drop=F]
glyceroFattyAcidMatrix<-lower.tri(set,diag=T)
colnames(glyceroFattyAcidMatrix)<-fattyAcids$FA
rownames(glyceroFattyAcidMatrix)<-fattyAcids$FA

############################# Need special function for TG ###################################
######################### (O) and (P) modifications cannot be applied to lyso species currently.

glyceroFattyAcidMatrix<-melt(glyceroFattyAcidMatrix)
glyceroFattyAcidMatrix<-glyceroFattyAcidMatrix[glyceroFattyAcidMatrix$value==TRUE,]
colnames(glyceroFattyAcidMatrix)<-c("SN2","SN1","headGroup")
eval(parse(text=paste("glyceroLipidMatrix<-rbind(",paste(rep("glyceroFattyAcidMatrix",sum(headGroupSet[1,])),collapse=","),")",sep="")))
glyceroLipidMatrix$headGroup<-rep(colnames(headGroupSet)[which(headGroupSet[1,]==1)],each=nrow(glyceroLipidMatrix)/sum(headGroupSet[1,]))
glyceroLipidMatrix$SN3<-rep(NA,nrow(glyceroLipidMatrix))
glyceroLipidMatrix$modification<-rep(NA,nrow(glyceroLipidMatrix))

for(i in which(vinylEtherGroupSet[1,]==1)){
  glyceroLipidMatrix<-rbind(glyceroLipidMatrix,data.frame(SN1=glyceroFattyAcidMatrix$SN1,SN2=glyceroFattyAcidMatrix$SN2,headGroup=rep(colnames(vinylEtherGroupSet)[i],nrow(glyceroFattyAcidMatrix)),SN3=rep(NA,nrow(glyceroFattyAcidMatrix)),modification=rep('(P)',nrow(glyceroFattyAcidMatrix))))
}
for(i in which(etherGroupSet[1,]==1)){
  glyceroLipidMatrix<-rbind(glyceroLipidMatrix,data.frame(SN1=glyceroFattyAcidMatrix$SN1,SN2=glyceroFattyAcidMatrix$SN2,headGroup=rep(colnames(etherGroupSet)[i],nrow(glyceroFattyAcidMatrix)),SN3=rep(NA,nrow(glyceroFattyAcidMatrix)),modification=rep('(O)',nrow(glyceroFattyAcidMatrix))))
}
for(i in which(lysoGroupSet[1,]==1)){
  glyceroLipidMatrix<-rbind(glyceroLipidMatrix,data.frame(SN1=fattyAcids$FA,SN2=rep(NA,nrow(fattyAcids)),headGroup=rep(colnames(lysoGroupSet)[i],nrow(fattyAcids)),SN3=rep(NA,nrow(fattyAcids)),modification=rep(NA,nrow(fattyAcids))))
}


printLipids<-function(X) {
  if(!is.na(X$SN2)) {
    if(is.na(X$modification)) {
      print(paste(X$headGroup," ",X$SN1,"/",X$SN2,sep=""))
    } else {
      print(paste(X$headGroup,X$modification," ",X$SN1,"/",X$SN2,sep=""))
    }
  } else {
    if(is.na(X$modification)) {
      print(paste(X$headGroup," ",X$SN1,"/0_0",sep=""))
    } else {
      print(paste(X$headGroup,X$modification," ",X$SN1,"/0_0",sep=""))
    }
  }
}
printLipids(glyceroLipidMatrix[1100,])

sphingoLipidMatrix<-t(lipidSet[2,,drop=F])%*%lipidSet[3,,drop=F]
sphingoLipidMatrix<-melt(sphingoLipidMatrix)
sphingoLipidMatrix<-sphingoLipidMatrix[which(sphingoLipidMatrix$value==1),]


library(R6)
simLipid<-R6Class("Lipid",
  private=list(
    acylSN1="18_0",
    acylSN2="18_1",
    acylSN3=NA,
    .headGroup="PC",
    backBone="Glycerol",
    modification=NA,
    isVinylEther=FALSE,
    isEther=FALSE,
    isLyso=FALSE,
    composition=c(1,0,0,1,0,0,0,0,0,0,0,0,18,0,18,1,0,0,0,0,1,1,0,1),
    elements=NA,
    isotopeDistribution=NA,
    
    generateComposition=function() {
      private$composition<-numeric(24)
      if(private$backBone=="Glycerol") {
        private$composition[1:3]<-c(1,0,0)
      } else if(private$backBone=="Cholesterol") {
        private$composition[1:3]<-c(0,1,0)
        if(!is.na(private$acylSN1)) {
          private$.headGroup=NA
        }
        private$acylSN2<-NA
        private$acylSN3<-NA
      } else if(private$backBone=="Sphingoid") {
        private$composition[1:3]<-c(0,0,1)
        private$acylSN3<-NA
      }
      if(!is.na(private$.headGroup)) {
        if(private$.headGroup=="PC") {
          private$composition[4:12]<-c(1,0,0,0,0,0,0,0,0)
          private$acylSN3<-NA
        } else if(private$.headGroup=="PE") {
          private$composition[4:12]<-c(0,1,0,0,0,0,0,0,0)
          private$acylSN3<-NA
        } else if(private$.headGroup=="PS") {
          private$composition[4:12]<-c(0,0,1,0,0,0,0,0,0)
          private$acylSN3<-NA
        } else if(private$.headGroup=="PI") {
          private$composition[4:12]<-c(0,0,0,1,0,0,0,0,0)
          private$acylSN3<-NA
        } else if(private$.headGroup=="OH") {
          private$composition[4:12]<-c(0,0,0,0,1,0,0,0,0)
          private$acylSN3<-NA
        } else if(private$.headGroup=="Pi") {
          private$composition[4:12]<-c(0,0,0,0,0,1,0,0,0)
        } else if(private$.headGroup=="Glu") {
          private$composition[4:12]<-c(0,0,0,0,0,0,1,0,0)
          private$acylSN3<-NA
        } else if(private$.headGroup=="Lac") {
          private$composition[4:12]<-c(0,0,0,0,0,0,0,1,0)
          private$acylSN3<-NA
        } else if(private$.headGroup=="Sulfatide") {
          private$composition[4:12]<-c(0,0,0,0,0,0,0,0,1)
          private$acylSN3<-NA
        } 
      } else {
        private$composition[4:12]<-c(0,0,0,0,0,0,0,0,0)
      }
      if(private$backBone=="Glycerol"||(private$backBone=="Cholesterol"&!is.na(private$acylSN1))) {
        private$composition[13]<-fattyAcids$length[which(fattyAcids$FA==private$acylSN1)]
        private$composition[14]<-fattyAcids$db[which(fattyAcids$FA==private$acylSN1)]
        if(!is.na(private$acylSN2)) {
          private$composition[15]<-fattyAcids$length[which(fattyAcids$FA==private$acylSN2)]
          private$composition[16]<-fattyAcids$db[which(fattyAcids$FA==private$acylSN2)]
        }
        if(!is.na(private$acylSN3)) {
          private$composition[17]<-fattyAcids$length[which(fattyAcids$FA==private$acylSN3)]
          private$composition[18]<-fattyAcids$db[which(fattyAcids$FA==private$acylSN3)]
        }
      } else if(private$backBone=="Sphingoid") {
        private$composition[13]<-fattyAcids$length[which(fattyAcids$FA==private$acylSN1)]
        private$composition[14]<-fattyAcids$db[which(fattyAcids$FA==private$acylSN1)]
        if(!is.na(acylSN2)) {
          private$composition[15]<-fattyAcids$length[which(fattyAcids$FA==private$acylSN2)]
          private$composition[16]<-fattyAcids$db[which(fattyAcids$FA==private$acylSN2)]
        }
      }
      if(!is.na(private$modification)) {
        if(private$modification=="(P)") {
          private$composition[19:20]<-c(1,0)
        } else if(private$modification=="(O)") {
          private$composition[19:20]<-c(0,1)
        }
      }
      if(!is.na(private$acylSN1)) {
        private$composition[21]<-1
      }
      if(!is.na(private$acylSN2)) {
        private$composition[22]<-1
      }
      if(!is.na(private$acylSN3)) {
        private$composition[23]<-1
      }
      private$composition[24]<-1
    },
    getIsotopeDistribution=function(value) {
      output.combined<-matrix(c(0,1),nrow=1,ncol=2)
      for(p in 1:length(elements)) {
        n<-value[p]
        mass1<-natural.isotopic.abundance[p,2]
        mass2<-natural.isotopic.abundance[p,3]
        abundance<-natural.isotopic.abundance[p,4]
        mass.dist<-matrix(0,nrow=n+1,ncol=2)
        for(k in 0:n) {
          mass.dist[k+1,1]<-mass1*k+(n-k)*mass2
          mass.dist[k+1,2]<-dbinom(k,n,abundance)
        }
        mass.dist<-mass.dist[which(mass.dist[,2]>0.0001),,drop=FALSE]
        eval(parse(text=paste("output.combined<-rbind(",paste(rep("output.combined",nrow(mass.dist)),collapse=","),")",sep="")))
        output.combined[,1]<-output.combined[,1]+rep(mass.dist[,1],each=nrow(output.combined)/nrow(mass.dist))
        output.combined[,2]<-output.combined[,2]*rep(mass.dist[,2],each=nrow(output.combined)/nrow(mass.dist))
      }
      output.combined<-output.combined[which(output.combined[,2]>0.0001),,drop=FALSE]
      #output.combined[,2]<-output.combined[,2]*100/max(output.combined[,2])
      output.combined<-output.combined[order(output.combined[,1]),]
      private$isotopeDistribution<-output.combined
    }
  ),
  public=list(
    initialize=function(...) {
      for(item in list(...)) {
        print(item)
      }
    },
    getIsotope=function() {
      private$generateComposition()
      private$elements<-private$composition%*%elementalCoef
      private$getIsotopeDistribution(private$elements)
      print(private$isotopeDistribution)
      plot(private$isotopeDistribution,type="h")
    }
  ),
  active=list(
    SN1=function(value) {
      if(missing(value)) {
        private$acylSN1
      } else if(value%in%fattyAcids$FA) {
        private$acylSN1<-value
      } else {
        print("Fatty acid not recognized")
      }
    },
    SN2=function(value) {
      if(missing(value)) {
        private$acylSN2
      } else if(value%in%fattyAcids$FA) {
        if(private$isLyso) {
          print("Cannot set SN2 for lyso-lipid")
        } else {
          private$acylSN2<-value
        }
      } else {
        print("Fatty acid not recognized")
      }
    },
    SN3=function(value) {
      if(missing(value)) {
        private$acylSN3
      } else if(value%in%fattyAcids$FA) {
        if(private$.headGroup=="TG") {
          private$acylSN3<-value
        } else {
          print("Cannot set SN3: Lipid is not a triacylglycerol")
        }
      } else {
        print("Fatty acid not recognized")
      }
    },
    headGroup=function(value) {
      if(missing(value)) {
        private$.headGroup
      }
    }
  )
)

PC34_1<-simLipid$new()
PC34_1$getIsotope()

