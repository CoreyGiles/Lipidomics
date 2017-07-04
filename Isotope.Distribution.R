##  Generate isotope distribution of lipids
##  Input is the number of C,H,N,O,P atoms of the molecule
##  Output is the isotope distribution based on natural elemental isotope abundance
##  Isotope abundance is calculated using binomial distributions


# Currently, isotope abundance are normalized to most abundant ion
# Currently, only returns the most abundant ions to 0.01 %

elements<-c("C","H","N","O","P")
mass1<-c(12,1.007825,14.003074,15.994915,30.973762)
mass2<-c(13.003355,2.014102,15.000109,17.99916,30.973762)
abundance<-c(0.9894,0.999855,0.996205,0.99757,1.0)

natural.isotopic.abundance<-data.frame(elements=elements,mass1=mass1,mass2=mass2,abundance=abundance)
compound<-c(39,74,1,8,1)


isoltope.distribution<-function(X) {
  output.combined<-matrix(c(0,1),nrow=1,ncol=2)
  for(p in 1:length(elements)) {
    n<-X[p]
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
  output.combined[,2]<-output.combined[,2]*100/max(output.combined[,2])
  output.combined<-output.combined[order(output.combined[,1]),]
  return(output.combined)
}
test<-isoltope.distribution(compound)
plot(test)


compound.list<-matrix(0,nrow=0,ncol=5,dimnames = list(NULL,elements))
compound.list<-rbind(compound.list,"PE34:2"=c(39,74,1,8,1))
compound.list<-rbind(compound.list,"PC34:1"=c(42,82,1,8,1))
compound.list<-rbind(compound.list,"PC36:0"=c(44,88,1,8,1))
compound.isotopes<-apply(compound.list,1,isoltope.distribution)

long.form<-matrix(0,nrow=0,ncol=2,dimnames = list(NULL,c("mass","abundance")))
for(i in 1:length(compound.isotopes)) {
  long.form<-rbind(long.form,compound.isotopes[[i]])
}
plot(long.form,type="h")
