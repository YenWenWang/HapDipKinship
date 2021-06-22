#' Kinship
#'
#' @description Estimates the kinships for all combinations of pairs of individuals.
#'
#' @param genotypematrix A matrix encoding genotypes with columns of individuals and rows of SNP sites. (0=homozygotic reference, 1=heterozygotic, 2=homozygotic alternative).
#'
#' @param ploidy A data frame with two columns. First column is id (same as column names of genotype matrix). Second column is ploidy (1 or 2). Guess ploidy if not provided.
#'
#' @param skipRelBased Binary variable for skipping relative-based kinship estimates or not. F if not provided. See 'Details.'
#'
#' @param RelBasedKinshipThreshold A kinship threshold for defining "relatives" in relative-based kinship. 0.1 if not provided. See 'Details.'
#'
#' @return A vector encoding the genotype of a haploid individual (0,1).
#' @export
#'
#' @examples
#' dip<-sample(c(0:2),10,replace=T)
#' meiosis(dip)
kinship<-function(genotypematrix,ploidy,skipRelBased=F,RelBasedKinshipThreshold=0.1){
  ##check if everyone has ploidy coded
  kinshipmatrix<-as.data.frame(t(combn(colnames(genotypematrix),2)))
  tempkin<-t(apply(kinshipmatrix,1,function(x)
      if(all(ploidy$ploidy[ploidy$id%in%x[1:2]]==2)){
        return(dipking(genotypematrix[,x[1]],genotypematrix[,x[2]]))
      }else if(all(ploidy$ploidy[ploidy$id%in%x[1:2]]==1)){
        return(c(kinship=NA,IBS0=NA))
      }else if(ploidy$ploidy[ploidy$id%in%x[2]]==1&ploidy$ploidy[ploidy$id%in%x[1]]==2){
        return(hapdipking(genotypematrix[,x[1]],genotypematrix[,x[2]]))
      }else if(ploidy$ploidy[ploidy$id%in%x[1]]==1&ploidy$ploidy[ploidy$id%in%x[2]]==2){
        return(hapdipking(genotypematrix[,x[2]],genotypematrix[,x[1]]))
      }else{
        return(c(kinship=NA,IBS0=NA))
      }))
  colnames(tempkin)<-c("kinship","IBS0")
  kinshipmatrix<-cbind(kinshipmatrix,tempkin)

  if(!skipRelBased){
    kinshipmatrix$kinshipRelBased<-NA
    kinshipmatrix$IBS0RelBased<-NA

    for(pair in 1:nrow(kinshipmatrix)){
      inds<-c(kinshipmatrix$V1[pair],kinshipmatrix$V2[pair])
      submatrix<-lapply(inds,function(x)kinshipmatrix[kinshipmatrix$V1==x|kinshipmatrix$V2==x,])
      related<-lapply(submatrix,function(x)which(x$kinship>RelBasedKinshipThreshold))
      related<-mapply(function(x,y)unique(c(y$V1[x][y$V1[x]%in%ploidy$id[ploidy$ploidy==2]],
                        y$V2[x][y$V2[x]%in%ploidy$id[ploidy$ploidy==2]])),related,submatrix)
      if(any(sapply(related,length)!=0)){
        if(all(ploidy$ploidy[ploidy$id%in%inds]==1)){
          kinshipmatrix[pair,c("kinshipRelBased","IBS0RelBased")]<-rowMeans(sapply(related,function(x){
            hapking(genotypematrix[,inds[1]],genotypematrix[,inds[2]],genotypematrix[,x])}))
        }else if(all(ploidy$ploidy[ploidy$id%in%inds]==2)){
          kinshipmatrix[pair,c("kinshipRelBased","IBS0RelBased")]<-rowMeans(sapply(related,function(x){
            dipking(genotypematrix[,inds[1]],genotypematrix[,inds[2]],genotypematrix[,x])}))
        }else{
          if(ploidy$ploidy[ploidy$id%in%inds[1]]==1){
            kinshipmatrix[pair,c("kinshipRelBased","IBS0RelBased")]<-rowMeans(sapply(related,function(x){
              hapdipking(genotypematrix[,inds[2]],genotypematrix[,inds[1]],genotypematrix[,x])}))

          }else{
            kinshipmatrix[pair,c("kinshipRelBased","IBS0RelBased")]<-rowMeans(sapply(related,function(x){
              hapdipking(genotypematrix[,inds[1]],genotypematrix[,inds[2]],genotypematrix[,x])}))

          }
        }

      }else{
        kinshipmatrix[pair,c("kinshipRelBased","IBS0RelBased")]<-NA
      }

    }
  }

  colnames(kinshipmatrix)[c(1:2)]<-c("ind1","ind2")
  return(kinshipmatrix)
}


#' @export
dipking<-function(dip1,dip2,dipref=NA){
  if(all(is.na(dipref))){
    set<-dip1%in%c(0,1,2)&dip2%in%c(0,1,2)
    dip1<-dip1[set]
    dip2<-dip2[set]
    N1<-sum(dip1==1)
    N2<-sum(dip2==1)
    AaAa<-sum(dip1==1&dip2==1)
    AAaa<-sum((dip1==0&dip2==2)|(dip1==2&dip2==0))
    return(c(kinship=1/2-(sum(N1,N2)/min(N1,N2))/4+(AaAa-2*AAaa)/2/min(N1,N2),
             IBS0=AAaa/min(N1,N2)))
  }else{
    if(!is.matrix(dipref)){
      dipref<-matrix(dipref)
    }
    results<-apply(dipref,2,function(x){
      set<-dip1%in%c(0,1,2)&dip2%in%c(0,1,2)&x%in%c(0,1,2)
      dip1<-dip1[set]
      dip2<-dip2[set]
      x<-x[set]
      N<-sum(x==1)
      N1<-sum(dip1==1)
      N2<-sum(dip2==1)
      AaAa<-sum(dip1==1&dip2==1)
      AAaa<-sum((dip1==0&dip2==2)|(dip1==2&dip2==0))
      return(c(kinship=1/2-1/4*((4*AAaa-2*AaAa+N1+N2)/N),IBS0=AAaa/N))})
    kinship<-median(results[1,])
    IBS0<-median(results[2,which(abs(results[1,] - kinship) == min(abs(results[1,] - kinship)))])
    return(c(kinship=kinship,IBS0=IBS0))

  }

}

#' @export
hapdipking<-function(dip,hap,dipref=NA){
  if(all(is.na(dipref))){
    set<-dip%in%c(0,1,2)&hap%in%c(0,1)
    dip<-dip[set]
    hap<-hap[set]
    N<-sum(dip==1)
    AAa<-sum((dip==0&hap==1)|(dip==2&hap==0))
    return(c(kinship=1/2-AAa/N,IBS0=AAa/N))
  }else{
    if(!is.matrix(dipref)){
      dipref<-matrix(dipref)
    }
    results<-apply(dipref,2,function(x){
      set<-dip%in%c(0,1,2)&hap%in%c(0,1)&x%in%c(0,1,2)
      dip<-dip[set]
      hap<-hap[set]
      x<-x[set]
      N<-sum(x==1)
      AAa<-sum((dip==0&hap==1)|(dip==2&hap==0))
      return(c(kinship=1/2-AAa/N,IBS0=AAa/N))})
    kinship<-median(results[1,])
    IBS0<-median(results[2,which(abs(results[1,] - kinship) == min(abs(results[1,] - kinship)))])
    return(c(kinship=kinship,IBS0=IBS0))
  }
}

#' @export
hapking<-function(hap1,hap2,dipref){
  if(!is.matrix(dipref)){
    dipref<-matrix(dipref)
  }
  results<-apply(dipref,2,function(x){
    set<-hap1%in%c(0,1)&hap2%in%c(0,1)&x%in%c(0,1,2)
    hap1<-hap1[set]
    hap2<-hap2[set]
    x<-x[set]
    N<-sum(x==1)
    Aa<-sum((hap1==0&hap2==1)|(hap1==1&hap2==0))
    return(c(kinship=1-Aa/N,IBS0=Aa/N))})

  kinship<-median(results[1,])
  IBS0<-median(results[2,which(abs(results[1,] - kinship) == min(abs(results[1,] - kinship)))])
  return(c(kinship=kinship,IBS0=IBS0))
}

