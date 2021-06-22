#' Meiosis
#'
#' @description Take a diploid individual and generate gamete
#'
#' @param dip A vector encoding the genotype of a diploid individual (0,1,2).
#'
#' @return A vector encoding the genotype of a haploid individual (0,1).
#' @export
#'
#' @examples
#' dip<-sample(c(0:2),10,replace=T)
#' meiosis(dip)
meiosis<-function(dip){
  out<-ifelse(dip==0,0,ifelse(dip==2,1,-1))
  out[out==-1]<-sample(c(0,1),sum(out==-1),replace = T)
  return(out)
}

#' Anastomosis
#'
#' @description Combine two gametes into a diploid individual
#'
#' @param hap1 A vector encoding the genotype of a haploid individual (0,1).
#'
#' @param hap2 A vector encoding the genotype of a haploid individual (0,1).
#'
#' @return A vector encoding the genotype of a diploid individual (0,1,2).
#' @export
#'
#' @examples
#' hap1<-sample(c(0:1),10,replace=T)
#' hap2<-sample(c(0:1),10,replace=T)
#' anastomosis(hap1,hap2,replace=T)
anastomosis<-function(hap1,hap2){
  return(hap1+hap2)
}

#' @export
PopulationSim<-function(sites,sampsize,thetak,pthresh=0.1,beta=1,popsizethresh=30){

  npop<-length(thetak)
  popsize<-c(round(gtools::rdirichlet(1,rep(1,npop))*sampsize))
  while(any(is.na(popsize))|any(popsize<popsizethresh)){   #redo if some pop too small
    popsize<-c(round(gtools::rdirichlet(1,rep(1,npop))*sampsize))
  }
  x<-0
  while(sum(popsize)<sampsize){   #add individuals if the total number<sampsize
    popsize[x%%npop+1]<-popsize[x%%npop+1]+1
    x<-1
  }
  ps<-rbeta(sites,shape1 = beta,shape2=beta)
  if(!is.na(pthresh)){
    rmsites<-ps<pthresh|ps>1-pthresh
    while(any(rmsites)){
      ps[rmsites]<-rbeta(sum(rmsites),shape1 = beta,shape2=beta)
      rmsites<-ps<pthresh|ps>1-pthresh
    }
  }
  subpopps<-sapply(thetak,function(z)t(sapply(ps,function(x){rbeta(1,x*(1-z)/z,(1-x)*(1-z)/z)})))
  return(subpopps)
}


#' @export
HapdipPedigreeSim<-function(popgenomatrix,pedigree,RelOfInterest,ancestry,ancestryprop=NA,RelBasedKinshipThreshold=0.1)
{
  #ancestry is a matrix with dirichlet params for popgenotypes (nrow(ancestry)==ncol(popgenomatrix))
  #ancestryprop is a vector describing the proportions of each populations in pedigree
  if(any((is.na(pedigree$mother)&!is.na(pedigree$father))|(!is.na(pedigree$mother)&is.na(pedigree$father))))
    stop("can't have single parent that's unknown")
  if(is.na(ancestryprop))
    ancestryprop<-rep(1/ncol(ancestry),ncol(ancestry))
  pedigreegeno<-matrix(-1,ncol=nrow(pedigree),nrow=nrow(popgenomatrix))
  colnames(pedigreegeno)<-pedigree$id

  ##unrelated samples
  unrelated<-pedigree$id[is.na(pedigree$mother)&is.na(pedigree$father)]
  ##get ancestry for unrelated samples
  ancindv<-t(apply(ancestry[,sample(c(1:ncol(ancestry)),length(unrelated),replace = T,prob=ancestryprop)],2,
        function(x)gtools::rdirichlet(1,x)))
  unrelatedprob<-t(apply(popgenomatrix,1,function(y)rowSums(ancindv*y)))
  pedigreegeno[,unrelated[unrelated%in%pedigree$id[pedigree$sex=="F"]]]<-
    apply(unrelatedprob[,1:sum(unrelated%in%pedigree$id[pedigree$sex=="F"])],2,
          function(x)sapply(x,function(y)rbinom(1,2,y)))
  pedigreegeno[,unrelated[unrelated%in%pedigree$id[pedigree$sex=="M"]]]<-
    apply(unrelatedprob[,(sum(unrelated%in%pedigree$id[pedigree$sex=="F"])+1):length(unrelated)],2,
          function(x)sapply(x,function(y)rbinom(1,1,y)))

  ##loop for mating and meiosis
  while(any(pedigreegeno[1,]==-1)){
    temp<-pedigreegeno
    unfinished<-names(which(temp[1,pedigree$id]==-1))
    pedigreesubset<-pedigree[pedigree$id%in%unfinished,]
    workable<-pedigreesubset$id[pedigreegeno[1,pedigreesubset$mother]!=-1&
                                  pedigreegeno[1,pedigreesubset$father]!=-1]
    for(sample in workable){
      if(pedigree$sex[pedigree$id==sample]=="F"){
        temp[,sample]<-
          anastomosis(meiosis(pedigreegeno[,pedigree$mother[pedigree$id==sample]]),
                      pedigreegeno[,pedigree$father[pedigree$id==sample]])
      }else{
        temp[,sample]<-
          meiosis(pedigreegeno[,pedigree$mother[pedigree$id==sample]])
      }
    }
    if(all(pedigreegeno==temp)){
      stop("stuck in loop, check pedigree!")
    }else{
      pedigreegeno<-temp
    }
  }
  ploidy<-data.frame(id=pedigree$id,ploidy=ifelse(pedigree$sex=="F",2,1))
  kins<-kinship(pedigreegeno,ploidy,RelBasedKinshipThreshold=RelBasedKinshipThreshold)
  ## dip pairs
  RelOfInterest$temptempidid<-apply(RelOfInterest,1,function(x)paste(sort(x[c("ind1","ind2")]),collapse = "_"))
  kins$temptempidid<-apply(kins,1,function(x)paste(sort(x[c("ind1","ind2")]),collapse = "_"))
  kins<-kins[,-c(1:2)]
  RelOfInterest<-merge(RelOfInterest,kins,by="temptempidid")[,-1]
  return(RelOfInterest)
}
