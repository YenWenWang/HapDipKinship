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
#' dip<-rbinom(1000,2,0.5)
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
#' hap1<-rbinom(1000,1,0.5)
#' hap2<-rbinom(1000,1,0.5)
#' anastomosis(hap1,hap2,replace=T)
anastomosis<-function(hap1,hap2){
  return(hap1+hap2)
}

#' Simulate Population
#'
#' @description Simulate population with substructure using Balding-Nichols model.
#'
#' @param sites Number of sites simulated.
#'
#' @param thetak An x-length vector of Fsts for x ancestral subpopulations (ancestral subpopulations vs total ancestral population).
#'
#' @param pthresh A threshold to avoid low frequency SNPs. 0.1 if not provided.
#'
#' @param beta
#' Shape parameter for beta distribution for simulating original allele frequency. 1 if not provided.
#'
#' @return A matrix of allele frequencies in ancestral subpopulations.
#'
#' @seealso [HapDipKinship::HapdipPedigreeSim()]
#'
#' @export
#'
#' @examples
#' PopulationSim(1000,c(0.05,0.15,0.25))
PopulationSim<-function(sites,thetak,pthresh=0.1,beta=1){
  npop<-length(thetak)
  ps<-rbeta(sites,shape1 = beta,shape2=beta)
  if(!is.na(pthresh)){
    rmsites<-ps<pthresh|ps>1-pthresh
    while(any(rmsites)){
      ps[rmsites]<-rbeta(sum(rmsites),shape1 = beta,shape2=beta)
      rmsites<-ps<pthresh|ps>1-pthresh
    }
  }
  ancestryps<-sapply(thetak,function(z)t(sapply(ps,function(x){rbeta(1,x*(1-z)/z,(1-x)*(1-z)/z)})))
  return(ancestryps)
}


#' Simulate Pedigree
#'
#' @description Simulate a haplodiploidy population pedigree and calculate kinships from a few subpopulations
#'
#' @param ancestrygenomatrix A matrix of allele frequencies in ancestral subpopulations.
#'
#' @param pedigree A data frame describing the simulation of a pedigree. See [HapDipKinship::pedigree].
#'
#' @param ancestry
#' A matrix including dirichlet params for the admixture of ancestral subpopulations.
#' Rows denote ancestral subpopulations. Columns denote different admixture schemes.
#'
#' @param ancestryprop
#' A vector with the length of number of column of ancestry
#' denoting the probability of drawing individuals from each admixed populations.
#' Average sampling if none give.
#'
#' @return Kinships of PairsOfInterest from simulated populations.
#' @seealso
#' [HapDipKinship::PopulationSim()],
#' [HapDipKinship::kinship()]
#'
#' @export
#'
#' @examples
#' ancestrygenomatrix<-PopulationSim(1000,c(0.05,0.15,0.25))
#' ancestry<-matrix(c(6,2,0.3,2,6,0.3),nrow=3)
#' HapdipPedigreeSim(ancestrygenomatrix,pedigree,ancestry)
HapdipPedigreeSim<-function(ancestrygenomatrix,pedigree,ancestry,ancestryprop=NA)
{
  if(any((is.na(pedigree$mother)&!is.na(pedigree$father))|(!is.na(pedigree$mother)&is.na(pedigree$father))))
    stop("can't have single parent that's unknown")
  if(any(is.na(ancestryprop)))
    ancestryprop<-rep(1/ncol(ancestry),ncol(ancestry))
  pedigreegeno<-matrix(-1,ncol=nrow(pedigree),nrow=nrow(ancestrygenomatrix))
  colnames(pedigreegeno)<-pedigree$id

  ##unrelated samples
  unrelated<-pedigree$id[is.na(pedigree$mother)&is.na(pedigree$father)]
  ##get ancestry for unrelated samples
  ancindv<-t(apply(ancestry[,sample(c(1:ncol(ancestry)),length(unrelated),replace = T,prob=ancestryprop)],2,
        function(x)gtools::rdirichlet(1,x)))
  unrelatedprob<-t(apply(ancestrygenomatrix,1,function(y)rowSums(ancindv*y)))
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
  return(pedigreegeno)
}

#' Extract pairs of interest
#'
#' @description A helper function to extract kinships of specific pairs of individuals.
#'
#' @param kinships
#' A data frame with at least two columns: ind1 and ind2
#' @param PairsOfInterest
#' A data frame describing focal relationships.
#' All ids in ind1 and ind2 should be found in kinships (the order is not important).
#' Report all if none provided. See [HapDipKinship::PairsOfInterest].
#'
#' @export
#' @return Reduced data frame of Kinship.
#' @examples
#' ancestrygenomatrix<-PopulationSim(1000,c(0.05,0.15,0.25))
#' ancestry<-matrix(c(6,2,0.3,2,6,0.3),nrow=3)
#' kins<-HapdipPedigreeSim(ancestrygenomatrix,pedigree,ancestry)
#' ExtractPairsOfInterest(kins,PairsOfInterest)
ExtractPairsOfInterest<-function(kinships,PairsOfInterest){
  PairsOfInterest$temptempidid<-apply(PairsOfInterest,1,function(x)paste(sort(x[c("ind1","ind2")]),collapse = "_"))
  kinships$temptempidid<-apply(kinships,1,function(x)paste(sort(x[c("ind1","ind2")]),collapse = "_"))
  kinships<-kinships[,-c(1:2)]
  PairsOfInterest<-merge(PairsOfInterest,kinships,by="temptempidid")[,-1]
  return(PairsOfInterest)
}


