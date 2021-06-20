#' @export
meiosis<-function(dip){
  out<-ifelse(dip==0,0,ifelse(dip==2,1,-1))
  out[out==-1]<-sample(c(0,1),sum(out==-1),replace = T)
  return(out)
}

#' @export
anastomosis<-function(hap1,hap2){
  return(hap1+hap2)
}

#' @export
admix_prop_completeadmixture<-function (labs, subpops = sort(unique(labs))) 
{
  ## bnpsd extension
  if (!all(labs %in% subpops)) 
    stop("provided `subpops` does not contain all labels in `labs`!")
  n_ind <- length(labs)
  k_subpops <- length(subpops)
  admix_proportions <- matrix(1/length(subpops), nrow = n_ind, ncol = k_subpops)
  colnames(admix_proportions) <- subpops
  admix_proportions
}

#' @export
admix_prop_subadmixture<-function (labs, admixtheme, dirichletparam, subpops = sort(unique(labs))) 
{
  ## bnpsd extension
  if (!all(labs %in% subpops)) 
    stop("provided `subpops` does not contain all labels in `labs`!")
  if (length(admixtheme)!=length(labs))
    stop("`admixtheme` does not contain exact same individuals in `labs`!")
  if(!is.matrix(dirichletparam))
    stop("`dirichletparam` is not a matrix of dirichlet parameters for admixtheme!")
  if (!all(unique(admixtheme)%in%colnames(dirichletparam)))
    stop("some `admixtheme` are not included in the columns of `admixtheme`!")
  if (!all(labs%in%row.names(dirichletparam)))
    stop("some labels in `labs` are not included in the rows of `admixtheme`!")
  
  ## dirichlet
  
  n_ind <- length(labs)
  k_subpops <- length(subpops)
  unsortedproportions<-list()
  x<-1
  while(x<length(table(admixtheme))+1){
    unsortedproportions<-c(unsortedproportions,list(gtools::rdirichlet(table(admixtheme)[x],dirichletparam[,x])))
    x<-x+1
  }
  names(unsortedproportions)<-names(table(admixtheme))
  x<-1
  while(x<length(unsortedproportions)+1){
    unsortedproportions[[x]]<-cbind(unsortedproportions[[x]],
                                    paste(names(unsortedproportions)[x],c(1:nrow(unsortedproportions[[x]]))))
    x<-x+1
  }
  unsortedproportions<-do.call(rbind,unsortedproportions)
  
  admixtheme_numbered<-paste(admixtheme,ave(admixtheme, admixtheme, FUN=seq_along))
  admix_proportions<-unsortedproportions[match(unsortedproportions[,ncol(unsortedproportions)],admixtheme_numbered),]
  admix_proportions<-admix_proportions[,-ncol(admix_proportions)]
  admix_proportions<-t(apply(admix_proportions,1,as.numeric))
  colnames(admix_proportions) <- subpops
  admix_proportions
}

#' @export
HapdipPedigreeSim<-function(genotypematrix,pedigree,RelOfInterest,RelBasedKinshipThreshold=0.1)
{
  if(any((is.na(pedigree$mother)&!is.na(pedigree$father))|(!is.na(pedigree$mother)&is.na(pedigree$father))))
    stop("can't have single parent that's unknown")
  pedigreegeno<-matrix(-1,ncol=nrow(pedigree),nrow=nrow(genotypematrix))
  colnames(pedigreegeno)<-pedigree$id
  
  ##unrelated samples
  femaleunrelated<-pedigree$id[is.na(pedigree$mother)&is.na(pedigree$father)&pedigree$sex=="F"]
  maleunrelated<-pedigree$id[is.na(pedigree$mother)&is.na(pedigree$father)&pedigree$sex=="M"]
  s<-sample(ncol(genotypematrix),sum(length(femaleunrelated),length(maleunrelated)))
  pedigreegeno[,femaleunrelated]<-genotypematrix[,s[1:length(femaleunrelated)]]
  pedigreegeno[,maleunrelated]<-apply(
    genotypematrix[,s[(length(femaleunrelated)+1):(length(femaleunrelated)+length(maleunrelated))]],2,meiosis)
  
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