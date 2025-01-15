# HapDipKinship
This is a package containing algorithms, exKING-robust and KIMGENS, to estimate kinship in haploid-diploid mixed populations. These algorithms are strongly inspired by KING-robust. 

## Installation
To install this package, run
```{r}
devtools::install_github("YenWenWang/HapDipKinship")
```

## Kinship estimation
This package includes algorithms, exKING-robust and KIMGENS, to estimate diploid-diploid, diploid-haploid and haploid-haploid kinships.  
Example will be given below in simulations.

## Simulations
This package also allows to simulate haplodiploidy ancestry (e.g. Hymenoptera). Alternation of generations (e.g. algae) is still under development.

We first simulate a population with three ancestral subpopulations. 
```{r}
ancestrygenomatrix<-PopulationSim(1000,c(0.05,0.15,0.25))
```

Say we want to draw individuals following either two of the Dirichlet distribution:
- Dir(6,2,0.3)
- Dir(2,6,0.3)

We run
```{r}
ancestry<-matrix(c(6,2,0.3,2,6,0.3),nrow=3)
```

Then, we can simulate the pedigree.  
We provide an example `pedigree`:  
<p>
    <img src="fig/pedigree.png" alt="drawing" width="600"/>  
</p>
(drawn with R package `kinship2`)  

Using the given pedigree, we can simulate with:
```{r}
pedigreegeno<-HapdipPedigreeSim(ancestrygenomatrix, pedigree, ancestry)
```

One can subset `pedigreegeno`, but we will use it directly.
```{r}
ploidy<-data.frame(id=pedigree$id,ploidy=ifelse(pedigree$sex=="F",2,1))
kins<-kinship(pedigreegeno,ploidy)
```


## Reference
Wang, Y-W and Ané, C. 2022. KIMGENS: A novel method to estimate kinship in organisms with mixed haploid diploid genetic systems robust to population structure. Bioinformatics 38(11):3044-3050.