# HapDipKinship
This is a package containing algorithms to estimate kinship in haploid-diploid mixed populations. These algorithms are strongly inspired by KING-robust. 

## Installation
To install this package, run
```{r}
devtools::install_github("YenWenWang/HapDipKinship")
```

## Kinship estimation
This package includes algorithms to estimate diploid-diploid, diploid-haploid and haploid-haploid kinships.  
For example,
```{r}
### need to provide some sort of data here
```

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

Using the given pedigree, we can simulate the pedigree with:
```{r}
HapdipPedigreeSim(ancestrygenomatrix, pedigree, ancestry, PairsOfInterest)
```
This will directly spit out the kinship estimates of the pairs of individuals in PairsOfInterest.


## Reference
Wang et al. in prep. Kinship algorithm for haploid-diploid mixed populations.