---
layout: page
title: Basic inference for microarray data
---




# Basic inference for microarray data

We have data for two strains of mice which we will refer to as strain 0 and 1. We want to know which genes are differentially expressed.  We extracted RNA from 12 randomely selected mice from each strain. In one experiment we pooled the RNA from all individuals from each strain and then created 4 replicate samples from this pool. 


```r
library(Biobase)
```

```
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```r
library(devtools)
```

```
## WARNING: Rtools is required to build R packages, but is not currently installed.
## 
## Please download and install Rtools 3.1 from http://cran.r-project.org/bin/windows/Rtools/ and then run find_rtools().
## 
## Attaching package: 'devtools'
## 
## The following objects are masked from 'package:utils':
## 
##     ?, help
## 
## The following object is masked from 'package:base':
## 
##     system.file
```

```r
install_github("dagdata", "genomicsclass")
```

```
## Installing github repo dagdata/master from genomicsclass
## Downloading master.zip from https://github.com/genomicsclass/dagdata/archive/master.zip
## Installing package from C:\Users\deeds\AppData\Local\Temp\Rtmpkvj1OR/master.zip
## Installing dagdata
## "C:/PROGRA~1/R/R-31~1.0/bin/x64/R" --vanilla CMD INSTALL  \
##   "C:\Users\deeds\AppData\Local\Temp\Rtmpkvj1OR\devtools87c2d7547a0\dagdata-master"  \
##   --library="I:/Documents/R/win-library/3.1" --install-tests
```

```r
library(dagdata)
data(maPooling)
e <- maPooling
head(pData(e))
```

```
##          a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a14 b2 b3 b5 b6 b8 b9 b10 b11
## a10       0  0  0  0  0  0  0  0   1   0   0   0  0  0  0  0  0  0   0   0
## a10a11    0  0  0  0  0  0  0  0   1   1   0   0  0  0  0  0  0  0   0   0
## a10a11a4  0  0  1  0  0  0  0  0   1   1   0   0  0  0  0  0  0  0   0   0
## a11       0  0  0  0  0  0  0  0   0   1   0   0  0  0  0  0  0  0   0   0
## a12       0  0  0  0  0  0  0  0   0   0   1   0  0  0  0  0  0  0   0   0
## a12a14    0  0  0  0  0  0  0  0   0   0   1   1  0  0  0  0  0  0   0   0
##          b12 b13 b14 b15
## a10        0   0   0   0
## a10a11     0   0   0   0
## a10a11a4   0   0   0   0
## a11        0   0   0   0
## a12        0   0   0   0
## a12a14     0   0   0   0
```

```r

# install_github('rafalib','ririzarr')
library(rafalib)
```

```
## Loading required package: RColorBrewer
```

```r
mypar()
flipt <- function(m) t(m[nrow(m):1, ])
myimage <- function(m, ...) {
    image(flipt(m), xaxt = "n", yaxt = "n", ...)
}

myimage(as.matrix(pData(e)), col = c("white", "black"), xlab = "experiments", 
    ylab = "individuals", main = "phenoData")
```

![plot of chunk unnamed-chunk-1](figure/basic_inference_microarray-unnamed-chunk-11.png) 

```r

individuals <- which(rowSums(pData(e)) == 1)
individuals
```

```
##   a10   a11   a12   a14    a2    a3 a3tr1 a3tr2    a4    a5    a6    a7 
##     1     4     5     8     9    12    13    14    15    17    19    21 
##    a8    a9   b10   b11   b12   b13   b14   b15    b2    b3 b3tr1 b3tr2 
##    22    24    29    32    33    36    38    39    40    43    44    45 
##    b5    b6    b8    b9 
##    46    48    50    52
```

```r

## remove replicates
names(individuals)
```

```
##  [1] "a10"   "a11"   "a12"   "a14"   "a2"    "a3"    "a3tr1" "a3tr2"
##  [9] "a4"    "a5"    "a6"    "a7"    "a8"    "a9"    "b10"   "b11"  
## [17] "b12"   "b13"   "b14"   "b15"   "b2"    "b3"    "b3tr1" "b3tr2"
## [25] "b5"    "b6"    "b8"    "b9"
```

```r
individuals <- individuals[-grep("tr", names(individuals))]

es <- e[, individuals]
myimage(as.matrix(pData(es)), col = c("white", "black"))
```

![plot of chunk unnamed-chunk-1](figure/basic_inference_microarray-unnamed-chunk-12.png) 

```r

es$group <- factor(as.numeric(grepl("b", colnames(es))))
es$group
```

```
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1
## Levels: 0 1
```


## Plots of gene expression across group

Let's look at 2 pre-selected genes for illustration, which are the same genes from the lecture.


```r
i = 11425
j = 11878
mypar(1, 2)
stripchart(split(exprs(es)[i, ], es$group), vertical = TRUE, method = "jitter", 
    col = c(1, 2), main = "Gene 1", xlab = "Group", pch = 15)
stripchart(split(exprs(es)[j, ], es$group), vertical = TRUE, method = "jitter", 
    col = c(1, 2), main = "Gene 2", xlab = "Group", pch = 15)
```

![plot of chunk unnamed-chunk-2](figure/basic_inference_microarray-unnamed-chunk-2.png) 


## Compute a t-test for each gene (row)


```r
# biocLite('genefilter')
library(genefilter)
```

```
## 
## Attaching package: 'genefilter'
## 
## The following object is masked from 'package:base':
## 
##     anyNA
```

```r
tt <- rowttests(exprs(es), es$group)
head(tt)
```

```
##            statistic        dm p.value
## 1367452_at  -1.14745 -0.092545  0.2635
## 1367453_at  -0.40971 -0.026868  0.6860
## 1367454_at  -0.03675 -0.003017  0.9710
## 1367455_at   1.30954  0.101727  0.2039
## 1367456_at   0.11782  0.006900  0.9073
## 1367457_at  -0.54473 -0.038318  0.5914
```

```r
head(tt, 1)
```

```
##            statistic       dm p.value
## 1367452_at    -1.147 -0.09254  0.2635
```

```r

mean(exprs(es)[1, es$group == 0]) - mean(exprs(es)[1, es$group == 1])
```

```
## [1] -0.09254
```

```r

simple.t <- t.test(exprs(es)[1, ] ~ es$group, var.equal = TRUE)
simple.t$p.value
```

```
## [1] 0.2635
```

```r

tt$p.value[i]
```

```
## [1] 0.08987
```

```r
tt$p.value[j]
```

```
## [1] 1.979e-07
```

```r

mypar(1, 1)
with(tt, plot(dm, -log10(p.value), xlab = "difference in means", main = "'Volcano' plot"))
tt[with(tt, identify(dm, -log10(p.value))), ]
```

![plot of chunk unnamed-chunk-3](figure/basic_inference_microarray-unnamed-chunk-3.png) 

```
## [1] statistic dm        p.value  
## <0 rows> (or 0-length row.names)
```


## Compare with non-parametric tests


```r
es2 <- es[, c(1, 2, 3, 13, 14, 15)]
mypar(1, 1)
stripchart(exprs(es2)[1, ] ~ es2$group, vertical = TRUE, method = "jitter", 
    col = c(1, 2), main = "three samples per group", xlab = "Group", ylab = "", 
    pch = 15)
```

![plot of chunk unnamed-chunk-4](figure/basic_inference_microarray-unnamed-chunk-41.png) 

```r
t.test(exprs(es2)[1, ] ~ es2$group)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  exprs(es2)[1, ] by es2$group
## t = -0.7016, df = 3.888, p-value = 0.5227
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.4688  0.2814
## sample estimates:
## mean in group 0 mean in group 1 
##           9.975          10.069
```

```r
wilcox.test(exprs(es2)[1, ] ~ es2$group)
```

```
## 
## 	Wilcoxon rank sum test
## 
## data:  exprs(es2)[1, ] by es2$group
## W = 4, p-value = 1
## alternative hypothesis: true location shift is not equal to 0
```

```r

y <- 1:6
x <- es2$group
stripchart(y ~ x, vertical = TRUE, method = "jitter", col = c(1, 2), main = "three samples per group", 
    xlab = "Group", ylab = "", pch = 15)
```

![plot of chunk unnamed-chunk-4](figure/basic_inference_microarray-unnamed-chunk-42.png) 

```r
t.test(y ~ x)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  y by x
## t = -3.674, df = 4, p-value = 0.02131
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -5.267 -0.733
## sample estimates:
## mean in group 0 mean in group 1 
##               2               5
```

```r
wilcox.test(y ~ x)
```

```
## 
## 	Wilcoxon rank sum test
## 
## data:  y by x
## W = 0, p-value = 0.1
## alternative hypothesis: true location shift is not equal to 0
```

```r

y <- c(1:3, 11:13)
stripchart(y ~ x, vertical = TRUE, method = "jitter", col = c(1, 2), main = "three samples per group", 
    xlab = "Group", ylab = "", pch = 15)
```

![plot of chunk unnamed-chunk-4](figure/basic_inference_microarray-unnamed-chunk-43.png) 

```r
t.test(y ~ x)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  y by x
## t = -12.25, df = 4, p-value = 0.0002552
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -12.267  -7.733
## sample estimates:
## mean in group 0 mean in group 1 
##               2              12
```

```r
wilcox.test(y ~ x)
```

```
## 
## 	Wilcoxon rank sum test
## 
## data:  y by x
## W = 0, p-value = 0.1
## alternative hypothesis: true location shift is not equal to 0
```


## Basic inference on microarray using lmFit from limma package


```r
# biocLite('limma')
library(limma)
```

```
## Error: there is no package called 'limma'
```

```r
# ?lmFit
design <- model.matrix(~es$group)
design
```

```
##    (Intercept) es$group1
## 1            1         0
## 2            1         0
## 3            1         0
## 4            1         0
## 5            1         0
## 6            1         0
## 7            1         0
## 8            1         0
## 9            1         0
## 10           1         0
## 11           1         0
## 12           1         0
## 13           1         1
## 14           1         1
## 15           1         1
## 16           1         1
## 17           1         1
## 18           1         1
## 19           1         1
## 20           1         1
## 21           1         1
## 22           1         1
## 23           1         1
## 24           1         1
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$`es$group`
## [1] "contr.treatment"
```

```r
fit <- lmFit(es, design)
```

```
## Error: could not find function "lmFit"
```

```r
names(fit)
```

```
## Error: object 'fit' not found
```

```r
head(coef(fit))
```

```
## Error: object 'fit' not found
```

```r
tt[1, ]
```

```
##            statistic       dm p.value
## 1367452_at    -1.147 -0.09254  0.2635
```

```r
# we will introduce the eBayes() function in a later module called
# 'hierarchical modeling' but we call it now as it is standard in microarray
# analysis
fit <- eBayes(fit)
```

```
## Error: could not find function "eBayes"
```

```r
names(fit)
```

```
## Error: object 'fit' not found
```

```r
fit$p.value[1, ]
```

```
## Error: object 'fit' not found
```

```r
fit$t[1, ]
```

```
## Error: object 'fit' not found
```

```r
tt[1, ]
```

```
##            statistic       dm p.value
## 1367452_at    -1.147 -0.09254  0.2635
```

```r
plot(-1 * tt$statistic, fit$t[, 2], xlab = "rowttests", ylab = "eBayes t")
```

```
## Error: error in evaluating the argument 'y' in selecting a method for function 'plot': Error: object 'fit' not found
```

```r
abline(0, 1, col = "red", lwd = 3)
```

```
## Error: plot.new has not been called yet
```

```r
head(topTable(fit, coef = 2, sort.by = "p"), 3)
```

```
## Error: could not find function "topTable"
```

