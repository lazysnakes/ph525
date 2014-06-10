# Basic Bioconductor infrastructure for genomics, microarray and NGS

Open Bioconductor with


```r
biocLite()
```

```
## Error: could not find function "biocLite"
```


if this fails, then install with 
source("http://bioconductor.org/biocLite.R")

## IRanges


```r
library(IRanges)
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
```

```r
ir <- IRanges(5, 10)
ir
```

```
## IRanges of length 1
##     start end width
## [1]     5  10     6
```

```r
start(ir)
```

```
## [1] 5
```

```r
end(ir)
```

```
## [1] 10
```

```r
width(ir)
```

```
## [1] 6
```

```r
`?`(IRanges)
```

```
## starting httpd help server ... done
```



```r
ir <- IRanges(start = c(3, 5, 17), end = c(10, 8, 20))
ir
```

```
## IRanges of length 3
##     start end width
## [1]     3  10     8
## [2]     5   8     4
## [3]    17  20     4
```

```r
ir <- IRanges(5, 10)
```



```r
`?`("intra-range-methods")
shift(ir, -2)
```

```
## IRanges of length 1
##     start end width
## [1]     3   8     6
```


Remeber, all of these commands can work on more than one range at once. Here we show the effects of the different methods using a single range:


```r
shift(ir, -2)
```

```
## IRanges of length 1
##     start end width
## [1]     3   8     6
```

```r
narrow(ir, start = 2)
```

```
## IRanges of length 1
##     start end width
## [1]     6  10     5
```

```r
narrow(ir, end = 5)
```

```
## IRanges of length 1
##     start end width
## [1]     5   9     5
```

```r
flank(ir, width = 3, start = TRUE, both = FALSE)
```

```
## IRanges of length 1
##     start end width
## [1]     2   4     3
```

```r
flank(ir, width = 3, start = FALSE, both = FALSE)
```

```
## IRanges of length 1
##     start end width
## [1]    11  13     3
```

```r
flank(ir, width = 3, start = TRUE, both = TRUE)
```

```
## IRanges of length 1
##     start end width
## [1]     2   7     6
```

```r
ir * 2
```

```
## IRanges of length 1
##     start end width
## [1]     6   8     3
```

```r
ir + 2
```

```
## IRanges of length 1
##     start end width
## [1]     3  12    10
```

```r
ir - 2
```

```
## IRanges of length 1
##     start end width
## [1]     7   8     2
```




```r
# set up a plotting window so we can look at range operations
plotir <- function(ir, i) {
    arrows(start(ir) - 0.5, i, end(ir) + 0.5, i, code = 3, angle = 30, lwd = 3)
}
plot(0, 0, xlim = c(0, 15), ylim = c(0, 11), type = "n", xlab = "", ylab = "", 
    xaxt = "n")
axis(1, 0:15)
abline(v = 0:30 + 0.5, col = rgb(0, 0, 0, 0.5))

# plot the original IRange
plotir(ir, 1)

# draw a red shadow for the original IRange
polygon(c(start(ir) - 0.5, start(ir) - 0.5, end(ir) + 0.5, end(ir) + 0.5), c(-1, 
    12, 12, -1), col = rgb(1, 0, 0, 0.2), border = NA)
plotir(shift(ir, -2), 2)
plotir(narrow(ir, start = 2), 3)
plotir(narrow(ir, end = 5), 4)
plotir(flank(ir, width = 3, start = TRUE, both = FALSE), 5)
plotir(flank(ir, width = 3, start = FALSE, both = FALSE), 6)
plotir(flank(ir, width = 3, start = TRUE, both = TRUE), 7)
plotir(ir * 2, 8)
plotir(ir + 2, 9)
plotir(ir - 2, 10)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



```r
`?`("inter-range-methods")
ir <- IRanges(start = c(3, 5, 17), end = c(10, 8, 20))
range(ir)
```

```
## IRanges of length 1
##     start end width
## [1]     3  20    18
```

```r
reduce(ir)
```

```
## IRanges of length 2
##     start end width
## [1]     3  10     8
## [2]    17  20     4
```

```r
gaps(ir)
```

```
## IRanges of length 1
##     start end width
## [1]    11  16     6
```

```r
disjoin(ir)
```

```
## IRanges of length 4
##     start end width
## [1]     3   4     2
## [2]     5   8     4
## [3]     9  10     2
## [4]    17  20     4
```


## GRanges and GRangesList

### GRanges


```r
library(GenomicRanges)
```

```
## Error: there is no package called 'GenomicRanges'
```

```r
gr <- GRanges("chrZ", IRanges(start = c(5, 10), end = c(35, 45)), strand = "+", 
    seqlengths = c(chrZ = 100L))
```

```
## Error: could not find function "GRanges"
```

```r
gr
```

```
## Error: object 'gr' not found
```

```r
shift(gr, 10)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'shift': Error: object 'gr' not found
```

```r
shift(gr, 80)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'shift': Error: object 'gr' not found
```

```r
trim(shift(gr, 80))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'trim': Error in shift(gr, 80) : 
##   error in evaluating the argument 'x' in selecting a method for function 'shift': Error: object 'gr' not found
```

```r
mcols(gr)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'mcols': Error: object 'gr' not found
```

```r
mcols(gr)$value <- c(-1, 4)
```

```
## Error: object 'gr' not found
```

```r
gr
```

```
## Error: object 'gr' not found
```


### GRangesList


```r
gr2 <- GRanges("chrZ", IRanges(11:13, 51:53))
```

```
## Error: could not find function "GRanges"
```

```r
mcols(gr)$value <- NULL
```

```
## Error: object 'gr' not found
```

```r
grl <- GRangesList(gr, gr2)
```

```
## Error: could not find function "GRangesList"
```

```r
grl
```

```
## Error: object 'grl' not found
```

```r
length(grl)
```

```
## Error: object 'grl' not found
```

```r
grl[[1]]
```

```
## Error: object 'grl' not found
```

```r
mcols(grl)$value <- c(5, 7)
```

```
## Error: object 'grl' not found
```

```r
grl
```

```
## Error: object 'grl' not found
```

```r
mcols(grl)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'mcols': Error: object 'grl' not found
```


### findOverlaps and %over%


```r
gr1 <- GRanges("chrZ", IRanges(c(1, 11, 21, 31, 41), width = 5))
```

```
## Error: could not find function "GRanges"
```

```r
gr2 <- GRanges("chrZ", IRanges(c(19, 33), c(38, 35)))
```

```
## Error: could not find function "GRanges"
```

```r
gr1
```

```
## Error: object 'gr1' not found
```

```r
gr2
```

```
## Error: object 'gr2' not found
```

```r
fo <- findOverlaps(gr1, gr2)
```

```
## Error: error in evaluating the argument 'query' in selecting a method for function 'findOverlaps': Error: object 'gr1' not found
```

```r
fo
```

```
## Error: object 'fo' not found
```

```r
queryHits(fo)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'queryHits': Error: object 'fo' not found
```

```r
subjectHits(fo)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'subjectHits': Error: object 'fo' not found
```

```r
gr1 %over% gr2
```

```
## Error: error in evaluating the argument 'query' in selecting a method for function 'overlapsAny': Error: object 'gr1' not found
```

```r
gr1[gr1 %over% gr2]
```

```
## Error: object 'gr1' not found
```


### Rle and Views


```r
r <- Rle(c(1, 1, 1, 0, 0, -2, -2, -2, rep(-1, 20)))
r
```

```
## numeric-Rle of length 28 with 4 runs
##   Lengths:  3  2  3 20
##   Values :  1  0 -2 -1
```

```r
str(r)
```

```
## Formal class 'Rle' [package "IRanges"] with 4 slots
##   ..@ values         : num [1:4] 1 0 -2 -1
##   ..@ lengths        : int [1:4] 3 2 3 20
##   ..@ elementMetadata: NULL
##   ..@ metadata       : list()
```

```r
as.numeric(r)
```

```
##  [1]  1  1  1  0  0 -2 -2 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
## [24] -1 -1 -1 -1 -1
```

```r
Views(r, start = c(4, 2), end = c(7, 6))
```

```
## Views on a 28-length Rle subject
## 
## views:
##     start end width
## [1]     4   7     4 [ 0  0 -2 -2]
## [2]     2   6     5 [ 1  1  0  0 -2]
```




## ExpressionSet and SummarizedExperiment

To install GEOquery run biocLite("GEOquery"). In fact, 
try biocLite for all missing packages


```r
library(Biobase)
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```r
library(GEOquery)
```

```
## Error: there is no package called 'GEOquery'
```

```r
geoq <- getGEO("GSE9514")
```

```
## Error: could not find function "getGEO"
```

```r
names(geoq)
```

```
## Error: object 'geoq' not found
```

```r
e <- geoq[[1]]
```

```
## Error: object 'geoq' not found
```


### ExpressionSet


```r
dim(e)
```

```
## Error: object 'e' not found
```

```r
exprs(e)[1:3, 1:3]
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'exprs': Error: object 'e' not found
```

```r
dim(exprs(e))
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'exprs': Error: object 'e' not found
```

```r

pData(e)[1:3, 1:6]
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'pData': Error: object 'e' not found
```

```r
dim(pData(e))
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'pData': Error: object 'e' not found
```

```r
names(pData(e))
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'pData': Error: object 'e' not found
```

```r
pData(e)$characteristics_ch1
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'pData': Error: object 'e' not found
```

```r

fData(e)[1:3, 1:3]
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'fData': Error: object 'e' not found
```

```r
dim(fData(e))
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'fData': Error: object 'e' not found
```

```r
names(fData(e))
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'fData': Error: object 'e' not found
```

```r
head(fData(e)$"Gene Symbol")
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'head': Error in fData(e) : 
##   error in evaluating the argument 'object' in selecting a method for function 'fData': Error: object 'e' not found
```

```r
head(rownames(e))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'head': Error in rownames(e) : 
##   error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'e' not found
```

```r

experimentData(e)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'experimentData': Error: object 'e' not found
```

```r
annotation(e)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'annotation': Error: object 'e' not found
```


### Summarized Experiment


```r
library(parathyroidSE)
```

```
## Error: there is no package called 'parathyroidSE'
```

```r
data(parathyroidGenesSE)
```

```
## Warning: data set 'parathyroidGenesSE' not found
```

```r
se <- parathyroidGenesSE
```

```
## Error: object 'parathyroidGenesSE' not found
```

```r
se
```

```
## Error: object 'se' not found
```




```r
dim(se)
```

```
## Error: object 'se' not found
```

```r
assay(se)[1:3, 1:3]
```

```
## Error: could not find function "assay"
```

```r
dim(assay(se))
```

```
## Error: could not find function "assay"
```

```r

colData(se)[1:3, 1:6]
```

```
## Error: could not find function "colData"
```

```r
dim(colData(se))
```

```
## Error: could not find function "colData"
```

```r
names(colData(se))
```

```
## Error: could not find function "colData"
```

```r
colData(se)$treatment
```

```
## Error: could not find function "colData"
```

```r

rowData(se)[1]
```

```
## Error: could not find function "rowData"
```

```r
class(rowData(se))
```

```
## Error: could not find function "rowData"
```

```r
length(rowData(se))
```

```
## Error: could not find function "rowData"
```

```r
head(rownames(se))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'head': Error in rownames(se) : 
##   error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'se' not found
```

```r
metadata(rowData(se))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'metadata': Error: could not find function "rowData"
```

```r

exptData(se)$MIAME
```

```
## Error: could not find function "exptData"
```

```r
abstract(exptData(se)$MIAME)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'abstract': Error: could not find function "exptData"
```

