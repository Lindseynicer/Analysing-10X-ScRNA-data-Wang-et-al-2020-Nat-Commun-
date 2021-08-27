Untitled
================

<style type="text/css">
.scroll-100 {
  max-height: 100px;
  overflow-y: auto;
  background-color: inherit;
}
</style>

``` r
#setwd("R_HelenHeZhu_StemCell")

###############################################################################
### Step01 removal of poor quality cells and contaminated non-epithelial cells 
###############################################################################

library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
#

data <- Read10X(data.dir = "GSE111429_RAW")
PG_all    <- CreateSeuratObject(count = data, min.cells = 0, min.genes = 0, project = "PG")
PG_all <- PercentageFeatureSet(PG_all, pattern = "^mt-", col.name = "percent.mt")
before <- VlnPlot(PG_all, features =c('nFeature_RNA','nCount_RNA','percent.mt'))

# Poor-quality cells with less than 1000 genes detected, less than 5000 UMIs or more than 5% UMI mapped to mitochondria genes were removed. 
PG_all <- subset(PG_all, subset = nFeature_RNA > 1400 & nCount_RNA > 5000 & percent.mt < 5)
after <- VlnPlot(PG_all, features =c('nFeature_RNA','nCount_RNA','percent.mt'))
ggarrange(before,after, ncol = 2, nrow = 1)
```

![](abc_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
