# Multiple Testing Tutorial (Differential Abundance)

## Background on the data
In this example I use the publicly available data from a study on colorectal cancer:

**Genomic analysis identifies association of Fusobacterium with colorectal carcinoma.** Kostic, A. D., Gevers, D., Pedamallu, C. S., Michaud, M., Duke, F., Earl, A. M., et al. (2012). **Genome research**, 22(2), 292-298.

Data source, from methods section in article:

The 16S gene data set consists of 454 FLX Titanium sequences spanning the V3 to V5 variable regions obtained for 190 samples (95 pairs). Detailed protocols used for 16S amplification and sequencing are available on the HMP Data Analysis and Coordination Center website.

#### Study Abstract:

The tumor microenvironment of colorectal carcinoma is a complex community of genomically altered cancer cells, nonneoplastic cells, and a diverse collection of microorganisms. Each of these components may contribute to carcino genesis; however, the role of the microbiota is the least well understood. We have characterized the composition of the microbiota in colorectal carcinoma using whole genome sequences from nine tumor/normal pairs. Fusobacterium sequences were enriched in carcinomas, confirmed by quantitative PCR and 16S rDNA sequence analysis of 95 carcinoma/normal DNA pairs, while the Bacteroidetes and Firmicutes phyla were depleted in tumors. Fusobacteria were also visualized within colorectal tumors using FISH. These findings reveal alterations in the colorectal cancer microbiota; however, the precise role of Fusobacteria in colorectal carcinoma pathogenesis requires further investigation.

##  Import data with phyloseq

Start by loading phyloseq.

``` r
library("magrittr")
library("phyloseq")
library("data.table")
library("ggplot2")
library("ggrepel")
library("plotly")
library("DESeq2")
```
Import the data that I already processed and organized for you ahead of time.

``` r
kostic = kosticUnmod = readRDS("Kostic2012StageII.RDS")
```
Cleanup OTU names

```r
taxa_names(kostic) <- taxa_names(kostic) %>% paste0("OTU_", .)
```
Cleanup sample names

```r
sample_names(kostic) <- sample_names(kostic) %>% paste0("sa_", .)
```

Remove OTUs that only appeared in those samples

```r
kostic <- filter_taxa(kostic, function(x){sum(x) > 0}, prune = TRUE)
ntaxa(kostic)
```

Require a prevalence threshold

```r
kosticPrevFilt <- filter_taxa(kostic, function(x){sum(x > 0) > 2}, prune = TRUE)
ntaxa(kosticPrevFilt)
```

## DESeq2

In this example I’m using the major sample covariate, DIAGNOSIS, as the study design factor. The focus of this study was to compare the microbiomes of pairs of healthy and cancerous tissues, so this makes sense. Your study could have a more complex or nested design, and you should think carefully about the study design formula, because this is critical to the test results and their meaning. You might even need to define a new factor if none of the variables in your current table appropriately represent your study’s design.

First load DESeq2.

```r
library("DESeq2")
```
Convert phyloseq object to DESeq2 object for the design we want to use.

```r
countData = (kosticPrevFilt %>% otu_table() %>% as("matrix") %>% t())
rownames(countData) <- NULL
colnames(countData) <- NULL
dds <- 
  DESeqDataSetFromMatrix(
    countData = countData,
    colData = (kosticPrevFilt %>% get_variable("DIAGNOSIS") %>% 
                 DataFrame(DIAGNOSIS = .)),
    design = ~DIAGNOSIS)
```
Define a function for computing geometric mean

```r
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

geoMeans = apply(counts(dds), 1, gm_mean, na.rm = FALSE)
ddsDiagnosis = estimateSizeFactors(dds, geoMeans=geoMeans)
```
Run DESeq2 (just one function!).

```r
ddsDiagnosis <- DESeq(ddsDiagnosis)
```

The following `results` function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the `ddsDiagnosis` object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display.

Prepare to join results data.table and taxonomy table.

```r
resdt <-
  ddsDiagnosis %>% 
  results(cooksCutoff = FALSE) %>% 
  as("data.frame") %>% 
  data.table(keep.rownames = "OTU")
resdt$OTU <- taxa_names(kosticPrevFilt)

taxdt <-
  kosticPrevFilt %>% 
  tax_table() %>% as("matrix") %>% data.frame() %>% 
  data.table(keep.rownames = "OTU")
```

Join results data.table and taxonomy table

```r
resdt <- taxdt[resdt, on = "OTU"]
setorder(resdt, pvalue)
resdt %>% head()
```

```r
alpha = 0.05
resdt[, Significant := padj < alpha]
resdt[!is.na(Significant)]
```

```r
resdt[(Significant)]
```

Let’s look at the OTUs that were significantly different between the two tissues. The following makes a nice ggplot2 summary of the results.

```r
pDESeq2 <-
  ggplot(
    data = resdt[!is.na(Significant)][(pvalue < 0.9)],
    mapping = aes(x = log2FoldChange,
                  y = -log10(pvalue),
                  color = Phylum,
                  label = OTU, label1 = Genus)) + 
  geom_point() + 
  geom_point(data = resdt[(Significant)], size = 6) + 
  geom_text_repel(
    size = 3,
    data = resdt[(Significant)], hjust = -0.1,
    mapping = aes(label = paste("Genus:", Genus)),
    color = "black") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  geom_hline(yintercept = -log10(alpha)) + 
  ggtitle("Volcano Plot: DESeq2. Kostic unpaired test.")
pDESeq2
```


