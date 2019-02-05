# BCB420.2019.encode
# 1 About the package
This package describes the workflow of how to download file data from Encode, how to map IDs to HGNC IDs and statistics of data matrix.


The package serves dual duty, as an RStudio project, as well as an R package that can be installed. Package checks pass without errors, warnings, or notes.

# 2 Encode database
The Encyclopedia of DNA Elements (ENCODE) is a public research project whose purpose is to identify functional elements in the human genome. Right now it has extended to other organisms like mouse and fly.
It analyzes data in an integrative fashion. Registry of candidate regulatory elements is one of the most important annotations in the integrative level.
One achievement from this project shows that some "useless" elements in genomes are not as useless as we thought before.
we can look at the enhancer region part and use Yufei's GTRD data to annotate the specific TFs that bind there, and identify co-regulated genes.

## 2.1 Data semantics

The Metedata organization contains six accessions to be reused in experimental protocol and computational analysis
- An assay: the replicates will be performed using the same method, on the same kind of biosample, and investigating the same target. Assays may contain one or more biological replicates.
- A biosample: An accessioned biosample refers to a tube or sample of that biological material that is being used.
- A strain or donor: Every strain (for model organisms) and donor (for humans) is given a donor accession. This accession allows multiple samples obtained from a single donor to be grouped together.
- An antibody lot:Each antibody lot(unique) is also associated with characterizations for its target in a species.
- A library: A unique library is accessioned to ensure the correct files are associated with the nucleic acid material that has been created from the biosample.
- A file: Each data file is accessioned. This accession is used as the file name, along with its file format as an extension.

# 3 Data Download and clean up
There are lots of ways to download source data from Encode and you can even access those data on other database like USUC or GO,etc.Here I introduced one method.

- Navigate to  [**Encode**](https://www.encodeproject.org/) , click "HUMAN" and then click filtered Data Matrix
- Now you will see experimental results listed in tables. In the left side you can choose any filter you want to achieve your own datasets.here we choose GRCH38 in genome assembly as the filter.
- Now it shows 8795 results, there are three buttons under the "8795 Results", click the first button and the results will be listed as one in each line
- Click the "Download" button
- you will now get a file.txt which contains a list of URLs to a file containing all experimental metedata and links to download the data.
- choose what you want.
- In reality, you can search for the specific experiment/assay/biosample you want on search portal only if you know their IDs or some keys, and look at that experiment.Those file data will be listed inside and download button is just beside their ids.
- Open the "file.txt", choose URL: https://www.encodeproject.org/files/ENCFF768GAH/@@download/ENCFF768GAH.tsv, (4.8MB) , download it and save it into a sister directory of your working directory which is called data. (It should be reachable with file.path("..", "data")).
- ENCFF768GAH file contains differential quantifications of gene expression in shRNA RNA-seq of HepG2 experiments.


# 4 Mapping ENSEMBL IDs to HGNC symbols
This assay uses ensembl ID as gene ID. To associate with HGNC symbol, we can map ENSG ID to HGNC symbols by biomart. However,there are some different ENSG IDs that can be mapped to the same HGNC symbol, because same region may have several versions of genome. We need to be careful in the mapping and notice that not all ENSG IDs can be mapped uniquely.

Preparations: packages, functions, files
- make sure the required packages are installed:

readr provides functions to read data which are particularly suitable for large datasets. 
```
if (! requireNamespace("readr")) {
  install.packages("readr")
}
```
biomaRt biomaRt is a Bioconductor package that implements the RESTful API of biomart, the annotation framwork for model organism genomes at the EBI. It is a Bioconductor package, and as such it needs to be loaded via the BiocManager,  
```
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
```
Next we source a utility function that we will use later, for mapping ENSP IDs to gene symbols if the mapping can not be achieved directly.

```
source("inst/scripts/recoverIDs.R")
```

Finally we fetch the HGNC reference data from GitHub. The recoverIDs() function requires the HGNC object to be available in the global namespace. (Nb. This shows how to load .RData files directly from a GitHub repository!)

 
```
myURL <- paste0("https://github.com/hyginn/",
                "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
load(url(myURL))  # loads HGNC data frame
```

## 4.1 which ID to map?
```
# Read the differential data, there are several columns: id , fold change, log fold change and pvalues,etc.

tmp <- read.delim("/Users/yuhanzhang/Downloads/BCB420.2019.encode/data/ENCFF768GAH.tsv", header=TRUE, sep="\t")
head(tmp)
  rowID                 id  baseMean baseMeanA baseMeanB foldChange log2FoldChange         pval         padj
1     1 ENSG00000000003.10 2693.0461 2627.1502 2758.9420  1.0501653     0.07061643 3.285189e-01 6.047569e-01
2     2  ENSG00000000005.5    0.0000    0.0000    0.0000         NA             NA           NA           NA
3     3  ENSG00000000419.8 1814.2725 1488.1768 2140.3682  1.4382486     0.52431311 3.483299e-10 2.346446e-09
4     4  ENSG00000000457.9  239.4773  305.9247  173.0298  0.5655962    -0.82215578 9.737404e-07 4.775592e-06
5     5 ENSG00000000460.12  354.1951  234.6976  473.6925  2.0183103     1.01314798 1.714617e-12 1.408435e-11
6     6  ENSG00000000938.8    0.0000    0.0000    0.0000         NA             NA           NA           NA
# id are Ensembl Gene id, each of them would be mapped to hgnc symbol.
# fold change, log2foldchange and pval and adjusted pvalue are the results of difference gene expression.
# some genes does not have values(NA), which need to be deleted 
tmp <- na.omit(tmp)

# delete suffix in ids
tmp$id <-gsub(".[0-9]+$","", tmp$id)

uENSG <- unique(tmp$id) ## 25645 unique ids
```

## 4.2 mapping via biomaRt
```
# Map ENSP to HGNC symbols: open a "Mart" object ..
myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
tmp2 <- biomaRt::getBM(filters = "ensembl_gene_id",
+                       attributes = c("ensembl_gene_id",
+                                      "hgnc_symbol"),
+                       values = uENSG,
+                       mart = myMart)
                                                                                                            
> head(tmp2)
  ensembl_gene_id hgnc_symbol
1 ENSG00000000003      TSPAN6
2 ENSG00000000419        DPM1
3 ENSG00000000457       SCYL3
4 ENSG00000000460    C1orf112
5 ENSG00000000971         CFH
6 ENSG00000001036       FUCA2

# defining mapping tool
ensp2sym <- tmp2$hgnc_symbol
names(ensp2sym) <- tmp2$ensembl_gene_id
> head(ensp2sym)
ENSG00000000003 ENSG00000000419 ENSG00000000457 ENSG00000000460 ENSG00000000971 ENSG00000001036 
       "TSPAN6"          "DPM1"         "SCYL3"      "C1orf112"           "CFH"         "FUCA2" 

```

There are three possible problems we may encounter
- There might be more than one value returned. The ID appears more than once in tmp2$ensembl_gene_id, with different mapped symbols.
```
sum(duplicated(tmp2$ensembl_gene_id)) ## 0 Luckily we don't have this kind of problem.
```
- There might be nothing returned for one ENSG ID. We have the ID in uENSG, but it does not appear in tmp2$ensembl_gene_id:
```
sum(! (uENSG) %in% tmp2$ensembl_gene_id)
[1] 1125
```
- There might be no value returned: NA, or "". The ID appears in `tmp2$ensembl_gene_id`, but there is no symbol in `tmp2$hgnc_symbol`.
```
sum(is.na(ensp2sym)) # 0, we don't need to consider this part.
```

Now we fix these problems.
First, we add the symbols that were not returned by biomaRt to the map. They are present in uENSG, but not in ensp2sym:

 
```
  sel <- ! (uENSG %in% names(ensp2sym))
  x <- rep(NA, sum( sel))
  names(x) <- uENSG[ sel ]

  # confirm uniqueness
  any(duplicated(c(names(x), names(ensp2sym))))  # FALSE

  # concatenate the two vectors
  ensp2sym <- c(ensp2sym, x)

  # confirm
  all(uENSG %in% names(ensp2sym))  # TRUE
  
```
Next we set the symbols for which only an empty string was returned to NA:
```
sel <- which(ensp2sym == "") # 4754 elements
  ensp2sym[head(sel)] # before ...
  ensp2sym[sel] <- NA
  ensp2sym[head(sel)] # ... after

  # Do we still have all ENSG IDs accounted for?
  all( uENSG %in% names(ensp2sym))  # TRUE

```
## 4.2.3 Additional symbols
A function for using biomaRt for more detailed mapping is in the file `inst/scripts/recoverIds.R`. We have loaded it previously, and use it on all elements of ensp2sym that are NA.

 
```
  # How many NAs are there in "ensp2sym" column?
  sum(is.na(ensp2sym))   # 5879

  # subset the ENSG IDs
  unmappedENSG <- names(ensp2sym)[is.na(ensp2sym)]

  # use our function recoverIDs() to try and map the unmapped ensp IDs
  # to symboils via other cross-references
  recoveredENSG <- recoverIDs(unmappedENSG)

  # how many did we find
  nrow(recoveredENSG)  # 11. Not much, but it's honest work.

  # add the recovered symbols to ensp2sym
  ensp2sym[recoveredENSP$ensp] <- recoveredENSP$sym

  # validate:
  sum(is.na(ensp2sym))  # 436 - 11 less than 447
```
