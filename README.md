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
- Now you will see experimental results listed in matrix. In the left side you can choose any filter you want ro achieve your own datasets.here we choose GRCH38 in genome assembly as the filter.
- Now it shows 8795 results, there are three buttons under the "8795 Results", click the first button and the results will be listed as one in each line
- Click the "Download" button
- you will now get a file.txt which contains a list of URLs to a file containing all experimental metedata and links to download the data.
- choose what you want.
- Inreality, you can search for the specific experiment/assay/biosample you want on search portal only if you know their IDs or some keys, and look at that experiment.Those file data will be listed inside and download button is just beside their ids.
- Open the "file.txt", choose URL: https://www.encodeproject.org/files/ENCFF768GAH/@@download/ENCFF768GAH.tsv, (4.8MB) , download it and save it into a sister directory of your working directory which is called data. (It should be reachable with file.path("..", "data")).
- ENCFF768GAH file contains differential quantifications of gene expression in shRNA RNA-seq of HepG2 experiments.


# 4 Mapping ENSEMBL IDs to HGNC symbols
This assay uses ensembl ID as gene ID. To associate with HGNC symbol, we can map ENSG ID to HGNC symbols by biomart. However,there are some different ENSG IDs that can be mapped to the same HGNC symbol, because same region may have several versions of genome. We need to be careful in the mapping and notice that not all ENSG IDs can be mapped uniquely.

Preparations: packages, functions, files
- make sure the required packages are installed:

readr provides functions to read data which are particularly suitable for large datasets. They are much faster than the built-in read.csv() etc. But caution: these functions return "tibbles", not data frames. (Know the difference.)
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

