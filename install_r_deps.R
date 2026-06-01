lib <- path.expand("~/R/library")
dir.create(lib, recursive=TRUE, showWarnings=FALSE)

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", lib=lib, repos="https://cloud.r-project.org/")

BiocManager::install(
    c("GenomicRanges", "S4Vectors", "IRanges",
      "GenomeInfoDb", "Biostrings", "DNAshapeR"),
    lib=lib, ask=FALSE, update=FALSE
)

cat("All R dependencies installed successfully.\n")
