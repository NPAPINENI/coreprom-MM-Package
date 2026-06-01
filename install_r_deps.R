# install_r_deps.R
# Run once after cloning:
# Rscript install_r_deps.R

# ── Python: rpy2 ──────────────────────────────────────────────────────────────
message("Installing rpy2...")

exit_code <- system("python -m pip install rpy2==3.5.16 --only-binary=:all: --quiet")

if (exit_code != 0) {
    message("rpy2==3.5.16 not available as binary, trying conda...")
    exit_code <- system("conda install -c conda-forge rpy2 -y --quiet")
}

if (exit_code != 0) {
    message("WARNING: Could not install rpy2 automatically.")
    message("Please install manually: conda install -c conda-forge rpy2")
}

# ── R packages ────────────────────────────────────────────────────────────────
message("Installing R dependencies...")
lib <- path.expand("~/R/library")
dir.create(lib, recursive=TRUE, showWarnings=FALSE)

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", lib=lib, repos="https://cloud.r-project.org/")

BiocManager::install(
    c("BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb",
      "GenomicRanges", "Biostrings", "DNAshapeR"),
    lib=lib, ask=FALSE, update=FALSE
)

message("All R dependencies installed successfully.")
message("If rpy2 failed above, run: conda install -c conda-forge rpy2")
