repo <- "https://cloud.r-project.org/"

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("multtest")

install.packages("jsonlite", repos = repo)
install.packages("Seurat", repos = repo)
install.packages("SoupX", repos = repo)
