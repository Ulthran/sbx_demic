if ("demic" %in% rownames(installed.packages())) {
    print("Already installed")
} else {
    install.packages("processx", repos = "https://packagemanager.rstudio.com/cran/latest")
    install.packages("fs", repos = "https://packagemanager.rstudio.com/cran/latest")
    install.packages("devtools", repos = "https://packagemanager.rstudio.com/cran/latest", clean=TRUE)
    library(devtools)

    devtools::install_github("Ulthran/DEMIC")
}

x <- data.frame()
write.table(x, file=snakemake@output[['out']], col.names=FALSE)