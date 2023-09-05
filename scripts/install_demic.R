if ("demic" %in% rownames(installed.packages())) {
    print("Already installed")
} else {
    install.packages("processx", repos = "https://packagemanager.rstudio.com/cran/latest")
    install.packages("fs", repos = "https://packagemanager.rstudio.com/cran/latest")
    install.packages("wrapr", repos = "https://packagemanager.rstudio.com/cran/latest")
    install.packages("vctrs", repos = "https://packagemanager.rstudio.com/cran/latest")
        
    install.packages("devtools", repos = "https://packagemanager.rstudio.com/cran/latest")
    library(devtools)
    devtools::install_github("Ulthran/DEMIC", ref = "4-thorough-logging-and-rewrite-of-complicated-bits")
}

x <- data.frame()
write.table(x, file=snakemake@output[['out']], col.names=FALSE)