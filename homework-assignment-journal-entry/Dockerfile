# Using the course base Docker image
FROM risserlin/bcb420-base-image:winter2025-arm64
# arm64 for m-series chips

# Installing the additional R packages
RUN R -e "install.packages(c('DESeq2', 'pheatmap', 'enrichplot'), repos='http://cran.rstudio.com/')"

# To only install Bioconductor packages if necessary
# RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
#          BiocManager::install(c('DESeq2', 'enrichplot'))"

# Set the default command to run RStudio Server
CMD ["/init"]
