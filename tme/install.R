installed <- rownames(installed.packages(.Library))
pkgs <- c(
'ggplot2',
'dplyr',
'tidyr',
'readr',
'stringr',
'Biobase',
'genefilter',
'limma',
'sva',
'purrr'
)
#'biomaRt'
options(repos=c("https://ftp.osuosl.org/pub/cran/", "https://cran.revolutionanalytics.com/"))

# Install BiocManager
if (!any(rownames(installed.packages()) == "BiocManager")){
  install.packages("BiocManager")
}

# Install pkgs
for(ii in pkgs)
{
if (!any(rownames(installed.packages()) == ii)){
	BiocManager::install(ii, update=FALSE)
    }
}

devtools::install_github("Sage-Bionetworks/CMSclassifier")
