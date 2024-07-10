check_package_available <- function(pkg, repo=c("CRAN", "Bioconductor")) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    repo=match.arg(repo)
    installmsg=switch(repo,
                      CRAN=paste0("install.packages('",pkg,"')"),
                      Bioconductor=paste0('if (!require("BiocManager", quietly = TRUE))',
                                          '\n  install.packages("BiocManager")',
                                          '\nBiocManager::install("',pkg,'")'))
    stop("Please install suggested package: ", pkg, " by doing\n",
         installmsg, call. = F)
  }
}
