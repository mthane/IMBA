#devtools::document(roclets = c('rd', 'collate', 'namespace', 'vignette'))

source("create_variables.R")
devtools::install(quick=TRUE)
#devtools::build_manual()
library(illav)

run_app()
devtools::install_github("https://github.com/mthane/IMBA",subdir = "IMBAvisualizer")
