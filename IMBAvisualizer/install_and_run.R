#devtools::document(roclets = c('rd', 'collate', 'namespace', 'vignette'))

source("create_variables.R")
devtools::install(quick=TRUE)
#devtools::build_manual()
library(illav)
illav::run_app()
