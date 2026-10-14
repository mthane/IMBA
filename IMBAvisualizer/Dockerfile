# IMBAvisualizer shiny_app. Reads lrvTrack CSV tables; it does not track.
# rocker/shiny-verse already contains Shiny and the tidyverse. The remaining
# packages are the library() calls in shiny_app/global.R.
FROM rocker/shiny-verse:4.5.3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake \
        gfortran \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libcairo2-dev \
        libmagick++-dev \
        libfftw3-dev \
        libgsl-dev \
        libudunits2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c( \
      "attempt", "config", "DT", "golem", "processx", "shinythemes", \
      "dashboardthemes", "shinyFiles", "shinyWidgets", "shinyjs", \
      "RcppRoll", "shinydashboard", "colourpicker", "webshot", "corrr", \
      "sortable", "caTools", "gtools", "gmodels", "outliers", "spant", \
      "correlation", "zoo", "ggfortify", "plotly", "viridis", "cvms", \
      "gridExtra", "yardstick", "GGally", "fields", "RColorBrewer", "ggtern", \
      "factoextra", "cluster", "mclust", "FactoMineR", "randomForest", "caret", \
      "corrplot", "ggcorrplot", "ggnewscale", "imagine", "tictoc", "EnvStats", \
      "thinkr" \
    ), repos = "https://cloud.r-project.org", Ncpus = parallel::detectCores())'

# Parallel install can skip a package while its dependency is still building.
# Install anything global.R needs that did not land, and fail the image if it is still missing.
RUN Rscript -e 'need <- c( \
      "sortable", "cvms", "factoextra", "FactoMineR", "spant", "ggtern", \
      "imagine", "dashboardthemes", "shinyFiles", "EnvStats", "mclust", \
      "caret", "ggcorrplot", "ggnewscale", "correlation", "outliers", \
      "gmodels", "webshot", "colourpicker", "shinyWidgets", "shinyjs", \
      "RcppRoll", "shinydashboard", "corrr", "caTools", "gtools", "zoo", \
      "ggfortify", "plotly", "viridis", "gridExtra", "yardstick", "GGally", \
      "fields", "RColorBrewer", "cluster", "randomForest", "corrplot", \
      "tictoc", "DT", "shinythemes", "data.table", "tidyverse" \
    ); inst <- rownames(installed.packages()); miss <- need[!need %in% inst]; \
    if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org"); \
    inst <- rownames(installed.packages()); still <- need[!need %in% inst]; \
    if (length(still)) stop(paste("missing:", paste(still, collapse = ", ")))'

WORKDIR /opt/imba/shiny_app
COPY shiny_app/ ./

ENV IMBA_DATA_ROOT=/data/
EXPOSE 3838
CMD ["Rscript", "-e", "shiny::runApp('.', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]
