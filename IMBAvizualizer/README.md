
<!-- README.md is generated from README.Rmd. Please edit that file -->

# illav

<!-- badges: start -->

<!-- badges: end -->

The goal of illav is to …

## Installation

You can install the released version of illav from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("illav")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mthane/illav_v0")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(illav)
#> Warning: replacing previous import 'shiny::dataTableOutput' by
#> 'DT::dataTableOutput' when loading 'illav'
#> Warning: replacing previous import 'shiny::renderDataTable' by
#> 'DT::renderDataTable' when loading 'illav'
#> Warning: replacing previous import 'data.table::last' by 'dplyr::last' when
#> loading 'illav'
#> Warning: replacing previous import 'data.table::first' by 'dplyr::first' when
#> loading 'illav'
#> Warning: replacing previous import 'data.table::between' by 'dplyr::between'
#> when loading 'illav'
#> Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when
#> loading 'illav'
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'illav'
#> Warning: replacing previous import 'stats::sd' by 'spant::sd' when loading
#> 'illav'
#> Warning: replacing previous import 'data.table::shift' by 'spant::shift' when
#> loading 'illav'
#> Registered S3 methods overwritten by 'parameters':
#>   method                           from      
#>   as.double.parameters_kurtosis    datawizard
#>   as.double.parameters_skewness    datawizard
#>   as.double.parameters_smoothness  datawizard
#>   as.numeric.parameters_kurtosis   datawizard
#>   as.numeric.parameters_skewness   datawizard
#>   as.numeric.parameters_smoothness datawizard
#>   print.parameters_distribution    datawizard
#>   print.parameters_kurtosis        datawizard
#>   print.parameters_skewness        datawizard
#>   summary.parameters_kurtosis      datawizard
#>   summary.parameters_skewness      datawizard
#> Warning: replacing previous import 'ggplot2::last_plot' by 'plotly::last_plot'
#> when loading 'illav'
#> Warning: replacing previous import 'stats::filter' by 'plotly::filter' when
#> loading 'illav'
#> Warning: replacing previous import 'scales::viridis_pal' by
#> 'viridis::viridis_pal' when loading 'illav'
#> Warning: replacing previous import 'shiny::validate' by 'cvms::validate' when
#> loading 'illav'
#> Warning: replacing previous import 'dplyr::combine' by 'gridExtra::combine' when
#> loading 'illav'
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2
#> Warning: replacing previous import 'ggfortify::unscale' by 'fields::unscale'
#> when loading 'illav'
#> Registered S3 methods overwritten by 'ggtern':
#>   method           from   
#>   grid.draw.ggplot ggplot2
#>   plot.ggplot      ggplot2
#>   print.ggplot     ggplot2
#> Warning: replacing previous import 'ggplot2::theme_dark' by 'ggtern::theme_dark'
#> when loading 'illav'
#> Warning: replacing previous import 'gridExtra::grid.arrange' by
#> 'ggtern::grid.arrange' when loading 'illav'
#> Warning: replacing previous import 'ggplot2::theme_void' by 'ggtern::theme_void'
#> when loading 'illav'
#> Warning: replacing previous import 'ggplot2::annotate' by 'ggtern::annotate'
#> when loading 'illav'
#> Warning: replacing previous import 'ggplot2::ggplot_build' by
#> 'ggtern::ggplot_build' when loading 'illav'
#> Warning: replacing previous import 'ggplot2::ggplot_gtable' by
#> 'ggtern::ggplot_gtable' when loading 'illav'
#> Warning: replacing previous import 'ggplot2::theme_classic' by
#> 'ggtern::theme_classic' when loading 'illav'
#> Warning: replacing previous import 'ggplot2::theme_gray' by 'ggtern::theme_gray'
#> when loading 'illav'
#> Warning: replacing previous import 'ggplot2::aes' by 'ggtern::aes' when loading
#> 'illav'
#> Warning: replacing previous import 'ggplot2::theme_bw' by 'ggtern::theme_bw'
#> when loading 'illav'
#> Warning: replacing previous import 'ggplot2::ggplot' by 'ggtern::ggplot' when
#> loading 'illav'
#> Warning: replacing previous import 'gridExtra::arrangeGrob' by
#> 'ggtern::arrangeGrob' when loading 'illav'
#> Warning: replacing previous import 'ggplot2::ggplotGrob' by 'ggtern::ggplotGrob'
#> when loading 'illav'
#> Warning: replacing previous import 'ggplot2::ggsave' by 'ggtern::ggsave' when
#> loading 'illav'
#> Warning: replacing previous import 'ggplot2::theme_linedraw' by
#> 'ggtern::theme_linedraw' when loading 'illav'
#> Warning: replacing previous import 'ggplot2::theme_minimal' by
#> 'ggtern::theme_minimal' when loading 'illav'
#> Warning: replacing previous import 'ggplot2::theme_light' by
#> 'ggtern::theme_light' when loading 'illav'
#> Warning: replacing previous import 'ggplot2::layer_data' by 'ggtern::layer_data'
#> when loading 'illav'
#> Warning: replacing previous import 'shiny::em' by 'mclust::em' when loading
#> 'illav'
#> Warning: replacing previous import 'outliers::outlier' by
#> 'randomForest::outlier' when loading 'illav'
#> Warning: replacing previous import 'ggplot2::margin' by 'randomForest::margin'
#> when loading 'illav'
#> Warning: replacing previous import 'gridExtra::combine' by
#> 'randomForest::combine' when loading 'illav'
#> Warning: replacing previous import 'yardstick::precision' by 'caret::precision'
#> when loading 'illav'
#> Warning: replacing previous import 'yardstick::recall' by 'caret::recall' when
#> loading 'illav'
#> Warning: replacing previous import 'yardstick::specificity' by
#> 'caret::specificity' when loading 'illav'
#> Warning: replacing previous import 'yardstick::sensitivity' by
#> 'caret::sensitivity' when loading 'illav'
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
