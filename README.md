# Vedrovice

Author: Petr Pajdla (<petr.pajdla@protonmail.com>)

This repository contains an article:

> Pajdla, Petr (2022): Spatial patterns and grave goods differences at the Cemetery of Vedrovice (Czech Republic): A resampling approach to the Early Neolithic identity markers.

The code is organized in individual scripts because some of the procedures, especially the randomization experiments, take a fair amount of time to process.

### Contents

`Makefile`: Creates the text of the article in a given format.

#### Code 

Code folder contains scripts necessary to reproduce the results.

i. `knit_text.R`: create text of the article  
i. `knit_text_word.R`: create text in word format  

1. Data preparation  
    * `01-dataset.R`: data preparation  
    * `02-spat_dataprep.R`: spatial data preparation   

2. Randomization    
    * `03-randomization.R`: randomization (simulation) of artifact co-occurrences  
        (includes parallel computations on randomly generated matrices)  
    
3. Spatial analysis  
    * `04-spat_neighbors.R`: similarity of neighboring burials based on Gabriel graph and Delaunay triangulation  
    * `05-spat_functions.R`: point pattern analysis  
    * `06-spat_perc.R`: percolation analysis  
    * `07-spat_perc_clusters.R`: exploring cluster structure at different percolation thresholds   

#### Data 

Data folder contains primary data necessary to reproduce the results.

#### Text

Text folder contains `RMarkdown` files producing the text of the article.

#### Packages

Employed packages (in alphabetic order):

`bookdown`, `broom`, 
`concaveman`, `corrplot`, 
`DescTools`, `devtools`, `dplyr`, 
`forcats`, 
`ggforce`, `ggplot2`,  `ggrepel`, `ggsflabel`, `ggspatial`, `gridExtra`, 
`here`, 
`igraph`, 
`knitr`, 
`parallel`, `percopackage`, `purrr`, 
`randomizr`, `readr`, `reshape2`, `rmarkdown`, 
`sf`, `spatstat`, `spdep`, `stats`, `stringr`, 
`tibble`, `tidyr`, `tidyverse`, 
`vegan`

#### License

Text and figures: [CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

Code: [MIT License](https://petr-pajdla.mit-license.org/)

Data: [CC-0](http://creativecommons.org/publicdomain/zero/1.0/) attribution requested in reuse
