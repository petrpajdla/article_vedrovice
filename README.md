# Vedrovice

`Makefile`: Creates the figures from the individual scripts and the text.

Scripts

0. `knit_text.Rmd`: create text of the article

1. Data preparation  
    1.1 `dataset.R`: data preparation  
    1.2 `local.R`: assign local vs non-local Sr range  

2. Randomization  
    2.1 `randomization.R`: randomization (simulation) of artefact co-occurences  
        (includes parallel computations on randomly generated matrices)  
    2.2 `excetionality_index.R`
    
3. Spatial analysis  
    3.1 `spatial_org.R`: Point pattern analysis  
    3.2 `spatial_kde.R`: KDE for different sexed bodies  
    3.3 `spatial_neighbors.R`: Similarity of neighboring burials based on Gabriel graph  
    3.4 `spatial_functions.R`: Point pattern analysis  
    3.5 `spatial _percolation.R`: Percolation analysis
<!--    3.4 `spatial_bufer.R`: Similarity of neighboring burials based on couts in buffer zones    -->

Employed packages:

`here`
`vegan`
`dplyr`
`readr`
`igraph`
`tidyr`
`ggplot2`
`corrplot`
`ggspatial`
`spdep`
`DescTools`
`MASS`
