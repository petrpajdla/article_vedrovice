# Vedrovice

`Makefile`: Creates the figures from the individual scripts and the text.

Scripts

0. `knit_text.Rmd`: create text of the article

1. Data preparation  
    1.1 `dataset.R`: data preparation

2. Randomization  
    2.1 `randomization.R`: randomization (simulation) of artefact co-occurences  
        (includes parallel computations on randomly generated matrices)
    
3. Spatial analysis  
    3.1 `spatial_org.R`: Point pattern analysis  
    3.2 `spatial_neighbors.R`: Similarity of neighboring burials

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
