# Vedrovice

`Makefile`: Creates the figures from the individual scripts and the text.

Scripts

i. `knit_text.R`: create text of the article  
ii. `source_scr.R`: source all the scripts creating temporary outputs, figures etc.  

1. Data preparation  
    1.1 `dataset.R`: data preparation  
        *output:* `vedrovice_dataset.RDS`  
    1.2 `local.R`: assign local vs non-local Sr range  
        *output:* `local.RDS`  

2. Randomization and exceptionality index    
    2.1 `randomization.R`: randomization (simulation) of artefact co-occurences  
        (includes parallel computations on randomly generated matrices)  
        *output:* `cooc_random_mat.RDS`, `v_statistics.csv`, `v_normalized.csv`, 
        `variable_clusters.csv`, `non_random_vars.txt`     
    2.2 `exceptionality_index.R`  
        *output:* `exceptionality.csv`   
    2.3 `exceptionality_explore.R`: explore exceptionality   
        *output:* `exc_gg_prop.csv`  
    
3. Spatial analysis  
    3.1 `spat_dataprep.R`: point pattern analysis  
        *output:* `layout.shp`, `window.shp`  
    3.2 `spatial_kde.R`: KDE for different sexed bodies  
    3.3 `spatial_neighbors.R`: similarity of neighboring burials based on Gabriel graph  
    3.4 `spatial_functions.R`: point pattern analysis  
    3.5 `spatial_percolation.R`: percolation analysis
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

