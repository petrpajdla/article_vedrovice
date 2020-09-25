# Vedrovice

`Makefile`: Creates the figures from the individual scripts and the text.

Scripts

i. `knit_text.R`: create text of the article  
ii. `source_scr.R`: source all the scripts creating temporary outputs, figures etc.  

1. Data preparation  
    1.1 `dataset.R`: data preparation  
        *output:* `vedrovice_dataset.RDS`  
    3.1 `spat_dataprep.R`: spatial data preparation   
        *output:* `layout.shp`, `window.shp`  
    <!-- 1.2 `local.R`: assign local vs non-local Sr range  
        *output:* `local.RDS`-->  

2. Randomization    
    2.1 `randomization.R`: randomization (simulation) of artefact co-occurences  
        (includes parallel computations on randomly generated matrices)  
        *output:* `cooc_random_mat.RDS`, `v_statistics.csv`, `v_normalized.csv`, 
        `variable_clusters.csv`, `non_random_vars.txt`     
    <!-- 2.2 `exceptionality_index.R`  
        *output:* `exceptionality.csv`   
    2.3 `exceptionality_explore.R`: explore exceptionality   
        *output:* `exc_gg_prop.csv`-->  
    
3. Spatial analysis  
    
<!-- 3.2 `spat_kde.R`: KDE for different sexed bodies-->  
    3.3 `spat_neighbors.R`: similarity of neighboring burials based on Gabriel graph and Delaunay triangulation  
    3.4 `spat_functions.R`: point pattern analysis  
    3.5 `spat_perc.R`: percolation analysis  
    3.6 `spat_perc_clusters.R`: exploring cluster structure at different percolation thresholds   
<!-- 3.4 `spatial_bufer.R`: Similarity of neighboring burials based on couts in buffer zones -->

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

