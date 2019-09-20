# cluster-func
This code can be used to identify earthquake clusters probabilistically within earthquake catalogues as described in [Bayliss, Naylor and Main, 2019](https://doi.org/10.1093/gji/ggz034).

There are three main components:

    1. code to identify nearest-neighbour events in a space-time-magnitude sense as in Zaliapin and Ben-Zion, 2013 (doi:10.1002/jgrb.50178)
    
    2. Code to fit a Weibull mixture model with MCMC (requires R package gtools (https://cran.r-project.org/web/packages/gtools/index.html))
    
    3. Code to construct probabilistic clusters given the nearest neighbour distances and MCMC fits
    
Code to plot probabilistic networks as in Bayliss et al, 2019, is also included, which uses R packages [igraph](https://igraph.org/r/), [intergraph](http://mbojan.github.io/intergraph/), [GGally](https://cran.r-project.org/web/packages/GGally/index.html), [ggnetwork](https://briatte.github.io/ggnetwork/) and [ggplot2](https://ggplot2.tidyverse.org/).    
    
