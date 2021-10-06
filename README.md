## serac <a href="https://github.com/rosalieb/serac" target="_blank"><img src="vignettes/figures/hex-serac.png" align="right" height="220" width="190" ></a>

an R package for ShortlivEd RAdionuclide Chronology of recent sediment cores.

To report a problem, email me or use the Github "Issues" tool.

#### Citation

Bruel, R., Sabatier, P., 2020. serac: an R package for ShortlivEd RAdionuclide chronology of recent sediment cores. <i>Journal of Environmental Radioactivity</i> <b>225</b>, 106449. https://doi.org/10.1016/j.jenvrad.2020.106449

#### Download

*serac* is not available on CRAN, but can be downloaded directly from this GitHub repository:

```
install.packages("devtools")
devtools::install_github("rosalieb/serac", build_vignettes = TRUE)
library(serac)
```

See the vignette (`vignette("serac")`) for a complete example of the functionalities of _serac_. We included in the package an example dataset for Lake Allos ([Wilhelm et al., 2012](https://www.sciencedirect.com/science/article/pii/S0033589412000294)), that allows you to reproduce the age-depth model for the core ALO09P12.

<img src="vignettes/figures/ALO_example.png" align="center">
Figure 1. Age-depth model for the core ALO09P12 built with the package _serac_ ([Wilhelm et al., 2012](https://www.sciencedirect.com/science/article/pii/S0033589412000294))

