
# DMSTAr <img src="man/figures/DMSTAr_hexish.png" align="right" alt="" width="120" />

DMSTAr is a deterministic hydrology and phosphorous mass-balance engine
replicating the VBA code semantics of the original Dynamic Model for
Stormwater Treatment Areas (DMSTA). For complete details regarding DMSTA
background and development see [DMSTA
webpage](https://wwwalker.net/dmsta/) and relevant publications below.

- Kadlec RH, Wallace SD (2009) Treatment wetlands, 2nd edn. CRC Press,
  Boca Raton, FL

- Walker WW, Kadlec RH (2011) Modeling Phosphorus Dynamics in Everglades
  Wetlands and Stormwater Treatment Areas. Critical Reviews in
  Environmental Science and Technology 41:430–446. doi:
  10.1080/10643389.2010.531225

## Citing package

``` r
citation('DMSTAr')
```

    ## 
    ## To cite package 'DMSTAr' in publications use:
    ## 
    ##   Paul Julian (2026). DMSTAr: Dynamic Model for Stormwater Treatment
    ##   Areas in R. R package version 0.1.0.
    ##   https://github.com/SwampThingPaul/DMSTAr
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {DMSTAr: Dynamic Model for Stormwater Treatment Areas in R},
    ##     author = {Paul Julian},
    ##     year = {2026},
    ##     note = {R package version 0.1.0},
    ##     url = {https://github.com/SwampThingPaul/DMSTAr},
    ##   }

## Installation <a name="install"></a>

Development version can be installed from this repo using:

``` r
install.packages("devtools");# if you do not have it installed on your PC
devtools::install_github("SwampThingPaul/DMSTAr")
```

To install the deployed version from CRAN use:

``` r
install.packages("DMSTAr")
```
