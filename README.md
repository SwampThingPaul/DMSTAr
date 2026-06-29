
# DMSTAr <img src="man/figures/DMSTAr_hexish.png" align="right" alt="" width="120" />

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/1157489208.svg)](https://doi.org/10.5281/zenodo.21037938)

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

DMSTAr is a deterministic hydrology and phosphorous mass-balance engine
replicating the VBA code semantics of the original Dynamic Model for
Stormwater Treatment Areas (DMSTA). For complete details regarding DMSTA
background and development see [DMSTA
webpage](http://wwwalker.net/dmsta/) and relevant publications below.

- Kadlec RH, Wallace SD (2009) Treatment wetlands, 2nd edn. CRC Press,
  Boca Raton, FL

- Walker WW, Kadlec RH (2011) Modeling Phosphorus Dynamics in Everglades
  Wetlands and Stormwater Treatment Areas. Critical Reviews in
  Environmental Science and Technology 41:430–446. doi:
  10.1080/10643389.2010.531225

This package was developed from information presented in `dmsta2e.xlsm`
version of the DMSTA model.

Disclaimer from the original DMSTA file -

**DISCLAIMER**

DMSTA is a modeling tool with a constrained range of applicability. It
has been developed and calibrated to information specific to South
Florida. It is intended for use in evaluating Everglades Protection
Project by individuals with experience in hydrologic & water quality
modeling. It should not be exercised in any situation without careful
examination of all features, assumptions and calibrations, as they
relate to a given application and to the supporting research upon which
the calibrations are based. When properly calibrated by the user, the
hydraulics portion of DMSTA is thought to generate predictions that are
adequate for the purpose of simulating phosphorus dynamics. The
hydraulic simulations should not be relied upon for designing flood
control measures, designing levees, for any other purposes in which life
and/or property may be at risk. The user assumes all risks associated
with using the model for designing treatment areas or any other purpose.

## Citing package

``` r
citation('DMSTAr')
#> To cite DMSTAr in publications, please use:
#> 
#>   Julian P (2026). _DMSTAr: An R implementation of the Dynamic Model
#>   for Stormwater Treatment Areas_. doi:10.5281/zenodo.21037939
#>   <https://doi.org/10.5281/zenodo.21037939>, R package version 0.1.1,
#>   <https://github.com/SwampThingPaul/DMSTAr>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {DMSTAr: An R implementation of the Dynamic Model for Stormwater Treatment Areas},
#>     author = {Paul Julian},
#>     year = {2026},
#>     note = {R package version 0.1.1},
#>     doi = {10.5281/zenodo.21037939},
#>     url = {https://github.com/SwampThingPaul/DMSTAr},
#>   }
#> 
#> Please cite the version-specific Zenodo DOI for the version of DMSTAr
#> used in your analysis to support reproducibility.
```

## Installation <a name="install"></a>

Development version can be installed from this repo using:

``` r
install.packages("devtools");# if you do not have it installed on your PC
devtools::install_github("SwampThingPaul/DMSTAr")
```

<!-- To install the deployed version from CRAN use: -->

## DMSTA version and test cases <a name="versions"></a>

Currently `DMSTAr` has been evaluated against published DMSTA results
from various planning efforts.

- Restoration Strategies
  - DMSTA Model Version 2c (Version Date: 7/29/2011)

  - Link to results and documentation - [Restoration Strategies -
    Negotiation
    Phase](https://smmsviewer.apps.sfwmd.gov/modelreport/1687)
- Central Everglades Planning Project, Post Authorization Change Report
  - DMSTA Model Version 2c2b (2010?),
    [CEPP](https://smmsviewer.apps.sfwmd.gov/modelreport/2327) used
    Version 2c2b

  - Link to results and documentation - [DMSTA_C240TSP -
    EAA](https://smmsviewer.apps.sfwmd.gov/modelreport/3426)
- Western Everglades Restoration Plan
  - DMSTA Model Version 2 (incomplete information, assumes version 2c)

  - Link to results and documentation - [Western Everglades Restoration
    Project (WERP)](https://smmsviewer.apps.sfwmd.gov/modelreport/3458)
- Lake Okeechobee System Operating Manual (LOSOM)
  - DMSTA Model Version 2c2b 2010v2b

  - Link to results and documentation - [Lake Okeechobee System
    Operating Manual
    (LOSOM)](https://smmsviewer.apps.sfwmd.gov/modelreport/3580)

<!-- ## *Indevelopment* -->

<!--  * Develop flume and alum based models  -->
