## Internal Data
## Input data used for Restoration Strategies Modeling Efforts
## from field DMSTA_Apr2012-Final\PROJECT_SFWMD_EC_01MAR2012_NET_EAA_STA1E.xls

series <- read.csv("./data-raw/Input_PROJECT_SFWMD_EC_01MAR2012_NET_EAA_STA1E.csv")
usethis::use_data(series,internal = FALSE, overwrite = TRUE)
