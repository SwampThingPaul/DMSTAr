## Internal Data
## Current DMSTA calibration values

df <- read.delim("clipboard",header = FALSE)

template <- data.frame(
  Set = character(),
  Descript = character(),
  C0 = numeric(),
  C1 = numeric(),
  C2 = numeric(),
  Ks = numeric(),
  Z1 = numeric(),  Z2 = numeric(),  Z3 = numeric(),
  K1 = numeric(),
  C0_NEWS_SF = numeric(),
  C1_Peri = numeric(),
  Ks_Peri = numeric(),
  Zx_Peri = numeric(),
  Sm = numeric(),  Sb = numeric(),
  MinDepth = numeric(),  MaxDepth = numeric(),
  MinQW = numeric(),MaxQW = numeric(),
  MinConc = numeric(), MaxConc = numeric(),
  MinFreqZ_LT10cm = numeric(),
  MaxFreqZ_LT10cm = numeric(),
  K_CV = numeric(),
  stringsAsFactors = FALSE
)

names(df) <- names(template)

dmsta_cals <- list()
dmsta_cals[[1]] <- data.frame(
  Set = "EMG_3", Descript = "2005 - Emergent ",
  C0 = 3L, C1 = 22L, C2 = 300L, Ks = 16.8, Z1 = 40L, Z2 = 100L,
  Z3 = 200L, K1 = NA, C0_NEWS_SF = NA, C1_Peri = NA, Ks_Peri = NA,
  Zx_Peri = NA, Sm = NA, Sb = NA, MinDepth = 35L, MaxDepth = 76L,
  MinQW = 26L, MaxQW = 210L, MinConc = 19.5, MaxConc = 800L,
  MinFreqZ_LT10cm = 0, MaxFreqZ_LT10cm = 0.09, K_CV = 0.2
)
dmsta_cals[[2]] <- data.frame(
  Set = "PEW_3", Descript = "2005 - Pre-Existent Wetland",
  C0 = 3L, C1 = 22L, C2 = 300L, Ks = 34.9, Z1 = 40L, Z2 = 100L,
  Z3 = 200L, K1 = NA, C0_NEWS_SF = NA, C1_Peri = NA, Ks_Peri = NA,
  Zx_Peri = NA, Sm = NA, Sb = NA, MinDepth = 38L, MaxDepth = 66L,
  MinQW = 69L, MaxQW = 276L, MinConc = 8, MaxConc = 110L, MinFreqZ_LT10cm = 0,
  MaxFreqZ_LT10cm = 0.13, K_CV = 0.21
)
dmsta_cals[[3]] <- data.frame(
  Set = "SAV_3", Descript = "2005 - SAV ", C0 = 3L,
  C1 = 22L, C2 = 300L, Ks = 52.5, Z1 = 40L, Z2 = 100L, Z3 = 200L,
  K1 = NA, C0_NEWS_SF = NA, C1_Peri = NA, Ks_Peri = NA, Zx_Peri = NA,
  Sm = NA, Sb = NA, MinDepth = 62L, MaxDepth = 87L, MinQW = 162L,
  MaxQW = 374L, MinConc = 14.9, MaxConc = 153L, MinFreqZ_LT10cm = 0,
  MaxFreqZ_LT10cm = 0, K_CV = 0.16
)
dmsta_cals[[4]] <- data.frame(
  Set = "PSTA_3", Descript = "2005 - PSTA ", C0 = 3L,
  C1 = 22L, C2 = 300L, Ks = 23.6, Z1 = NA_integer_, Z2 = 100L,
  Z3 = 200L, K1 = NA, C0_NEWS_SF = NA, C1_Peri = NA, Ks_Peri = NA,
  Zx_Peri = NA, Sm = NA, Sb = NA, MinDepth = 13L, MaxDepth = 60L,
  MinQW = 3L, MaxQW = 132L, MinConc = 5.9, MaxConc = 56L, MinFreqZ_LT10cm = 0,
  MaxFreqZ_LT10cm = 0.38, K_CV = 0.22
)
dmsta_cals[[5]] <- data.frame(
  Set = "RES_3", Descript = "2005 - Reservoirs - with depth penalty",
  C0 = 3L, C1 = 150L, C2 = NA_integer_, Ks = 5, Z1 = 40L, Z2 = 100L,
  Z3 = 400L, K1 = NA, C0_NEWS_SF = NA, C1_Peri = NA, Ks_Peri = NA,
  Zx_Peri = NA, Sm = NA, Sb = NA, MinDepth = 90L, MaxDepth = 304L,
  MinQW = 68L, MaxQW = 1135L, MinConc = 50.3, MaxConc = 1144L,
  MinFreqZ_LT10cm = 0, MaxFreqZ_LT10cm = 0, K_CV = 0.45
)

dmsta_cals <- do.call(rbind,dmsta_cals)
usethis::use_data(dmsta_cals,internal = FALSE, overwrite = TRUE)


names(df) <- names(template)
dmsta_cals_obsolete <- data.frame(
  Set = c("EMG_3", "PEW_3", "SAV_3", "PSTA_3", "RES_3",
          "SAV_4", "NEWS", "NEWS_2", "WCA2A", "EWS2", "K123", "RES_2",
          "EMERG", "PSTA", "SAV", "NEWS", "SAV_C4", "RESERV_2", "EMERG",
          "PSTA", "SAV", "NEWS", "SAV_C4", "RES_3L", "SAV_4", "SAV_3L",
          "SAV_5", "NEWS_2", "NEWS_2a", "NEWS_2b", "SAV_Cz8", "SAV_C*6",
          "SAV_C*2", "SAV_C*2z60", "SAV_Z60", "NEWS_2old", "NEWS"),
  Descript = c("2005 - Emergent ", "2005 - Pre-Existent Wetland", "2005 - SAV ", "2005 - PSTA ",
               "2005 - Reservoirs - with depth penalty", "2009- SAV - modified depth penalty",
               "Non Emergent Wtld System ( SAV -> Periph. as C--> 0 ); Recommended for Designing SAV or PSTA Systems",
               "NEWS Tuned to Modern Data", "WCA2A Emergent--> PEW", "Alternative WCA2A",
               "Adjusted K1,k2,k2 vs. Z,c", "2005 - Reservoirs", "March 2002 EMERG Calibration",
               "Periphyton", "Submersed Macrophytes, Best Fit for All SAV Datasets",
               "Non Emergent Wtld System ( SAV -> Periph. as C--> 0 ); Recommended for Designing SAV or PSTA Systems",
               "SAV, Best Fit for ENRP Cell 4 Only, 1998-1999, Over-predicts Concentrations at C > 30 ppb",
               "2005 - Reservoirs", "March 2002 EMERG Calibration", "Periphyton",
               "Submersed Macrophytes, Best Fit for All SAV Datasets", "Non Emergent Wtld System ( SAV -> Periph. as C--> 0 ); Recommended for Designing SAV or PSTA Systems",
               "SAV, Best Fit for ENRP Cell 4 Only, 1998-1999, Over-predicts Concentrations at C > 30 ppb",
               "2005 - Reservoirs - with depth penalty conservative", "2012- Placeholder; C*=2,  Z1=60",
               "2005 - SAV Calibration Reduced by 20% (Lower CL of 2005 Calib)",
               "Updated in 2012  to Improve Simulation in Low P Range (Placeholder)",
               "NEWS with 2005 SAV & PSTA Calibrations", "NEWS with 2005 SAV & PSTA Calibrations",
               "NEWS with 2005 SAV & PSTA Calibrations", "Updated in 2012  to Improve Simulation in Low P Range (Placeholder)",
               "Updated in 2012  to Improve Simulation in Low P Range (Placeholder)",
               "Updated in 2012  to Improve Simulation in Low P Range (Placeholder)",
               "Updated in 2012  to Improve Simulation in Low P Range (Placeholder)",
               "Updated in 2012  to Improve Simulation in Low P Range (Placeholder)",
               "NEWS with 2005 SAV & PSTA Calibrations",
               "Non Emergent Wtld System ( SAV -> Periph. as C--> 0 ); Recommended for Designing SAV or PSTA Systems"),
  C0 = c(3L, 3L, 3L, 3L, 3L, 3L, 12L, 3L, 3L, 3L, 3L, 3L, 4L,
         4L, 12L, 12L, 4L, 3L, 4L, 4L, 12L, 12L, 4L, 3L, 2L,
         3L, 3L, 2L, 2L, 2L, 3L, 6L, 2L, 2L, 3L, 3L, 12L),
  C1 = c(22L, 22L, 22L, 22L, 150L, 22L, 22L, 22L, 22L, 11L, 22L, 150L,
         22L, 22L, 22L, 22L, 22L, 150L, 22L, 22L, 22L, 22L, 22L, 150L,
         22L, 22L, 22L, 22L,22L, 22L, 22L, 22L, 22L, 22L, 22L, 22L, 22L),
  C2 = c(300L, 300L, 300L, 300L, NA, 200L, NA, NA, 50L, 50L, 300L, NA,
         NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 300L, 300L,
         300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L,
         NA),
  Ks = c(16.8, 34.9, 52.5, 23.6, 5, 60, 128.7, 52.5, 50, 35, 34.93573193,
         3.1, 15.7,23.8, 128.7, 128.7, 80.1, 3.1, 15.7, 23.8, 128.7, 128.7,
         80.1, 2.5, 45, 42, 150, 100, 90, 80, 130, 160, 45, 45, 55, 80, 128.7),
  Z1 = c(40L, 40L, 40L, NA, 40L, 80L, 60L, 60L, 40L, 40L, 40L, 40L, 60L,
         0L, 60L, 60L, 60L, 40L, 60L, 0L, 60L, 60L, 60L, 40L, 60L, 40L,
         60L, 60L, 60L, 60L, 40L, 40L, 40L, 60L, 60L, 40L, 60L),
  Z2 = c(100L, 100L, 100L, 100L, 100L, 100L, NA, NA, 100L, 100L,
         100L, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 100L, 100L,
         100L, 100L, 100L, 100L, 100L, 100L, 100L, 100L, 100L, 100L, 100L,
         NA),
  Z3 = c(200L, 200L, 200L, 200L, 400L, 200L, NA, NA, 200L,
         200L, 200L, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 400L,
         200L, 200L, 200L, 200L, 200L, 200L, 200L, 200L, 200L, 200L, 200L,
         200L, NA),
  K1 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
         NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
         NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  C0_NEWS_SF = c(NA, NA, NA, NA, NA, NA, 4, 4, NA, NA, NA, NA
                 , NA, NA, NA, 4, NA, NA, NA, NA, NA, 4, NA,
                 NA, NA, NA, 8, 2, 2, 2, 8, NA, NA, NA, NA, 3, 4),
  C1_Peri = c(NA, NA, NA, NA, NA, NA, 22L, 22L, NA, NA,
              NA, NA, NA, NA, NA, 22L, NA, NA, NA, NA, NA, 22L, NA, NA, NA,
              NA, NA, 22L, 22L, 22L, NA, NA, NA, NA, NA, 22L, 22L),
  Ks_Peri = c(NA, NA, NA, NA, NA, NA, 23.8, 23.60209871, NA, NA, NA, NA, NA,
              NA,  NA, 23.8, NA, NA, NA, NA, NA, 23.8, NA, NA, NA, NA, NA, 30,
              30, 30, NA, NA, NA, NA, NA, 25, 23.8),
  Zx_Peri = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 40L, 40L, 40L,
              NA, NA, NA, NA, NA, 40L, NA),
  Sm = c(NA, NA, NA, NA, NA, NA, 400L, 400L, NA, NA,
         NA, NA, NA, NA, NA, 400L, NA, NA, NA, NA, NA, 400L, NA, NA, NA,
         NA, NA, 1200L, 1400L, 1000L, NA, NA, NA, NA, NA, 800L, 400L),
  Sb = c(NA, NA, NA, NA, NA, NA, 80L, 80L, NA, NA, NA, NA,
         NA, NA, NA, 80L, NA, NA, NA, NA, NA, 80L, NA, NA, NA, NA,
         NA, 80L, 80L, 80L, NA, NA, NA, NA, NA, 100L, 80L),
  MinDepth = c(35,
               38, 62, 13, 90, 62, 9, 9, 35.33941113, NA, NA, 0, 18, 9,
               30, 9, 30, 0, 18, 9, 30, 9.000000358, 30, 90, 90, 62, 62,
               13, 13, 13, 62, 62, 62, 62, 62.05025819, 12.50717716, 9.000000358),
  MaxDepth = c(76, 66, 87, 60, 304, 87, 109, 109, 76.38453248,
                  NA, NA, 0, 98, 70, 109, 109, 109, 0, 98, 70, 108.9720366,
                  108.9720366, 108.9720366, 304, 304, 87, 87, 87, 87, 87, 87,
                  87, 87, 87, 87.37593854, 87.37593854, 108.9720366),
  MinQW = c(26,
            69, 162, 3, 68, 162, 0, 0, 25.70540263, NA, NA, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 68, 68, 162, 162, 3, 3, 3, 162,
            162, 162, 162, 161.8217985, 2.616137163, 0),
  MaxQW = c(210,
            276, 374, 132, 1135, 374, 216, 216, 209.9445217, NA, NA,
            0, 345, 168, 216, 216, 216, 0, 345, 168, 215.8873472, 215.8873472,
            215.8873472, 1135, 1135, 374, 374, 374, 374, 374, 374, 374,
            374, 374, 373.7701997, 373.7701997, 215.8873472),
  MinConc = c(19.5,
              8, 14.9, 6, 50, 14.9, 6.3, 6.3, 19.53353075, NA, NA, 0, 12,
              6, 13, 6, 30, 0, 12, 6, 13.36412239, 6.276977539, 30, 50.3,
              50.3, 14.9, 14.9, 5.9, 5.9, 5.9, 14.9, 14.9, 15, 15, 14.85045332,
              5.865762869, 6.276977539),
  MaxConc = c(800, 110, 153, 56,
              1144, 153, 120, 120, 300, NA, NA, 0, 154, 25, 120, 120, 50,
              0, 154, 25, 120.458973, 120.458973, 50, 1144, 1144, 153,
              153, 153, 153, 153, 153, 153, 153, 153, 153.3227583, 153.3227583,
              120.458973),
  MinFreqZ_LT10cm = c("0%", "0%", "0%", "0", "0",
                      "0%", "", "", "0", "", "", "0", "", "", "", "", "", "0",
                      "", "", "", "", "", "0%", "0%", "0%", "0%", "0%", "0%", "0%",
                      "0%", "0%", "0%", "0", "0", "0", ""),
  MaxFreqZ_LT10cm = c("9%",
                      "13%", "0%", "0", "0", "0%", "", "", "0.089552239", "", "",
                      "0", "", "", "", "", "", "0", "", "", "", "", "", "0%", "0%",
                      "0%", "0%", "38%", "38%", "38%", "0%", "0%", "0%", "0", "0",
                      "0.379310345", ""),
  K_CV = c("0.20", "0.21", "0.16", "0.22",
           "0.45", "0.16", "0.20", "0.20", "0.20", "", "", "0.20", "0.20",
           "0.20", "0.20", "0.20", "0.20", "0.20", "0.20", "0.20", "0.2",
           "0.2", "0.2", "0.45", "0.16", "0.16", "0.16", "0.22", "0.22",
           "0.22", "0.16", "0.16", "16%", "0.16", "0.16", "0.22", "0.20"
  )

)
