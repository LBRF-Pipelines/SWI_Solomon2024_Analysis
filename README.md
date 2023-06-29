# SWI_Solomon2023_Analysis

This repository contains the supplemental material, data and analysis code for Solomon et al. 2023 examining the reliance of motor imagery on perceptual and motor systems using the size weight illusion.

## Experiment Code

The experiment code used to collect the data will be make available XX XXX, 2023.

## Analysis Requirements

All dependencies for these scripts can be installed by running the following line:

The scripts were run in R 4.1.2 for publication. They may work with older versions of R but are not guaranteed to function correctly.

All dependencies for these scripts can be installed by running the following line:
    
```r
install.packages(c("tidyverse", "gsignal", "slider", "brms", "tidybayes", "emmeans", "parametes", "model"))
```

## Analysis Code Usage

The raw data for the project can be found on the [Open Science Framework](https://osf.io/v45pq/). UPDATE THIS!

1. Extract the (`labview.zip`) in a created directory `SWI_Analysis/_Data/`.
2. Download the (`SWI_datasheet.csv` & `masterlist.csv`) in the directory `SWI_Analysis/_Data/`.
3. Open a new R session and set the working directory to the root `SWI_Analysis/` folder (or whatever you've renamed it to) using `setwd()` or the RStudio menu.
4. Run one of the following commands in the R terminal:

```r
source('./_Scripts/0_Import.R') # imports all the project data
source('./_Scripts/1_Analysis.R') # runs filters and analysis on data
source('./_Scripts/2_Troubleshooting.R') # An optional script used to visualize data at different and facilicate data cleaning
source('./_Scripts/3_Descriptives.R') # summarizes desctiptive results and demographics and creates relevant visualizations
source('./_Scripts/4_Modelling.R') # computes Bayesian models used to identify credible effects
```

Running the modelling script will run the import, analysis and descriptives scripts, so in most cases you just want to run the fifth line.

The code is an exact copy used for the publication. Future commits will be made to optimize the code's performance.

