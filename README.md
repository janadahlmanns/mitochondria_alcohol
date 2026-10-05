# mitochondria_alcohol: mitochondria positioning

Image-analysis code, written in R, for analyzing the localization of mitochondria in fluorescence microscopy images of neurons under different conditions (including an alcohol condition).

## Interactive analysis tool

`Analysis/distance_intensity_profile/` contains a **Shiny app** that guides the experimenter through the analysis of one image at a time:

1. Load a multi-channel `.tif` image.
2. Click into the image to mark the nucleus and three section lines through the cell body (the longest axis through the soma, and two lines at about +/-45 degrees).
3. The app extracts the **fluorescence intensity profile along the distance from the nucleus** for each line and channel (`EBImage`). Brightness and contrast sliders change only the display, never the analyzed values.
4. Assign the cell to an experimental condition and save the result.

Results are stored per cell as `.rds` files in `Results/Mitosections/`.

## Branch-point counts

`Analysis/branch_points/counts.csv` holds counts of branch points with and without mitochondria per slide and condition, and `Results/branch_points/plot_results.R` plots the ratio of mitochondria-positive branch points per condition.

## Run

```r
install.packages(c("shiny", "shinythemes", "shinyjs", "tidyverse"))
BiocManager::install("EBImage")
shiny::runApp("Analysis/distance_intensity_profile")
```

`Results/old_*` keeps earlier versions of the analysis for traceability.
