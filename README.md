# TubAR

TubAR (Tuber Analysis in R) provides potato researchers with a simple means of quantitative phenotyping replacing previously used imprecise qualitative scales. Users collect images of sample tubers using a light box, and then use the package to quantify: tuber length, width, roundness, skinning percentage, redness, and lightness. Finally, the data is prepared for export and use in mapping, genomic selection, etc.

# Installing the package

Use the `install_github` function from the devtools package in R to install TubAR.
TubAR makes use of the EBImage and Biobase packages, which are available through BiocManager, not CRAN, so you will also have to install those packages before TubAR.
```{r, include=T, eval=F}
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("EBImage")
BiocManager::install("Biobase")
devtools::install_github("shannonlabumn/TubAR")
```
More information on using the package is available in the package vignette.
