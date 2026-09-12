## WaterStrideR: a R package for morphological feature extraction on Gerroidea-like arthropods

This package enables automated extraction of biological traits associated with growth, reproduction, and fighting behaviour in water striders from large amounts image data directly in R. The entire process is fully integrated in a single tunable pipeline combining traditional computer vision, statistical modelling and supervised classification.

#### Main features:
- For water striders and similar insects: segmentation, landmarking, and measurement of the body and the hind legs segments.
- For Microvelia longipes: prediction of sex and wing presence.
- Automatic scale acquisition on images that contain a piece of red graph paper.

## Requirements
This package requires **R ≥ 3.5** and the following packages:

- **Imports:** dplyr, magrittr, imager, mclust, Momocs, pracma, MASS, progress, rootSolve  
- **Suggests:** knitr, rmarkdown, testthat (≥ 3.0.0)

## Installing
WaterStrideR is available on **GitHub**. To install it, you will first need to have the **remotes** package installed. Run all the commands below for a first installation:

```console
install.packages("remotes")
remotes::install_github("arthurgc-sci/WaterStrideR")
library(WaterStrideR)
```

## Quickstart
Run the entire analysis pipeline with a single tunable function with custom outputs as showcased below.
You will first need to input the path to either your image **or** the folder containing all the images you want to analyse.

```R
results <- gRunPipeline(your_path,
                        write_output = TRUE,
                        return_df = TRUE,
                        auto_scale = FALSE,
                        predict_sex_wing = TRUE) #for Microvelia longipes only
summary(results)
```
*The above line creates an output folder that will include diagnostic .png figures and a .csv file containing the data with prediction of sex and presence of wing. You will first be prompted to perform interactive manual scaling in a pop-up window.*

## Documentation

Workflow, use case, and tuning are detailed in the vignettes of the package:
#### Quickstart guide, setup and installation:
```R
vignette("getting_started", package = "WaterStrideR")
```
#### Output options and interpretation 
```R
vignette("interpreting-outputs", package = "WaterStrideR")
```
#### Tuning parameters for your specific acquisition protocol
```R
vignette("parameter_tuning", package = "WaterStrideR")
```

## Contributing
Contributions are welcome. WaterStrideR was initially developed for a specific species and imaging setup, but it may provide a very useful starting point for other high-throughput arthropod phenotyping projects. If you are working with large image datasets and would like to use WaterStrideR as a starting point for developing a mass-phenotyping pipeline, please get in touch!
- Read contribution guidelines [here](https://github.com/arthurgc-sci/WaterStrideR/tree/main/CONTRIBUTING.md)

## Testing

Tests are not built during installation.  
Run them locally from the source directory:

```r
devtools::test()
```

## Citation

If you use WaterStrideR in your research, please cite:
```
[Waiting JOSS publication]
```

## Acknowledgments

Development of this package was supported by Nicolas Goudemand and Abderrahman Khila, both from IGFL, ENS de Lyon.
Example datasets were provided by Abderrahman Khila (IGFL) and Claudia Pruvôt (IGFL).

