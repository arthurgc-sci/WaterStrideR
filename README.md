## WaterStrideR: a R package for morphological feature extraction on Gerroidea-like arthropods
Allows for acquisition of biologically relevant water strider morphological features from large amounts of image data directly in R. Treatment of pictures with lots of individuals in scientific contexts, e.g. laboratory breeding experiments or large field data.
#### Main features:
- For water striders or similar insects: hind leg segmentation and measurements with joint landmarking, using classic computer vision methods.
– For Microvelia longipes: sex and presence of wings prediction, using machine learning methods.
- Automatic scale acquisition on images following simple red graph paper protocol.

## Requirements
This package requires **R ≥ 3.5** and the following packages:

- **Imports:** dplyr, magrittr, imager, mclust, Momocs, pracma, MASS, progress, rootSolve  
- **Suggests:** knitr, rmarkdown, testthat (≥ 3.0.0)

## Installing
WaterStrideR is available on **GitHub**. To install, you will first need to have the **remotes** package installed. Run all the commands below for a first installation:

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
                        predict_sex_wing = FALSE)
summary(results)
```
*options for creation of an output folder including diagnostic figures and data as .csv file, manual interactive scaling, and no prediction of sex and presence of wings*

## Documentation

Workflow, use case, and tuning are detailed in the packages vignettes:
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
Contributions are welcome! However, this package is provided as-is for a specific research use case.  
Bug reports and pull requests are welcome but may not be actively reviewed or integrated.
- Read contribution guidelines [here](https://github.com/arthurgc-sci/WaterStrideR/CONTRIBUTING.md)

## Testing

Run package tests in R with:
```R
devtools::test(pkg = "WaterStrideR")
```

## Citation

If you use WaterStrideR in your research, please cite:
```
[Waiting JOSS publication]
```

## Acknowledgments

Development of this package was supported by Abderrahman Khila from Institut de Génomique Fonctionnelle de Lyon (IGFL), ENS de Lyon.
Example datasets were provided by Abderrahman Khila (IGFL) and Claudia Pruvôt (IGFL).

