## How to contribute:
## Bug 🪳

### Is it really a bug?
If you think you've encountered a bug (no pun intended) it might be that your weird or absent outputs are simply the outcome of inadequate parametrization of `gRunPipeline()`. You can check the [parameter tuning vignette](https://github.com/arthurgc-sci/WaterStrideR/blob/main/vignettes/parameter-tuning.Rmd):
```R
vignette("parameter-tuning", package = "WaterStrideR")
```
You can find documentation about advanced tuning options in R in the parameter description of the `gPipeline()` function:
```R
help(gPipeline, package = "WaterStrideR")
#or more simply:
?gPipeline
```

### It really is a bug...
First of all thank you for the time spent checking the section above. Please note that this project is maintained on a limited basis. If your issue has never been reported in a [previous issue](https://github.com/arthurgc-sci/WaterStrideR/issues), open a [new issue](https://github.com/arthurgc-sci/WaterStrideR/issues/new) with a clear title and description along with as much detail as possible:
- Code samples to reproduce the issue
- Console error messages
- Input data (preferably as `.rds`) that reproduces the issue help greatly.
- Your R version and package version ( with `sessionInfo()` )

#### You fixed the bug
Directly open a [new pull request](https://github.com/arthurgc-sci/WaterStrideR/pulls) with a detailed description outlining how you implemented the solution.

## New Feature

### Suggestion
To request a feature or suggest an improvement, feel free to open a [new issue](https://github.com/arthurgc-sci/WaterStrideR/issues/new) describing your needs to offer a platform for discussions. Note that this package was designed for a specific research use case and that we currrently do not plan on developping new features.

### Feature addition
If you've created a new feature that you would like to see implemented in `WaterStrideR`: thanks a lot! You can consider:
- Opening a [new issue](https://github.com/arthurgc-sci/WaterStrideR/issues/new) to allow for public discussion and critiques.
- Then, open a [new pull request](https://github.com/arthurgc-sci/WaterStrideR/pulls).

New functions should be compatible with related components of the package and come with clear documentation and unit tests (if appropriate).
