# DV-axis projection

This repository contains a small package designed to assess the relation between
an observation's position along the dorsal-ventral axis in the spinal cord
and the value of a feature of interest (FOI), as presented in the manuscript
_"Decoding the development of the human spinal cord by multi-omics"_.

## Installation
In order to install the package, we recommend you to use the
[devtools](https://www.r-project.org/nosvn/pandoc/devtools.html) library. Once
`devtools` is installed simply execute:

```R
devtools::install_github("almaan/axis-projection")
```

And you should be able to load the package using the command:
```R
library(axis.projection)
```

## Usage
The package is designed to be compatible with objects created using the
[STUtility](https://ludvigla.github.io/STUtility_web_site/) package, which is
used throughout the analysis of the data in the aforementioned manuscript.
However, many of the functions can be used for stand-alone purposes if that's
desired.

### Examples
We provide a Rmarkdown file (both the raw `.Rmd` and a knitted `.html`) file
that illustrates how the package can be used and applied to the data we present
in the manuscript. These files are found in the `examples` folder, they are 
named `example.ext`  where `ext` is either of `Rmd` and `html`.

## Documentation
Documentation is basically non-existent at the moment, as this package is
released more with the intention to make the analysis reproducible rather than
with the aim to provide a versatile tool to be implemented in all workflows.


## Contributions
All code was written by [almaan](https://github.com/almaan/axis-projection), if
you have any problems using the package or questions pertaining to its use,
please feel free to open a GitHub issue on this page.
