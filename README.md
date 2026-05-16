# Sliced Wasserstein distances and their application to regression problems
My Master's thesis in Statistics and Data Science, including R scripts.
See https://github.com/janhove/slicer for details about the `slicer` package.

## Thesis
The source (`.Rnw`) files for my thesis as well as the compiled PDF
are in the root directory. 

The `unified.bst` and `unified.csl` files in this directory implement
the _Unified Style Sheet for Linguistics_ referencing style, which I like.
They were created by Mark Dingsemanse and shared under a 
[CC BY-SA 3.0 license](https://creativecommons.org/licenses/by-sa/3.0/).

The other files in the root directory are licensed under the CC BY 4.0
license, see `LICENSE-thesis`.

## Scripts
* `toy_example.R` generates the data for the toy example of Chapter 4,
  computes the (sliced) Wasserstein distances, and fits the Gaussian process
  models.
* `toy_example_figures.R` takes the output of `toy_example.R` and generates
  the figures from Chapter 4.
* `anisotropy_distances.R` generates the data for the examples in Section 5.5.2
  and computes the distances.
* `anisotropy_models.R` fits models using these distances and draws the figures
  with the results.
* `brains.R` runs the analyses for Section 5.6.

The scripts are licensed under the MIT license, see `LICENSE-scripts`.

## Data
The data from Eckermann et al.'s (2021) study that are used for the examples
in Section 5.6 are available as part of a 64.8 GB zipped directory
from [zenodo](https://doi.org/10.5281/zenodo.5658994).
I extracted the relevant Excel files from this directory and made them
available here under `data/eckerman/` as allowed under the [CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/).

## References
Eckermann, Marina, Bernhard Schmitzer, Franziska van der Meer, 
Jonas Franz, Ove Hansen, Christine Stadelmann & Tim Salditt. 
2021. 
Three-dimensional virtual histology of the human hippocampus based on 
phase-contrast computed tomography. PNAS 118(48). e2113835118.
[doi:10.1073/pnas.2113835118](https://doi.org/10.1073/pnas.2113835118).
