## Overview

`fit_doublets.R` uses a Poisson-based maximum likelihood estimation (MLE) model to estimate the global doublet rate from the observed singlets and pairwise doublets counts from the genotype-based demultiplexing results.

## Prerequisites

* R (v3.5.0 or later)
* `dplyr`
* `tidyr`

## Usage

Run the script from the command line, pointing it to a compressed `.tsv` file containing your cell barcode assignments.

```bash
./fit_doublets.R path/to/assignments.tsv.gz
