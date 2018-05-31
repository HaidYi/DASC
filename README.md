
# DASC

[![Build Status](https://travis-ci.org/HaidYi/DASC.svg?branch=master)](https://travis-ci.org/HaidYi/DASC)

## Overview

DASC is an R package used for identifying batches and classifying samples into different batches in a high dimensional gene expression dataset. The batch information can be further used as a covariate in conjunction with other variables of interest among standard bioinformatics analysis like differential expression analysis. 

## Installation

**_Dependencies_**

Before installing the package, please make sure the prerequisite packages described in `DESCRIPTION` file have been installed successfully.

To install the DASC package, start the terminal and enter:

On linux/MacOs
```bash
git clone https://github.com/HaidYi/DASC.git
R CMD build DASC
R CMD INSTALL DASC_*.tar.gz
```

## Usage

For concrete usages, please refer the *.html* vignettes in `inst/doc` file.

## Citation info

If you use DASC for your analysis, please cite our paper as here below.

```
@article{Yi2018Detecting,
    title={Detecting hidden batch factors through data-adaptive adjustment for biological effects},
    author={Yi, H. and Raman, A. T. and Zhang, H. and Allen, G. I. and Liu, Z.},
    journal={Bioinformatics},
    volume={34},
    number={7},
    pages={1141},
    year={2018},
}
```

## Questions ?

If you have any questions or issues? Please ask in the Issues or email me directly at HaidYi@mail.nankai.edu.cn


