# microranger

`microranger` contains wrapper functions for random forest classification from the [`ranger` package](https://github.com/imbs-hl/ranger).
It is intended to be used to perform classification of microbial communities. It allows to train algorithms using several parameter values, and aggregating the microbial community table to different taxonomic levels to identify the best set of conditions to make the classification.

You can install the package using [`devtools`](https://devtools.r-lib.org/):

```
devtools::install_github("marccamb/microranger")
```
