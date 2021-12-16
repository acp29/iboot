# statistics-bootstrap

## A statistics package with a variety of bootstrap resampling tools

This package of functions can be used to estimate uncertainty (confidence intervals) and test hypotheses (*p*-values) using bootstrap. Variations of the bootstrap are included that improve the accuracy of bootstrap statistics for small samples and samples with complex dependence structures. The package is forked from the GitHub repository [iboot](https://github.com/acp29/iboot).

## Installation
`pkg install "https://github.com/gnu-octave/statistics-bootstrap/archive/refs/heads/master.zip"`
 
 Find out basic information about th package by typing: `pkg describe -verbose statistics-bootstrap`

## Usage
`pkg load statistics-bootstrap`

### Functions

* `bootstrp` performs bootstrap resampling 
* `bootci` calculates confidence intervals using bootstrap resampling
* `ibootci` calculates confidence intervals (calibrated) by iterated bootstrap resampling 
* `ibootp` calculates a two-tailed *p*-value for hypothesized value of the statistic using bootstrap
* `iboottest` convenience function to perform one sample or paired-sample bootstrap tests (two-tailed)
* `iboottest2` convenience function to perform two-sample, unpaired bootstrap test (two-tailed)
* `bootperm` calculates a *p*-value for a one-sample bootstrap variant of the permutation test.
* `bootperm2` calculates a *p*-value for a two-sample (unpaired) bootstrap variant of the permutation test.
* `bootmode` uses bootstrap to evaluate the likely number of real modes in a distribution
* `plotboot` plots an overlay of a histogram, kernel density estimate and interval limits from bootstrap statistics

At the Octave command prompt, type `help function-name` for more information about the function and it's usage.
