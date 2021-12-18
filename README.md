# statistics-bootstrap

## Package maintainer
Andrew Penn (andy.c.penn@gmail.com)

## Package contributors
Andrew Penn

## Citations
If you use this package please include the following citation(s):

* Penn, Andrew Charles. (2020). iboot: Iterated bootstrap for small samples or samples with complex dependence structures. Zenodo. https://doi.org/10.5281/zenodo.3992392  

## A statistics package for Octave/Matlab providing a variety of bootstrap resampling tools

This package of functions can be used to estimate uncertainty (confidence intervals) and test hypotheses (*p*-values) using bootstrap. Variations of the bootstrap are included that improve the accuracy of bootstrap statistics for small samples and samples with complex dependence structures. The statistics-bootstrap package is forked from the GitHub repository [iboot](https://github.com/acp29/iboot).

## Requirements and dependencies

This package is known to be compatible with versions of Octave v3.2.4+ and Matlab v7.4.0+. It may be compatible with some earlier versions of Octave. Most features of of the package do not depend on other packages. An exception is the parallel computing options, which require either the Parallel Computing MATLAB Toolbox or the OCTAVE forge parallel package.

## Installation
 
To install (or test) the statistics-bootstrap package at it's existing location in either Octave or Matlab, follow these steps: 
 
 * Download the package. If it is a compressed file, decompress it
 * Open Octave or Matlab
 * `cd` to the package directory. (The directory contains a file called 'install.m')
 * type `install` at the octave command prompt
 
 To uninstall, follow the above steps but use the `uninstall` command
 
 Alternatively, users of more recent versions of Octave can install the package automatically with the following command:
 
 `pkg install "https://github.com/gnu-octave/statistics-bootstrap/archive/refs/heads/master.zip"`
 
 In Octave, you can find out basic information about the package by typing: `pkg describe -verbose statistics-bootstrap`  

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


