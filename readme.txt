*** The iboot package for Octave/Matlab *** 

This is a set of functions to estimate uncertainty (confidence intervals) and test hypotheses (p-values). Variations of the bootstrap are included that improve the accuracy of bootstrap statistics for small samples and samples with complex dependence structures.

INSTALLATION

Add to the MATLAB or OCTAVE path:
/path/to/iboot/
/path/to/iboot/helper/
/path/to/iboot/param/

Note: In Octave, for parallel implementation the path of the above folders must be saved to the octaverc file (e.g. with 'savepath'), not just added temporarily to the path using 'addpath'. If not, then the the CPU workers do not have access to the parboot function

Current version: 2.8.9.0

Please cite as:

Penn, Andrew Charles. iboot: Iterated Bootstrap for Small Samples and Samples with Complex Dependence Structures [https://github.com/acp29/iboot]. Zenodo, 2020, doi:10.5281/zenodo.3992392
