# LOCATA I/O Framework

# Modified by Shuai Sun for DoA estimation using LAP algorithm

MATLAB framework to read datasets, run baseline algorithm, and write results to file.

## Final release of the LOCATA dev & eval datasets

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3630471.svg)](https://doi.org/10.5281/zenodo.3630471)

## See also: LOCATA Evaluation Framework

MATLAB framework for evaluation of participants' submissions:
https://github.com/cevers/sap_locata_eval.git


Shuai SUN
One issue is the audio data is recorded much faster than the ground truth is provided. So I 
need to properly find what the delay is when a group truth point is given. Not properly aligned
1. Extract_... was used to extract audio data for Chris to have a test
2. LS_azimuth provides an example
3. ls_dynamic contains the final code used for the paper
4. dynamic_... are some intermidiate scripts during the testing