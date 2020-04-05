# Transcriptional Bursting #

The codes here underpin the data analysis and theoretical arguments in

>J. Bagnall, W. Rowe, N. Alachkar, J. Roberts, H. England, C. Clark, M Platt, D. A .Jackson, M. R. Muldoon and P. Paszek (2020),*Variability of the inducible gene expression is constrained by reciprocal bursting characteristics along gene-specific linear trends.*

To investigate models of transcriptional regulation we calculated exact, time-dependent distributions of the number mRNA transcripts per cell using the Chemical Master Equation (CME) and an approach similar to that in

>M. Gómez-Schiavon, L.-F. Chen, A. E. West, and N. E. Buchler (2017), BayFish: Bayesian inference of transcription dynamics from population snapshots of single-molecule RNA FISH in single cells, *Genome Biology*, **18**:164. DOI: [10.1186/s13059-017-1297-9](https://doi.org/10.1186/s13059-017-1297-9)


### Outline of the computations ###
We estimated model parameters using a genetic algorithm based on MATLAB's `ga()` function in combination with `expmv()`, a MATLAB implementation of a fast matrix exponentiation algorithm developed by the authors of

>A. Al-Mohy and N. Higham (2011), Computing the action of the matrix exponential, with an application to exponential integrators, *SIAM Journal on Scientific Computing*, **33**:488–511. DOI: [10.1137/100788860](https://doi.org/10.1137/100788860)

The CME assigns nonzero (though eventually very small) probabilities to arbitrarily high mRNA counts, so we worked with distributions truncated at 2000 transcripts, a level that far exceeds any counts seen in our data.  The GA minimised the mean absolute difference between the theoretical (CME) and measured cumulative distribution functions (CDFs). CDFs were fitted using MATLAB's `fitdist()` function with an Epanechnikov kernel.

The 50 best model fits from independent GA runs were obtained for each condition (using a GA population size of 200, elite count of 2, 0.6 crossover factor, and the tournament selection function). Gene activation rates were constrained to lie below 0.2 1/min value, while the degradation rate was for TNFalpha transcripts was constrained to lie between 0.006 and 0.07 1/min (half-life between about 10 and 115 mins), while the degradation rate for IL1$\beta$ was constrained between 0.002 and 0.006 1/min (half-life between 115 and around 350 mins). We assumed two independent alleles per gene with the transcription rate no greater that 30 transcripts per minute per allele. This is in agreement with, for example

> N. Molina, D. M. Suter, R. Cannavo, B. Zoller, I. Gotic, and F. Naef (2013), Stimulus-induced modulation of transcriptional bursting in a single mammalian gene, *PNAS*, **110**:20563–20568. DOI: [10.1073/pnas.1312310110](https://dx.doi.org/10.1073/pnas.1312310110)

> B. Schwanhäusser, D. Busse, N. Li, G. Dittmar, J. Schuchhardt, J. Wolf, W. Chen, and M. Selbach (2011), Global quantification of mammalian gene expression control, *Nature*, **473**:337–342. DOI: [10.1038/nature10098](https://dx.doi.org/10.1038/nature10098)

>S. O. Skinner, H. Xu, S. Nagarkar-Jaiswal, P. R. Freire, T. P. Zwaka, and I. Golding (2016), Single-cell analysis of transcription kinetics across the cell cycle, *eLife*, **5**:e12175. DOI: [10.7554/eLife.12175](https://dx.doi.org/10.7554/eLife.12175)

> D. M. Suter, N. Molina, D. Gatfield, K. Schneider, U. Schibler, and F. Naef (2011),  Mammalian genes are transcribed with widely different bursting kinetics, *Science*, **332**:472–474. DOI: [10.1126/science.1198817](https://dx.doi.org/10.1126/science.1198817)


where rates as high as 2 to 10 transcripts min$\,^{-1}$ were reported for some genes. 

### Contents of the repository ###
The codes here allow one to fit smFISH count distributions using a genetic algorithm in MATLAB. There are two folders:

* `2_state_model_fitting`, which holds data and tools to fit a standard , two-state telegraph model to smFISH data about TNF$\alpha$ and
* `3_state_model_fitting`, which provides similar materials for a three-state model applied to smFISH data for IL1$\beta$.

In addition to the MATLAB files that specify the models, both folders include

 * A spreadsheet with a name of the form `*_smFISH.xlsx` that lists measured smFISH counts across a variety of experimental conditions.
 * A MATLAB script with a name of the form `checking_dist_*.m` that fits and plots a kernel density estimate of the CDF of the measured distribution of mRNA counts: such CDFs specify the objective function of the GA.
* A function with a name of the form `buildrateMat*.m` that takes a set of parameters for a model of transcription and constructs the associated (sparse) matrix of state-to-state transition rates.
* A script `ga_run.m` that runs the GA.
* A folder called `HighamAlMohly` that provides the MATLAB function `expmv()` and is essentially a copy of the GitHub repository [https://github.com/higham/expmv](https://github.com/higham/expmv).

Model-specific codes include;

* In the folder `2_state_model_fitting`, a MATLAB function `dist_fit.m` that takes a set of parameters for a model of transcription and computes the GA's objective function.
* In the folder `3_state_model_fitting`, a function `state3_IL1.m` that calculates the objective function for the 3-state model of transcription.





