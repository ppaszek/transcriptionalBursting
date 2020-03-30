
The codes here underpin the data analysis and theoretical arguments in

>J. Bagnall, W. Rowe, N. Alachkar, J. Roberts, H. England, C. Clark, M Platt, D. A .Jackson, M. R. Muldoon and P. Paszek (2020),*Variability of the inducible gene expression is constrained by reciprocal bursting characteristics along gene-specific linear trends.*


### Outline of the computations ###

Inference of mean-variance relationship was performed using a dataset from BMDCs incorporating 29 scRNA-seq experiments (each corresponding to a single Fluidigm C1 experiment with up to 96 cells) on the response time-course (at 0, 1, 2, 3 and 6 h) as well as additional perturbations such as treatment with IFNb, inhibition of paracrine secretion (chemical or physical on chip) or cell knockout for ifnar1 and STAT1 expression from (Shalek et al., 2014). 

>Shalek AK, Satija R, Shuga J, Trombetta JJ, Gennert D, Lu D, Chen P, Gertner RS, Gaublomme JT, Yosef N, Schwartz S, Fowler B, Weaver S, Wang J, Wang X, Ding R, Raychowdhury R, Friedman N, Hacohen N, Park H, May AP, Regev A. Single-cell RNA-seq reveals dynamic paracrine control of cellular variation.Nature. 2014 Jun 19;510(7505):363-9. doi: 10.1038/nature13437. Epub 2014 Jun 11.

We considered 813 genes that were induced by at least two-fold (compared to unstimulated cells) at the population level at any time point during the LPS stimulation [as identified in (Shalek et al., 2014)]. Visual inspection of the data revealed outliers in the linear regression fit, therefore, outlier removal method with Mahalanobis distance was used (with 0.05 threshold for outlier detection).
After removing low abundant genes (maximum mean expression <100 read counts) this resulted in 290 genes for the core TLR dataset (time-course) and 323 for the combined dataset (including perturbations). 
Bursting characteristics (based on moment estimators) for individual data points were fitted using linear regression and power functions (in semi log scale) when appropriate and presented as smooth curves. 
Robust regression (excluding data points with corresponding fitted residuals > 1.5 std(residuals)) was used to either remove noisy data (as expected in the scRNA-seq measurement) or remove individual data-point that did not affect the overall trend. 
Equations, fitted parameters, corresponding correlation coefficients and highlighted outliers are included in the encosed files. 
Fitting protocols were implemented in Python using R kernel, individual gene graphs were produced in MATLAB R2014a. 

### Contents of the repository ###

1)	Folder ‘_pre-processing’ includes:

•	Calculation of summary statistics (mean, variance, moment estimators of burst size and frequency, etc) for all TLR-induced genes and datasets using ‘summary_stats.ipynb’; raw data (read counts per gene per C1 chip) is stored in ‘induced_gene_raw_data.csv’, summary statistics are stored in ‘induced_genes_processed_data_all_genes_all_conditions.csv’.   

•	Application of Mahalonobis distance to remove measurement outliers (per gene) using ‘mahalonobis_1.ipynb’; output is stored in ‘induced_genes_outliers_removed_all_812_4runs.csv’
 
2)	Folder ‘core TLR analysis’ includes regression fitting of mean-variance relationship for core genes in the dataset (TLR2, TL3 and TLR4 stimulation). Fitting is performed using ‘my_mean_var_regression.ipynb’ using pre-processed data ‘induced_genes_outliers_removed_all_812.csv’ (after removing conditions with zero counts). Code in Matlab  (‘validation.m’) provides a detailed display of fitted data using robust regression.

3)	Folder ‘core TLR with perturbation’ includes curve fitting of different characteristics in the full dataset. These involve applications of robust regression as well as comparison with the core TLR dataset for mean-variance, mean-burst-size, mean-frequency as well as burst size variance and frequency-variance. Each subfolder is organised in analogous way to ‘core TLR analysis’.

