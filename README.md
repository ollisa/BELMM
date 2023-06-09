# BELMM: Bayesian Estimation of Latent Mixture Models

### Description
BELMM is a framework for analyzing, clustering, and modelling time-series data in a Bayesian setting. The basic idea of the approach is to treat the individual mixture components as realizations of latent random walk processes. The models are estimated using [Stan](https://mc-stan.org/docs/stan-users-guide/index.html) and the model selection is done using Reversible jump Markov chain Monte Carlo implemented in [CU-MSDSp](https://github.com/jtchavisIII/CU-MSDSp). 

### Purpose
This is a repository for supplementary material for "BELMM; Bayesian model selection and random walk smoothing in time-series clustering" (Olli Sarala, Tanja Pyhäjärvi, Mikko J. Sillanpää). For details, please refer to the paper and the associated attachment. Code for data simulation, data processing, and analysing and plotting the results, see *BELMM_R_file.Rmd*. The Stan implementations of the models can be found under *Stan_models* folder. 

### Example figures
Posterior estimates of the centres and the realized assignments:
<img src="figs/github_readme_fig.png" align="center" width="1440" />

A Rectified model posterior distribution:
<img src="figs/rectified_mpd.png" align="center" width="400" />

### References
[Paper](link when available). 
