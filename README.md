[![Build Status](https://travis-ci.com/f-puig/R.ComDim.svg?branch=main)](https://travis-ci.com/f-puig/R.ComDim)

# ComDim for evaluating replicate cross-platform variability and batch influence

## Introduction: the ComDim method
ComDim (also known as CCSWA) is an unsupervised multi-block method that aims to simultaneously consider multiple data tables to find the latent components that are common to all the tables as well as those that are specific to each data table, along with the contribution of each of the tables to each of these components. ComDim determines a common space describing the dispersion of the samples in all the blocks, each block having a specific weight (__salience__) associated with each dimension in this common space. Significant differences in the saliences for a given dimension reflect the fact that the dimension contains different amounts of information coming from each block. In addition to the saliences, __Local loadings__ for each analyzed block and two different sets of scores are obtained. The first set corresponds to the __Local scores__ for each analyzed block while the second set is composed of the __Global scores__, common to all the blocks.

## Functions
To successfully extract all the potential of the ComDim method, 3 functions coded in R are proposed:
* __BuildMultiBlock.R__: To merge several single data-blocks into a multi-block data set. 
* __splitRB.R__: To split one or more blocks into several smaller blocks, corresponding each new block to one batch.
* __comdim_PCA.R__: This function applies the ComDim algorithm on the multi-block object resulting from __BuildMultiBlock.R__ or from __splitRB.R__.

## Install and load R.ComDim package
```{r install, echo = TRUE}
  if (!require("devtools")) install.packages("devtools")
  library("devtools")
  install_github("f-puig/R.ComDim")
  
  # Load R.ComDim
  library("R.ComDim")
```

## References
* Puig-Castellví, F.; Jouan-Rimbaud Bouveresse, D.; Mazéas, L.; Chapleur, O.; Rutledge, D. N. Rearrangement of incomplete multi-omics datasets combined with ComDim for evaluating replicate cross-platform variability and batch influence. Under revision. 
* Qannari, E. M.; Courcoux, P.; Vigneau, E. Common Components and Specific Weights Analysis Performed on Preference Data. Food Qual. Prefer. 2001, 12 (5–7), 365–368. https://doi.org/10.1016/S0950-3293(01)00026-X.
* Mazerolles, G.; Hanafi, M.; Dufour, E.; Bertrand, D.; Qannari, E. M. Common Components and Specific Weights Analysis: A Chemometric Method for Dealing with Complexity of Food Products. Chemom. Intell. Lab. Syst. 2006, 81 (1), 41–49. https://doi.org/10.1016/J.CHEMOLAB.2005.09.004.
* Claeys-Bruno, M.; Béal, A.; Rutledge, D. N.; Sergent, M. Use of the Common Components and Specific Weights Analysis to Interpret Supersaturated Designs. Chemom. Intell. Lab. Syst. 2016, 152, 97–106. https://doi.org/10.1016/j.chemolab.2016.01.014.

