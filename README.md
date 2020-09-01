# An apparent paradox

A work in progress in semi-supervised generative classification. 

Based on the paper [An Apparent Paradox: A Classifier Trained from a Partially Classified Sample May Have Smaller Expected Error Rate Than That If the Sample Were Completely Classified](https://arxiv.org/abs/1910.09189), by Daniel Ahfock and Geoff Mclachlan.

### About the SAR Data

The SAR data products used as input in the script [SAR fit.R](SAR fit.R) are [cov_chol_d10_m1_cc.Rds](cov_chol_d10_m1_cc.Rds) and [labsp_d10.Rds](labsp_d10.Rds).

Original SAR Data is available on [this Zenodo archive](https://zenodo.org/record/4008883). From there, the 2x2 log-Cholesky transformed 1-lag autocovariance matrix of each pixel (over time) in the sequence of SAR images was obtained, after down-sampling using mean aggregation within 10x10 patches of pixels. Each 2x2 log-Chol matrix was then transformed to a 3-tuple, using equation (6) from the paper [here]( https://arxiv.org/pdf/2008.03454.pdf). This processing was done within the script [R/SAR_application.R](https://github.com/frycast/kmspd/blob/master/R/SAR_application.R) in the repo [frycast/kmspd](https://github.com/frycast/kmspd).

### Tools

[`{rsar}`](https://github.com/frycast/rsar) :package: provides tools for working with radar images in `SAR_matrix` format.

[`{EMMIX}`](https://github.com/suren-rathnayake/EMMIXmfa) :package: provides tools for fitting finite mixture models via the EM algorithm.
