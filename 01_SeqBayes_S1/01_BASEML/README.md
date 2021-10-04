# Running `BASEML`
Due to the large molecular data available for the 72 mammal taxa, using the exact likelihood
to compute the divergence times in a Bayesian framework would be computationally very expensive.
Instead, we can use an approximation to calculate the likelihood implemented in the `MCMCtree`
dating software ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)),
which is computationally more efficient.

Before proceeding with Bayesian inference, we need to calculate the gradient and the Hessian
for each of the partitions we previously generated when ordering the genes from slow- to fast-evolving. 
This step is required because the gradient and the Hessian are used to approximate the likelihood in `MCMCtree`.

We have divided this step into two main tasks:  

   1. **Preparing the rooted and calibrated trees**   
   The steps followed to generate the rooted and calibrated trees can be found in the [`01_trees`](01_trees) directory.   
   2. **Obtaining the gradient and the Hessian**   
   The steps followed to obtain the gradient and Hessian for each of the four partitions under each of the seven tree
   hypothesis can be found in the [`02_Hessian`](02_Hessian) directory.

	 