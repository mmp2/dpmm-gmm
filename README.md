# dpmm-gmm
Dirichlet Process Mixtures for Generalized Mallows Models
Efficient C/Matlab MCMC sampling for Dirichlet Process Mixtures of
Generalized Mallows Models described in

Marina Meila and Harr Chen
Bayesian non-parametric clustering of ranking data
IEEE-TPAMI, to appear 2015

Please cite this paper if you are using this code.

Authors: Harr Chen harr@gmail.com
	 Marina Meila mmp@stat.washington.edu

Generalities
============
v4/     the Slice-Gibbs and Beta-Gibbs sampler
v5mmp/  the Exact-Beta-Gibbs sampler

Beta-Gibbs is fast and exact for top-t rankings where t < n-10,
	    n being the length of a complete permutation

	    for t closer to n, Beta-Gibbs is an approximate sampler

Exact-Beta-Gibbs is an exact sampler for the cases when Beta-Gibbs is not,
		 when the parameter t0 is large enough. Experimentally,
		 we found t0 = 11 to be sufficient.
		 
		 It is much slower than Beta-Gibbs, sometimes even slower 
		 than Slice-Gibbs. This is measured per iteration. But mixing
		 is much better than for Slice-Gibbs, so it should require 
		 fewer iterations.

Slice-Gibbs     exact, slower in time per iteration, slower mixing. Not really
		recommended in any situation.

		Was introduced as the "straw-man" for the other two.


Set up
======
1. download the code (v4 or v5mmp)

2. compile

mex -largeArrayDims sample_model.c
mex -largeArrayDims compute_pi_R.c

for v5mmp replace sample_model.c with sample_model_t0.c

To eke a bit more performance do:

mex -largeArrayDims -v CFLAGS='$CFLAGS -O4 -mtune=native -march=native -pipe -ffast-math -mfpmath=sse -msse4 -m64 -ansi -D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread' sample_model.c


3. choose running parameters 

define.h --> make numbers larger to see less test

4. call the code

The main entry point function for the sampler is inference.m. The comments 
should describe basically how to run it, but there's also an example script
 small_test.m that shows a very simple toy example. 

Note: the permutations in this system are 0-based, not 1-based
      the theta parameters from the paper are denoted rho

***************************************************************************
Bug reports by mmp:

1. only in v5mmp, t0 should NOT be larger than tmax, or memory overflows will occur in resample_rho_beta_t0 and possibly other places.

2. resample_sigma_rho_t0 (inherited from resample_sigma_rho)

rho is sampled up to g_max_t in every cluster, even if in that cluster
none of the data have length g_max_t. Haven't explored the
consequences yet. If g_S[ idx+i ] is defined and 0 for the extra ranks, nothing drastic will happen, as rho will be sampled from the prior. But the extra rho's will affect the sampling of sigma 

v4
-----------------
Bug report by mmp: sample_model.c, around line 400, where *randLengthValue is calculated. I think that for BETA_SLICE g_sigma_gibbs should be g_rho_slices

