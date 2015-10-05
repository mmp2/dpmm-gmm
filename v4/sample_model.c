#include "mex.h"
#include "matrix.h"
#include "defines.h"
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* -------------------------------------------------------------- */
/* See defines.h for definitions of variables below */

/* Globals for permutations */
unsigned int* g_inv_pi;
unsigned int* g_t;
unsigned int  g_n;

/* Sample information */
unsigned int* g_c;
unsigned int* g_sigma;
double*       g_rho;

/* Caches */
unsigned int* g_nc;
unsigned int* g_S;

/* Globals for priors */
double*       g_r_0;
double        g_nu_0;
double        g_alpha_0;

/* Precomputed values */
unsigned int* g_pi_R;
double*       g_beta_table;

/* Model specification */
int           g_model;

/* -------------------------------------------------------------- */
/* Global helper variables */

/* Number of permutations total */
mwSize        g_N;

/* Maximum value of g_t */
unsigned int  g_max_t;

/* Computed prior = log(alpha_0 * (n - i)!/n!) as for i = 1:max_t */
/* Note indexing starts at zero - so g_dp_0 = log(alpha_0 * (1/n)), etc. */
double*       g_dp_0;

/* Inverse permutations of g_sigma */
unsigned int* g_inv_sigma;

/* A list view of g_nc > 0 */
unsigned int* g_valid_clusters;
unsigned int g_num_valid_clusters;

/* Scratch space for s calculations */
unsigned int* g_current_s;

/* -------------------------------------------------------------- */
/* For c resampler */

/* Probabilities of new c's */
double* g_c_prob;

/* --------------------------------------------------------------- */
/* For sigma and rho resamplers */

/* List of data points in a cluster */
unsigned int* g_cluster_elements;
unsigned int g_num_cluster_elements;

/* R_j matrices (n * n * max_t entries) */
unsigned int* g_R_j;

/* R = \sum R_j */
double* g_R;

/* Probabilities of new sigma values for stagewise sampler */
double* g_R_column_sums;
double* g_sigma_prob;

/* List of ranks that are not in a sigma, bit (char) vector of items that are not in a sigma */
unsigned int* g_empty_ranks;
char* g_is_item_missing;

/* Variables for MATLAB betarnd callback */
mxArray* g_alpha_array;
double*  g_alpha_param;
mxArray* g_beta_array;
double*  g_beta_param;
mxArray* g_betarnd_params[2];

/* Variables for MATLAB resample_rho_slice callback */
mxArray* g_coeff_array;
double*  g_coeff_param;
mxArray* g_n_minus_j_array;
double*  g_n_minus_j_param;
mxArray* g_nu_array;
double* g_nu_param;
mxArray* g_last_rho_array;
double*  g_last_rho_param;
mxArray* g_burnin_array;
mxArray* g_resample_rho_slice_params[5];

/* --------------------------------------------------------------- */
/* General sampler parameters */

unsigned int g_iterations;
unsigned int g_startiterations;
unsigned int g_sigma_gibbs;
unsigned int g_rho_slices;

/* NOTE: temperature is broken right now and will trip an assert; SET TO ONE */
double g_inverseTemperature;

/* Log probability calculations */
mxArray* g_logProbArray;
double*  g_logProb;

/* Random values from MATLAB */
double* g_unifrand;
unsigned int g_randIdx;

/* Helper files */
#include "utility.c"
#include "permutation.c"
#include "resample_sigmarho.c"
#include "resample_c.c"
/*#include "compute_log_prob.c"*/

/* Constructs nc cache for the first time */
mxArray* construct_nc(void)
{
	unsigned int i, j, c;

	mwSize dims[2];
	mxArray* nc_array;
	unsigned int* nc;

	dims[0] = g_max_t;
	dims[1] = g_N;
	nc_array = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
	nc = (unsigned int*)mxGetData(nc_array);
	
	for (i = 0; i < g_N; ++i)
	{
		c = g_c[i];

		for (j = 0; j < g_t[i]; ++j)
		{
			++nc[index2d(j, g_max_t, c)];
		}
	}

	return nc_array;
}

/* Constructs S cache for the first time */
mxArray* construct_S(void)
{
	unsigned int i, j, c;

	mwSize dims[2];
	mxArray* S_array;
	unsigned int* S;

	dims[0] = g_max_t;
	dims[1] = g_N;
	S_array = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
	S = (unsigned int*)mxGetData(S_array);
	
	for (i = 0; i < g_N; ++i)
	{
		c = g_c[i];

		relative_s(g_current_s, g_sigma + index2d(0, g_n, c), g_inv_pi + index2d(0, g_max_t, i), g_t[i]);
		for (j = 0; j < g_t[i]; ++j)
		{
			S[index2d(j, g_max_t, c)] += g_current_s[j];
		}
	}

	return S_array;
}

int initialize(int lhsCount, mxArray* lhs[], int rhsCount, const mxArray* rhs[])
{
	unsigned int i;

	log_start("initialize");

	/* Data */
	g_inv_pi = (unsigned int*)mxGetData(rhs[R_INV_PI]);
	g_t = (unsigned int*)mxGetData(rhs[R_T]);
	g_n = (unsigned int)mxGetScalar(rhs[R_N]);

	g_N = mxGetN(rhs[R_INV_PI]);
	g_max_t = findmax(g_t, g_N);

	/* Data integrity checks */
	assert(g_N == mxGetM(rhs[R_T]));
	assert(g_max_t < g_n);
	assert(g_max_t == mxGetM(rhs[R_INV_PI]));

#ifdef DEBUG
	/* Verify all t's are < n (fully specified is n - 1) */
	for (int i = 0; i < g_N; ++i)
	{
		assert(g_t[i] < g_max_t && g_t[i] <= n - 1);
	}
#endif

	/* Samples */
	g_c = (unsigned int*)mxGetPr(lhs[L_C] = mxDuplicateArray(rhs[R_C]));
	g_sigma = (unsigned int*)mxGetPr(lhs[L_SIGMA] = mxDuplicateArray(rhs[R_SIGMA]));
	g_rho = mxGetData(lhs[L_RHO] = mxDuplicateArray(rhs[R_RHO]));

	/* Needed for cache construction */
	g_current_s = mxMalloc(g_max_t * sizeof(unsigned int));

	/* Caches */
	if (mxIsEmpty(rhs[R_NC]))
	{
		/* Construct nc cache for the first time */
		g_nc = (unsigned int*)mxGetData(lhs[L_NC] = construct_nc());
	}
	else
	{
		g_nc = (unsigned int*)mxGetData(lhs[L_NC] = mxDuplicateArray(rhs[R_NC]));
	}
	if (mxIsEmpty(rhs[R_S]))
	{
		/* Construct nc cache for the first time */
		g_S = (unsigned int*)mxGetData(lhs[L_S] = construct_S());
	}
	else
	{
		g_S = (unsigned int*)mxGetData(lhs[L_S] = mxDuplicateArray(rhs[R_S]));
	}

	/* Priors */
	g_r_0 = mxGetPr(rhs[R_R_0]);
	g_nu_0 = mxGetScalar(rhs[R_NU_0]);
	g_alpha_0 = mxGetScalar(rhs[R_ALPHA_0]);

	/* Precomputed values */
	g_pi_R = (unsigned int*)mxGetData(rhs[R_PI_R]);
	g_beta_table = mxGetPr(rhs[R_BETA_TABLE]);

	/* Model spec */
	g_model = saferound(mxGetScalar(rhs[R_MODEL]));

	/* Log probability of model */
	g_logProbArray = lhs[L_LOGPROB] = mxCreateDoubleScalar(0.0);
	g_logProb = mxGetPr(g_logProbArray);

	/* Precomputed term for new cluster for DP = log(alpha_0 * (n - i)!/n!) (with i = 1 corresponding to index 0) */
	g_dp_0 = mxMalloc(g_max_t * sizeof(double));
	g_dp_0[0] = log(g_alpha_0 / (double)g_n);
	for (i = 1; i < g_max_t; ++i)
	{
		g_dp_0[i] = g_dp_0[i - 1] - log((double)g_n - (double)i);
	}

	/* Sampler view variables */
	g_valid_clusters = mxMalloc(g_N * sizeof(unsigned int));
	g_num_valid_clusters = 0;
	for (i = 0; i < g_N; ++i)
	{
		if (g_nc[index2d(0, g_max_t, i)] > 0)
		{
			g_valid_clusters[g_num_valid_clusters++] = i;
		}
	}
	g_inv_sigma = mxMalloc(g_n * g_N * sizeof(unsigned int));
	for (i = 0; i < g_num_valid_clusters; ++i)
	{
		invert_pi(g_inv_sigma + index2d(0, g_n, g_valid_clusters[i]), g_sigma + index2d(0, g_n, g_valid_clusters[i]), g_n);
	}

	/* Other variables for resamplers */
	g_c_prob = mxMalloc(g_N * sizeof(double));
	g_cluster_elements = mxMalloc(g_N * sizeof(unsigned int));
	g_R_j = mxMalloc(g_n * g_n * g_max_t * sizeof(unsigned int));
	g_R = mxMalloc(g_n * g_n * sizeof(double));
	g_R_column_sums = mxMalloc(g_n * sizeof(double));
	g_sigma_prob = mxMalloc(g_n * sizeof(double));
	g_empty_ranks = mxMalloc(g_n * sizeof(unsigned int));
	g_is_item_missing = mxMalloc(g_n * sizeof(char));

	/* For calling betarnd */
	if (g_model == BETA_GIBBS)
	{
		g_alpha_array = mxCreateDoubleMatrix(g_max_t, 1, mxREAL);
		g_alpha_param = mxGetPr(g_alpha_array);
		g_beta_array = mxCreateDoubleMatrix(g_max_t, 1, mxREAL);
		g_beta_param = mxGetPr(g_beta_array);
		g_betarnd_params[0] = g_alpha_array;
		g_betarnd_params[1] = g_beta_array;
	}

	/* For calling resample_rho_slice */
	if (g_model == SLICE_GIBBS)
	{
		g_coeff_array = mxCreateDoubleScalar(0.0);
		g_coeff_param = mxGetPr(g_coeff_array);
		g_n_minus_j_array = mxCreateDoubleScalar(0.0);
		g_n_minus_j_param = mxGetPr(g_n_minus_j_array);
		g_nu_array = mxCreateDoubleScalar(0.0);
		g_nu_param = mxGetPr(g_nu_array);
		g_last_rho_array = mxCreateDoubleScalar(0.0);
		g_last_rho_param = mxGetPr(g_last_rho_array);
		g_burnin_array = mxCreateDoubleScalar(mxGetScalar(rhs[R_RHO_SLICES]));
		g_resample_rho_slice_params[0] = g_coeff_array;
		g_resample_rho_slice_params[1] = g_n_minus_j_array;
		g_resample_rho_slice_params[2] = g_nu_array;
		g_resample_rho_slice_params[3] = g_last_rho_array;
		g_resample_rho_slice_params[4] = g_burnin_array;
	}

	/* Stuff specific to this run of the sampler */
	g_iterations = saferound(mxGetScalar(rhs[R_ITERATIONS]));
	g_startiterations = saferound(mxGetScalar(rhs[R_STARTITERATIONS]));
	g_inverseTemperature = mxGetScalar(rhs[R_TEMPERATURE]);
	g_inverseTemperature = (g_inverseTemperature == 0.0 ? DBL_MAX : 1.0 / g_inverseTemperature);
	g_sigma_gibbs = (unsigned int)mxGetScalar(rhs[R_SIGMA_GIBBS]);
	g_rho_slices = (unsigned int)mxGetScalar(rhs[R_RHO_SLICES]);

	/* Temporary variables for permutations */
	g_perm_length = mxCreateDoubleScalar(0.0);
	g_perm_length_data = mxGetPr(g_perm_length);

	log_end("initialize");
}

void cleanup(void)
{
	mxFree(g_dp_0);
	mxFree(g_valid_clusters);
	mxFree(g_inv_sigma);
	mxFree(g_current_s);
	mxFree(g_c_prob);
	mxFree(g_cluster_elements);
	mxFree(g_R_j);
	mxFree(g_R);
	mxFree(g_R_column_sums);
	mxFree(g_sigma_prob);
	mxFree(g_empty_ranks);
	mxFree(g_is_item_missing);
	if (g_model == BETA_GIBBS)
	{
		mxDestroyArray(g_alpha_array);
		mxDestroyArray(g_beta_array);
	}
	if (g_model == SLICE_GIBBS)
	{
		mxDestroyArray(g_coeff_array);
		mxDestroyArray(g_n_minus_j_array);
		mxDestroyArray(g_nu_array);
		mxDestroyArray(g_last_rho_array);
		mxDestroyArray(g_burnin_array);
	}
}

void mexFunction(int lhsCount, mxArray* lhs[], int rhsCount, const mxArray* rhs[])
{
	unsigned int i, iters;
	unsigned int* order;

	/* Random numbers */
	/* Stuff for calling unifrnd */
	mxArray* unifrandArray;
	mxArray* zero = mxCreateDoubleScalar(0.0);
	mxArray* one = mxCreateDoubleScalar(1.0);
	mxArray* randLength = mxCreateDoubleScalar(0.0);
	double* randLengthValue = mxGetPr(randLength);
	mxArray* unifrndParams[4];

	unifrndParams[0] = zero;
	unifrndParams[1] = one;
	unifrndParams[2] = randLength;
	unifrndParams[3] = one;

	/*	fast_print("Entering MEX function sample_model\n");*/

	initialize(lhsCount, lhs, rhsCount, rhs);

	mexSetTrapFlag(1);

	for (iters = 0; iters < g_iterations; ++iters)
	{
		/* Random numbers for this iteration */
		*randLengthValue = (double)g_N; /* for c */
		if (g_model == BETA_GIBBS)
		{
			*randLengthValue += (double)(2 * g_max_t * g_N); /* for exact sigma sampling, including possibility of being called from c */
		}
		*randLengthValue += (double)(g_sigma_gibbs * g_max_t * g_N); /* for approx. sigma sampling, including all iterations */
		if (g_model == SLICE_GIBBS)
		{
			*randLengthValue += (double)(g_sigma_gibbs * g_max_t * g_N); /* for approx. sigma sampling of new clusters */
		}
		mexCallMATLAB(1, &unifrandArray, 4, unifrndParams, "unifrnd"); 
		g_unifrand = mxGetPr(unifrandArray);
		g_randIdx = 0;

		/* RESAMPLE c */
		log_start("resample_c");
		order = randperm(g_N);
		for (i = 0; i < g_N; ++i)
		{
			resample_c(order[i]);
			if ((i + 1) % C_NOTIFICATION_PERIOD == 0)
			{
				mexPrintf("Finished permutation %d\n", i + 1);
				log_end("resample_c");
				log_start("resample_c");
			}
		}
		mxFree(order);
		log_end("resample_c");

		/* RESAMPLE sigma */
		log_start("resample_sigmarho");
		order = randperm(g_num_valid_clusters);
		for (i = 0; i < g_num_valid_clusters; ++i)
		{
			resample_sigmarho(g_valid_clusters[order[i]]);
			if ((i + 1) % SIGMA_NOTIFICATION_PERIOD == 0)
			{
				mexPrintf("Finished cluster %d\n", i + 1);
				log_end("resample_sigmarho");
				log_start("resample_sigmarho");
			}
		}
		mxFree(order);
		log_end("resample_sigmarho");

		/* CLEANUP */
		mxDestroyArray(unifrandArray);

		/* LOG PROB */
		/* TODO: eventually... */
		/**g_logProb = compute_log_prob();
		mexPrintf("Logprob: %e\n", *g_logProb);*/
		*g_logProb = 0.0;

		assert(g_randIdx <= *randLengthValue);

		mexPrintf("Iteration complete: %u\n", (g_startiterations + iters + 1));
		mexPrintf("Number of clusters: %d\n", g_num_valid_clusters);
		drawnow();
	}
	
	mxDestroyArray(zero);
	mxDestroyArray(one);
	mxDestroyArray(randLength);
	cleanup();

	/*	fast_print("Leaving MEX function sample_model\n");*/
}
