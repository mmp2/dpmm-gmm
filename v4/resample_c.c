#include "mex.h"
#include "matrix.h"
#include <limits.h>
#include <float.h>

/* Adds or removes the given permutation (indexed by idx) from the global arrays of sufficient statistics */
inline void modify_suff_stats(unsigned int idx, bool add)
{
	unsigned int i, j;
	unsigned int c = g_c[idx];
	size_t base = index2d(0, g_max_t, c);
	int modifier = (add ? 1 : -1);

	/* Compute s vector of this permutation */
	relative_s(g_current_s, g_sigma + index2d(0, g_n, g_c[idx]), g_inv_pi + index2d(0, g_max_t, idx), g_t[idx]);

	/* Adjust nc and S_j values accordingly */
	for (i = 0; i < g_t[idx]; ++i)
	{
		g_nc[base + i] += modifier;
		g_S[base + i] += modifier * (int)g_current_s[i];
		assert(g_nc[base + i] < 0xEFFFFFFF);
		assert(g_S[base + i] < 0xEFFFFFFF);
	}

	/* Drop or add a cluster from valid list if necessary */
	if (!add && g_nc[base] == 0)
	{
		for (i = 0; i < g_num_valid_clusters; ++i)
		{
			if (g_valid_clusters[i] == c)
			{
				for (j = i; j < g_num_valid_clusters - 1; ++j)
				{
					g_valid_clusters[j] = g_valid_clusters[j + 1];
				}
				--g_num_valid_clusters;
				break;
			}
		}
	}
	else if (add && g_nc[base] == 1)
	{
		g_valid_clusters[g_num_valid_clusters++] = c;
	}
}

/* Resamples P(c_idx | everything else) */
void resample_c(unsigned int idx)
{
	unsigned int i, j, c, clust;
	double r_j_post, N_j_post;
	double rho;

	/* Current inverse permutation */
	unsigned int* inv_pi = g_inv_pi + index2d(0, g_max_t, idx);
	unsigned int t = g_t[idx];

	/* First empty cluster */
	unsigned int emptyC;

	/* Remove idx from caches and sufficient statistics */
	modify_suff_stats(idx, false);

	/* Compute all possible cluster assignments to pre-existing clusters */
	for (c = 0; c < g_num_valid_clusters; ++c)
	{
		clust = g_valid_clusters[c];

		/* Compute s vector against this centroid */
		relative_s(g_current_s, g_sigma + index2d(0, g_n, clust), inv_pi, t);

		/* Probability contribution from DP */
		g_c_prob[c] = log(g_nc[index2d(0, g_max_t, clust)]);

		/* The rest of the term is as follows for marginalized case:
		   r_j_post = S_j + r_j * nu
		   N_j_post = N_j + nu + 1
		   Prob = \sum_{j = 1 to t} log(N_j_post) - log(N_j_post + r_j_post + s_j) +
					  \sum_{i = 0 to s_j - 1} (log(r_j_post + i) - log(r_j_post + N_j_post + i))

		   For explicit case it's just the GMM probability */
		if (g_model == BETA_GIBBS)
		{
			for (j = 0; j < t; ++j)
			{
				/* Add ratio of Beta functions here */
				r_j_post = (double)g_S[index2d(j, g_max_t, clust)] + g_r_0[j] * g_nu_0;
				N_j_post = (double)g_nc[index2d(j, g_max_t, clust)] + g_nu_0 + 1.0;

				g_c_prob[c] += log(N_j_post / (N_j_post + r_j_post + g_current_s[j]));
				for (i = 0; i < g_current_s[j]; ++i)
				{
					g_c_prob[c] += log((r_j_post + (double)i) / (r_j_post + N_j_post + (double)i));
				}
			}
		}
		else
		{
			for (j = 0; j < t; ++j)
			{
				rho = g_rho[index2d(j, g_max_t, clust)];

				/* Normalization factor */
				g_c_prob[c] += log1p(-exp(-rho)) - log1p(-exp(-((double)g_n - j) * rho));

				/* GMM term */
				g_c_prob[c] -= rho * (double)g_current_s[j];
			}
		}
	}

	/* Find first empty cluster */
	for (c = 0; c < g_N; ++c)
	{
		if (g_nc[index2d(0, g_max_t, c)] == 0)
		{
			emptyC = c;
			break;
		}
	}
	assert(emptyC != UINT_MAX);

	/* Compute probability of new cluster */
	/* Probability contribution from DP and marginalized GMM */
	g_c_prob[g_num_valid_clusters] = g_dp_0[t - 1];

	/* Sample cluster assignment */
	c = sample_normalized_log(g_c_prob, NULL, g_num_valid_clusters + 1, g_inverseTemperature, g_unifrand[g_randIdx++]);
	g_c[idx] = (c < g_num_valid_clusters ? g_valid_clusters[c] : emptyC);

	/* Insert back into caches */
	modify_suff_stats(idx, true);

	/* Is this a new cluster? If so, then we need to re-estimate a centroid */
	if (g_c[idx] == emptyC)
	{
		resample_sigmarho(emptyC);
	}
}