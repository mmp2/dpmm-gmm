#include "mex.h"
#include "matrix.h"
#include <math.h>

/* Random permutation stuff */
mxArray* g_rand_perm_array;
double* g_rand_perm;
mxArray* g_perm_length;
double* g_perm_length_data;

/* Inverts permutation into dest */
inline void invert_pi(unsigned int* dest, unsigned int* src, size_t length)
{
	unsigned int i;
	for (i = 0; i < length; ++i)
	{
		assert(src[i] < length);
		dest[src[i]] = i;
	}
}

/* Computes values of s vector with respect to given center permutation */
/* n is length of sigma, t is length of inv_pi */
inline void relative_s(unsigned int* s, unsigned int* sigma, unsigned int* inv_pi, size_t t)
{
	unsigned int i, j, rank;
	
	for (i = 0; i < t; ++i)
	{
		/* Compute s_i == rank of inv_pi[i] after removing inv_pi[0, ..., i - 1] */
		rank = sigma[inv_pi[i]];
		s[i] = rank;

		/* Account for already observed entries */
		for (j = 0; j < i; ++j)
		{
			if (sigma[inv_pi[j]] < rank)
			{
				--s[i];
			}
		}
	}
}

/* Generates a random permutation of the specified length */
/* Return value should be destroyed */
inline unsigned int* randperm(size_t length)
{
	/* Random permutation */
	size_t i;
	unsigned int* perm = mxMalloc(length * sizeof(unsigned int));
	double* perm_double;

	if (length == 1)
	{
		perm[0] = 0;
	}
	else
	{
		*g_perm_length_data = (double)length;
		mexCallMATLAB(1, &g_rand_perm_array, 1, &g_perm_length, "randperm");

		perm_double = mxGetPr(g_rand_perm_array);
		for (i = 0; i < length; ++i)
		{
			perm[i] = saferound(perm_double[i]) - 1;
		}
		mxDestroyArray(g_rand_perm_array);
	}
	return perm;
}

