#include "mex.h"
#include "matrix.h"
#include "defines.h"
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "utility.c"

/*
EQUIVALENT MATLAB HEADER:

function pi_R = compute_pi_R(inv_pi, t, n)
% Produces compact sufficient statistics for the given permutations
% inv_pi:   permutations in ordering format (zero-based), a uint32(t, N)
%           matrix
% n:        number of items
% Output:   uint32(n - 1, t, N) matrix;
%           each column specifies the linear index of a 1 in the R matrix
%           specified by the next two columns;
%           indices are zero-based, and x is variable

*/

void mexFunction(int lhsCount, mxArray* lhs[], int rhsCount, const mxArray* rhs[])
{
	mxArray* pi_R;
	unsigned int* pi_R_data;
	unsigned int* pi_R_j_data;
	unsigned int pi_R_j_idx;
	unsigned int* inv_pi = (unsigned int*)mxGetData(rhs[0]);
	size_t max_t = mxGetM(rhs[0]);
	size_t N = mxGetN(rhs[0]);
	unsigned int* t = (unsigned int*)mxGetData(rhs[1]);
	unsigned int n = (unsigned int)mxGetScalar(rhs[2]);
	char* already_ranked = mxMalloc(n * sizeof(char));
	unsigned int* this_inv_pi;

	unsigned int i, j, k;
	unsigned int row;

	mwSize dims[3];
	dims[0] = n - 1;
	dims[1] = max_t;
	dims[2] = N;

	pi_R = mxCreateNumericArray(3, dims, mxUINT32_CLASS, mxREAL);
	pi_R_data = (unsigned int*)mxGetData(pi_R);
	lhs[0] = pi_R;

	for (i = 0; i < N; ++i)
	{
		this_inv_pi = inv_pi + index2d(0, max_t, i);
		memset(already_ranked, 0, n * sizeof(char));
		for (j = 0; j < t[i]; ++j)
		{
			/* Get this sparse vector of R */
			pi_R_j_data = pi_R_data + index3d(0, n - 1, j, max_t, i);
			pi_R_j_idx = 0;

			/* Fill vector with proper entries */
			/* Only row will be inv_pi[j] */
			row = this_inv_pi[j];
			already_ranked[row] = 1;
			for (k = 0; k < n; ++k)
			{
				/* If not already ranked, and not current item being ranked, mark in vector */
				if (already_ranked[k] == 0)
				{
					pi_R_j_data[pi_R_j_idx++] = k * n + row;
				}
			}
			assert(pi_R_j_idx == n - j - 1);
		}
	}

	mxFree(already_ranked);
}