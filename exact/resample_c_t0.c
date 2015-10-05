#include "mex.h"
#include "matrix.h"
#include <limits.h>
#include <float.h>


/* This function has numerical errors when a,b too large (eb
   ~1000). Now it is called from dlog_btilde which calculates the
   difference in log Btilde()'s avoiding numerical errors.
*/
double betatilde( double a, double b, double n ) {
  /* Computes ~Beta( a, b, n ) by numerical integration.
     a, b parameters
     g_grid01 = uniform grid on [0, 1] of g_ngridBeta intervals

     It follows that grid[ 0 ] = 0, grid[ ngrid ] = 1 and value at 1 is
     added outside the loop. The integrand = 0 at 0, = 1/(n+1)^(b-1) at 1
  */
  double bt = 0;  /* this accumulates the ~Beta value */
  int i;
  for( i = 1; i < g_ngridBeta; ++i ) {
    double x = g_grid01[ i ]; 
    double y = pow((1-x)/(1-pow(x,n+1)),b-1);
    bt += pow(x,a-1)*y;
    /*if( isdebug) printf( " %.3g ", bt );*/
  }
  bt += pow(n+1,-b+1)/2.;  /* includes extremities /2, trapeze method */
  bt *= g_grid01[ 1 ];     
  /* this is g_grid01[ 1 ]-g_grid01[ 0 ], since the second term is 0 */
  /*if( isdebug ) printf( "%.3g\n", bt );*/
  return( bt );
}

inline double aloga( double a ) { return( a*log(a)); }

/* Computes 
   log( betatilde( a+s,b+1,1)/betatilde( a,b,1 ) for a,b large by Laplace's 
   method.
   REQUIRES b>=a>0  and a,b LARGE (otherwise it's a bad approximation)
   POSTCONDITION: return value is always <=0 for correct inputs
*/
inline double dlog_btil1( int s, double a, double b ) {
  double b1 = b-1.;
  double a1 = a-1.;
  double dd = 0;
  if( b1>2*a1 ) {
    if( s == 0 )
      dd = -aloga(b)+aloga(b-a1)-log(b-a1)-log(1./a1-1./b)*0.5
	+aloga(b1)-aloga(b-a)+log(b-a)+log(1./a1-1./b1)*0.5;
    else
      dd = aloga(a)-aloga(b)-log(1./a-1./b)*0.5
	-aloga(a1)+aloga(b1)+log(1./a1-1./b1)*0.5;
    return dd;
  }
  else{
    /** bad gaussian approximation, not used
    double aha = 4.*a1-b1;
    if( s==0 )
      dd = (b*b/(aha-1)-b1*b1/aha)/8.-0.5*log(1.-1./aha)-log(2.)-0.375;
    else
      dd = (b*b/(aha+3.)-b1*b1/aha)/8.-0.5*log(1+3./aha)-log(2)+0.125;
    */
    if( s==0 ) 
      dd = (1.+1./(2*a1-b))*.5;
    else
      dd = (1.-1./(2*a-b))*.5;
    return log( dd ); /* linear approximation of e^f, must return log(dd) */
  }
}
/** approximates log(btilde(a,b,n)) by gaussian around maximum.
    REQUIRES: max f in (0,1). Otherwise the calling function dlog_btilde
    provides an approximation by derivative.
 */
double log_btiln( double a, double b, double n ) {
  double b1 = b-1.;
  int i;
  /* find zero of derivative of f where
     f=a*ln(u)-(b-1)[ln(1-u^(n+1))-ln(1-u)] */
  int niter = ceil( log( 1e8 )/log( 2 ));  /* error is halfed each step so number of iterations is ~log_2(tolerance) */
  double umin = 0., umax = 1., u = 0.5;
  double un;
  for( i=0; i<niter; i++) {
    un = pow(u,n);
    double fprim = a/u-b1*(1./(1.-u)-(n+1)*un/(1.-un*u));
      if( fprim > 0 ) 
	umin = u;
      else
	umax = u;
      u = (umin+umax)/2.;
  }   /* now u is the location of the maximum of f */
  /*printf( "u=%3g ", u );*/
  /*assert( u<1 ); /* otherwise, this function should not be called */
  /* -f''(u) */
  double neg_fsecund = a/u/u+b*(1./(1.-u)/(1.-u)-(n+1.)*(n*un/u+un*un)/(1-un*u)/(1-un*u));
  /* f(u) */
  double f = a*log(u)-(b1)*(log(1-un*u)-log(1-u));
  /* calculate log(Betatilde) */
  return( f-log(neg_fsecund/2./3.1416)*0.5 );
}

/* Calculates 
 log( betatilde( a+s,b+1,n)/betatilde( a,b,n) ) in a numerically  stable way.
 REMARKS:
   1. Switching to approximation is done curently at Hu=-703, but practically I have  seen betatilde() work up to values twice as negative.
   2. The approximation is closed form, so it's much faster than computing the
   integral. May also be more accurate, (when the max is in (0,1)). One could
   consider lowering Hu instead of making it as high as it can get.
*/

double dlog_btilde( double s, double a, double b, double n ) {
  double Hu = aloga(a)+aloga(b)-aloga(a+b); /* Hu<0 */ 
  /*printf( "Hu=%.0f %.0f ", Hu, log( DBL_MIN )+5.0 );*/
  /* breakpoint currently just over a=b=500 */
  if((a<1) || (Hu > log( DBL_MIN )+5.0 ))  /* direct integration OK*/
    return( log( betatilde( a+s,b+1.,n) )-log( betatilde( a,b,n )));
  else { /* make approximations */
    /*printf(".");*/
    if( (int)n==1 )
      return dlog_btil1( (int)s, a, b );
    else {
      double fprim1 = 2*(a-1.)-(b-1.)*n;
      if( fprim1 > 0 )  /* max at u=1, use derivative approximation */
        return( -log( n+1 )-log(1.+2*s/fprim1));
      else
	return( log_btiln( a+s, b+1, n )-log_btiln( a, b, n ));
    }
  }
}

/*------------------end of functions added by mmp--------------------------*/
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
void resample_c_t0(unsigned int idx)
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
		/*if(idx < 5)
		  printf( "tilb i=%d clust %d: ", idx, c );*/

		for (j = 0; j < t; ++j)
		  {
		    /* Add ratio of Beta functions here */
		    r_j_post = (double)g_S[index2d(j, g_max_t, clust)] + g_r_0[j] * g_nu_0;
		    N_j_post = (double)g_nc[index2d(j, g_max_t, clust)] + g_nu_0 + 1.0;
		    
		    /*if( 1 ) {  /* temporary for debugging */
		    if( j < g_t0 ) { /* do standard BETA_GIBBS */
		      g_c_prob[c] += log(N_j_post / (N_j_post + r_j_post + g_current_s[j]));
		      for (i = 0; i < g_current_s[j]; ++i)
			{
			  g_c_prob[c] += log((r_j_post + (double)i) / (r_j_post + N_j_post + (double)i));
			}
		    }
		    else {  /* compute ~Beta by numerical integration*/
		      double nj = (double)g_n - j;
		      /*if( idx < 5 ) printf("%.0f %.0f %.0f %.0f::", (double)g_current_s[j],r_j_post,N_j_post+1, nj);
		      double bsus = betatilde((double)g_current_s[j]+r_j_post,N_j_post+1, nj);
		      double bjos = betatilde( r_j_post,N_j_post, nj);*/
		      g_c_prob[c] += dlog_btilde( g_current_s[j], r_j_post, N_j_post+1, nj);
		    }
		  }
		/*if( idx < 5 ) printf("\n");*/
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

	/*
	if( idx < 5 ) {
	  printf( "Cluster probabilities for point %d:", idx );
	  for( c=0; c<g_num_valid_clusters; c++ )
	    printf("%f ", g_c_prob[ c ] ); 
	  printf("\n");
	  }*/

	/* Sample cluster assignment */
	c = sample_normalized_log(g_c_prob, NULL, g_num_valid_clusters + 1, g_inverseTemperature, g_unifrand[g_randIdx++]);
	g_c[idx] = (c < g_num_valid_clusters ? g_valid_clusters[c] : emptyC);

	/* Insert back into caches */
	modify_suff_stats(idx, true);

	/* Is this a new cluster? If so, then we need to re-estimate a centroid */
	if (g_c[idx] == emptyC)
	{
		resample_sigmarho_t0(emptyC);
	}
}
