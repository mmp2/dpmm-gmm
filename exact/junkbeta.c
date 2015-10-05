#include "mex.h"
#include "matrix.h"
#include <limits.h>
#include <float.h>

#include "mex.h"
#include "matrix.h"
#include <limits.h>
#include <float.h>


/* todo's: call only when a large, b small
           inline?
   This function has numerical errors when a,b too large. Now it is called from
   dlog_btilde which calculates the difference in log Btilde()'s 
   avoiding numerical errors.
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
    bt += pow(x,a-1)*(pow(1-x,b-1)/pow((1-pow(x,n+1)),b-1));
  }
  bt += pow(n+1,-b+1)/2.;  /* includes extremities /2, trapeze method */
  bt *= g_grid01[ 1 ];     
  /* this is g_grid01[ 1 ]-g_grid01[ 0 ], since the second term is 0 */
  return( bt );
}
/* Calculates 
 log( betatilde( a+s,b+1,n)/betatilde( a,b,n) ) in a (hopefully numerically
 stable way.

 notes: 1. the approximation by Beta can be replaced with a gradient descent,
         to be implemented
  2. other refinements possible
  3. speedups/inlining possible
*/

double dlog_btilde( double s, double a, double b, double n ) {
  double u0 = (s+a)/(s+a+b);
  double Hu = a*log(a)+b*log(b);
  if( (a+b+s)*Hu < -log( MINDOUBLE )-5.0 )  /* direct integration OK*/
    return( log( betatilde( a+s,b+1.,n) )-log( betatilde( a,b,n )));
  else { /* make approximations */
    if( n==1 )
      return dlog_btil1( iround(s), a, b );
    else if ( n==2 )
      return log_btil2( a+s, b+1 )-log_btil2( a, b );
    else { /* approximation by Beta */
      int i;
      double dl = 0;
      for (i = 0; i < s; ++i)
	 dl += log((a + (double)i)/(a + b + (double)i));
      return dl;
    }
  }
}
inline double fintegrand_btilde( double u, double a, double b, double n ) {
  /* to pass in b-1 */
  return( a*log( u )-(b)*(log(1.-pow(u,n+1))-log(1-u)) );
}

inline double aloga( double a ) { return a*log(a); }

/* Computes 
 log( betatilde( a+s,b+1,1)/betatilde( a,b,1 ) for a,b large by Laplace's 
 method.
 REQUIRES b>a>0
 POSTCONDITION: return value is always <=0 for correct inputs
*/
inline double dlog_btil1( double* value, int s, double a, double b ) {
  double b1 = b-1.;
  double d_once = aloga(b1)-aloga(b)+log(1./a+1./b1)*0.5;
  return( s == 0 ? (d_once+-aloga(b1-a)+aloga(b-a)-log(1./a+1./b)*0.5) : (d_once + aloga( a+1. )-aloga(a) -log(1./(a+1.)+1./b))*0.5));
}

/** approximates log(btilde(a,b,n)) by gaussian around maximum
 */
double log_btiln( double a, double b, double n ) {
  double b1 = b-1.;
  int i;
  /* find zero of derivative of f where
     f=a*ln(u)-(b-1)[ln(1-u^(n+1))-ln(1-u)] */
  int niter = ceil( log( 1e8 )/log( 2 ));  /* error is halfed each step so number of iterations is ~log_2(tolerance) */
  double umin = 0, umax = 1, u = 0.5;
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

  /* -f''(u) */
  double neg_fsecund = a/u/u+b*(1./(1.-u)/(1.-u)-(n+1.)*(n*un/u+un*un)/(1-un*u)/(1-un*u));
  
  /* f(u) */
  double f = a*log(u)-(b1)*(log(1-un*u)-log(1-u));
  /* calculate log(Betatilde) */
  return( f-log(neg_fsecund/2./3.1416)*0.5 );
}
 
/*------------------end of functions added by mmp--------------------------*/
