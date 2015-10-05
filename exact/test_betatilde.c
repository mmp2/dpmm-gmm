#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <assert.h>

int g_ngridBeta;
double* g_grid01;
int isdebug = 0;

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
    //if( isdebug) printf( " %.3g ", bt );
  }
  bt += pow(n+1,-b+1)/2.;  /* includes extremities /2, trapeze method */
  bt *= g_grid01[ 1 ];     
  /* this is g_grid01[ 1 ]-g_grid01[ 0 ], since the second term is 0 */
  //if( isdebug ) printf( "%.3g\n", bt );
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
/** approximates log(btilde(a,b,n)) by gaussian around maximum
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
  //printf( "u=%3g ", u );
  //assert( u<1 ); /* otherwise, this function should not be called */
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
  double Hu = aloga(a)+aloga(b)-aloga(a+b); // Hu<0 
  //printf( "Hu=%.0f %.0f ", Hu, log( DBL_MIN )+5.0 );
  /* breakpoint currently just over a=b=500 */
  //  if((a<1) || (Hu > log( DBL_MIN )+5.0 ))  /* direct integration OK*/
  if( 0 )
    return( log( betatilde( a+s,b+1.,n) )-log( betatilde( a,b,n )));
  else { /* make approximations */
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

/* not used */
inline double fintegrand_btilde( double u, double a, double b, double n ) {
  /* to pass in b-1 */
  return( a*log( u )-(b)*(log(1.-pow(u,n+1))-log(1-u)) );
}

old_main( int argc, char* argv[] ) { //tests only betatilde
  
  if( argc == 1 ) {
    printf( "Usage\n  test_betatilde n grid01 amax bmax  \n amax > 1, bmax>= 3, recommended amax < 3 bmax\n" );
    exit(1);
  }
  int n = atoi( argv[ 1 ] );
  int ngrid01 = atoi( argv[ 2 ] );
  int iamax = atoi( argv[ 3 ] );
  int ibmax = atoi( argv[ 4 ] );

  int ia;
  int ib;
  double dn = (double)n;
  g_ngridBeta = ngrid01;
  int i;
  g_grid01 = (double *)malloc( ngrid01*sizeof( double ));
  for( i = 0; i < ngrid01+2; i++ ) 
    g_grid01[ i ] = (double)i/ngrid01;

  for( ib = 3; ib <= ibmax; ib++ ) {
    double b = (double)ib;
    for( ia = 1; ia <= iamax; ia++ ) {
      double a = (double)ia;
      double xx = betatilde( a, b, dn );
      printf( "%f\t", xx );
    }
    printf( "\n" );
  }
}

main( int argc, char* argv[] ) {
  
  if( argc == 1 ) {
    printf( "Usage\n  test_betatilde n grid01 a b s  \n a > 1, b>= 3, isdebug recommended a < 3 b\n" );
    exit(1);
  }
  int n = atoi( argv[ 1 ] );
  int ngrid01 = atoi( argv[ 2 ] );
  int ia = atoi( argv[ 3 ] );
  int ib = atoi( argv[ 4 ] );
  if( argc>5 )
    isdebug = atoi( argv[ 5 ] );

  double dn = (double)n;
  g_ngridBeta = ngrid01;
  int i;
  g_grid01 = (double *)malloc( ngrid01*sizeof( double ));
  for( i = 0; i < ngrid01+2; i++ ) 
    g_grid01[ i ] = (double)i/ngrid01;

  double a = (double)ia;
  double b = (double)ib;

  int is; 
  double z=0, zap = 0, bex = 0;
  for( is = 0; is <=n; is ++ ){
    double s = (double)is;
    double bexs = betatilde( a+s, b+1, dn );
    bex = betatilde( a, b, dn );
    double dlbet = dlog_btilde( s, a, b, dn );
    z += bexs;
    zap += exp( dlbet );

    printf( "%d   %f\t %g \n", is, log(bexs)-log(bex), dlbet );

  }
  printf( "z=%g bex=%g zap=%f (must have z = bex, zap = 1)\n", z, bex, zap );
}
