/**
   @file cruel_runnings.c
   @brief adaptive cash-carp running of the QCD strong coupling
 */
#include <stdio.h>
#include <math.h>

// choose the solver you want
#define RK4ADAPTIVE

// defined in the on shell scheme ...
#define MC 1.6 //1.29 // charm mass  OS Scheme :: 1.6
#define MB 4.7 //4.19 // bottom mass OS Scheme :: 4.7
#define MT 172.9 // Top mass
#define MZ 91.1876 // z-mass
#define HALLEY_TOL 1.E-14 // tolerance for Halley's method

// some common definitions
#define _4Pi (0.079577471545947673)
#define _16Pi2 (0.006332573977646111)
#define _64Pi3 (0.00050393022551874201)
#define _256Pi4 (4.0101493182360688e-05)

// These are the embedded cash-karp factors, their actual values
// can be found in NRC or wikipedia or something.
static const double b21 = 0.2 ;
static const double b31 = 0.075 ;
static const double b32 = 0.225 ;
static const double b41 = 0.3 ;
static const double b42 = -0.9 ;
static const double b43 = 1.2 ;
static const double b51 = -0.20370370370370369 ;
static const double b52 = 2.5 ;
static const double b53 = -2.5925925925925926 ;
static const double b54 = 1.2962962962962963 ;
static const double b61 = 0.029495804398148147 ;
static const double b62 = 0.341796875 ;
static const double b63 = 0.041594328703703706 ;
static const double b64 = 0.40034541377314814 ;
static const double b65 = 0.061767578125 ;
static const double c1 = 0.097883597883597878 ;
static const double c3 = 0.40257648953301128 ;
static const double c4 = 0.21043771043771045 ;
static const double c6 = 0.28910220214568039 ;
static const double dc1 = -0.0042937748015873106 ;
static const double dc3 = 0.018668586093857853 ;
static const double dc4 = -0.034155026830808066 ;
static const double dc5 = -0.019321986607142856 ;
static const double dc6 = 0.039102202145680387 ;

// OK cool, write the adaptive bit
const double ADAPTIVE_EPS = 1E-12 ;
// Standard shrink and factor from NRC
const double ADAPTIVE_SHRINK = -0.25 ;
// Standard growth and factor from NRC
const double ADAPTIVE_GROWTH = -0.20 ;
// define adaptive safe
const double ADAPTIVE_SAFE = 0.9 ;
// adaptive error conserving
#define ADAPTIVE_ERRCON pow( 5./ADAPTIVE_SAFE , 1./ADAPTIVE_GROWTH )

// fmin and fmax functions for the adaptive
static inline double adaptfmax( const double a , const double b ) { return ( b < a ? a : b ) ; }
static inline double adaptfmin( const double a , const double b ) { return ( a < b ? a : b ) ; }

// cash-karp is my favourite
#define CASHKARP

static double beta0 = 0.0 , beta1 = 0.0 , beta2 = 0.0 , beta3 = 0.0 ;

// run using the RK4
#define one_6 (0.16666666666666666)

// need to call this each time we pass through a threshold
static void init_betas( const int nf )
{
  register const double nf_3 = nf / 3.0 ;
  beta0 = ( 11. - 2. * nf_3 ) * _4Pi ;
  beta1 = ( 102.0 - 38. * nf_3 ) * _16Pi2 ;
  beta2 = ( 1428.5 + ( nf_3 / 3.0 ) * ( -2516.5 + nf_3 * 162.5 ) ) * _64Pi3 ;
  beta3 = ( 29242.96414 + nf * ( -6946.289617 + nf * ( +405.0890405 + nf * 1.499314129 ) ) ) * _256Pi4 ;
}

// all hail glorious beta function!
static inline double
beta_function( const double alpha , const int nf , const int loops )
{
  // try with gotos
  if( loops > 1 ) goto loop2 ;
  return -alpha * alpha * beta0 ;
 loop2 :
  if( loops > 2 ) goto loop3 ;
  return -alpha * alpha * ( beta0 + alpha * beta1 ) ;
 loop3 : 
  if( loops > 3 ) goto loop4 ;
  return -alpha * alpha * ( beta0 + alpha * ( beta1 + alpha * beta2 ) ) ;
 loop4 :
  return -alpha * alpha * ( beta0 + alpha * ( beta1 + alpha * ( beta2 + alpha * beta3 ) ) ) ;
}

///////////////// Integration Schemes live here 
static double
RK4_step( const double mu , double *__restrict alpha , 
	  const int nf , const int loops ,
	  const double dlq , const double exp_step )
{
  const double F1 = beta_function( *alpha , nf , loops ) ;
  const double F2 = beta_function( *alpha + 0.5 * dlq * F1 , nf , loops ) ;
  const double F3 = beta_function( *alpha + 0.5 * dlq * F2 , nf , loops ) ;
  const double F4 = beta_function( *alpha + dlq * F3 , nf , loops ) ; 
  *alpha += dlq * ( F1 + 2.0 * F2 + 2.0 * F3 + F4 ) * one_6 ;
  return mu * exp_step ;
}

// this is the cash-karp embedded adaptive runge-kutta step
static void
RK4_adapt( const double mu , double *__restrict alpha , 
	   const int nf , const int loops ,
	   const double dlq , 
	   double *__restrict err )
{
  const double K1 = beta_function( *alpha , nf , loops ) ;
  const double K2 = beta_function( *alpha + dlq*b21*K1 , nf , loops ) ;
  const double K3 = beta_function( *alpha + dlq*(b31*K1 + b32*K2) , nf , loops ) ;
  const double K4 = beta_function( *alpha + dlq*(b41*K1 + b42*K2 + b43*K3) , nf , loops ) ;
  const double K5 = beta_function( *alpha + dlq*(b51*K1 + b52*K2 + b53*K3 + b54*K4) , nf , loops ) ;
  const double K6 = beta_function( *alpha + dlq*(b61*K1 + b62*K2 + b63*K3 + b64*K4 + b65*K5) , nf , loops ) ;
  // compute the new alpha
  *alpha += dlq*( c1*K1 + c3*K3 + c4*K4 + c6*K6 ) ; 
  *err = dlq*( dc1*K1 + dc3*K3 + dc4*K4 + dc5*K5 + dc6*K6 ) ;
  return ;
}

// driver for the adaptive algorithm
static double
adaptive( const double mu , double *__restrict alpha ,
	  const int nf , const int loops , 
	  double *__restrict dlq )
{
  double errmax = 10.0 , a1 ;
  while( errmax > 1.0 ) {
    // initialise our attempts
    a1 = *alpha ;
    RK4_adapt( mu , &a1 , nf , loops , *dlq , &errmax ) ;

    errmax = fabs(errmax) / ADAPTIVE_EPS ;
    if( errmax < 1.0 ) { break ; } 

    // set up a tolerance s.t del_temp doesn't go too crazy, only really used when starting guess is bad
    register const double del_temp =  ADAPTIVE_SAFE * (*dlq) * pow( errmax , ADAPTIVE_SHRINK ) ; 
    register const double tol = 0.1 * (*dlq) ; 
    *dlq = ( 0.0 < del_temp ? adaptfmax( del_temp , tol ) : adaptfmin( del_temp , tol ) ) ;
  }
  // and set the variables to their new values
  *alpha = a1 ;
  const double step = mu * exp( 0.5 * (*dlq) ) ;
  *dlq = errmax > ADAPTIVE_ERRCON ? ADAPTIVE_SAFE * (*dlq) * pow( errmax , ADAPTIVE_GROWTH ) : ADAPTIVE_SAFE * 5.0 * (*dlq) ;
  return step ;
}

//////////////////////////////////// Running codes live here 

// run from mu to mu-prime
static double
RUN( double mu , const double alpha_mu , const double muprime , 
     const int nf , const int loops )
{
  // these get initialised if we need them, because the integration is constant in
  // the step size "mu" we need only compute the exponentials here.
  // However, for the adaptive this is not the case.
  double alpha_mup = alpha_mu ;
  int flag = 0 ;

  // perfectly sensible starting point
  double dlq = 0.01 ;
  while( mu < muprime ) {
    mu = adaptive( mu , &alpha_mup , nf , loops , &dlq ) ;
    flag = 1 ;
  }
  // run backwards ...
  double mdlq = -0.01 ;
  while( mu > muprime && flag == 0 ) {
    mu = adaptive( mu , &alpha_mup , nf , loops , &mdlq ) ;
  }

  // final RK4 step to go to exactly the scale we wish to be
  register const double murat = muprime/mu ;
  const double diff = 2.0 * log( murat ) ; 
  mu = RK4_step( mu , &alpha_mup , nf , loops , diff , murat ) ;

  return alpha_mup ;
}

// I invert this numerically opposed to using the inverse series so I have numerical forwards
//  -> backwards symmetry
static double
match_up_OS( const double Alpha_Ms , const double Mh , const int nf , const int loops )
{
  if( loops <= 2 )
    return Alpha_Ms ; // nice result, threshold effects only kick in at two loops
                      // but we always match the running to one lower order in the matching ...
  double h2_term = 0. , h3_term = 0. , h4_term = 0. ;
  // switch for the scheme we are matching to
  if( loops > 2 ) { h3_term = -0.02955201190 ; }
  if( loops > 3 ) { h4_term = -0.1717036285 + ( nf - 1 ) * 0.008465086429 ; }

  //set up guess and new_guess
  double guess = Alpha_Ms , correction = 1.0 ;
  
  // don't need to do so many as halley's method converges like x^3
  int counter = 0 ;
  while( fabs( correction ) > HALLEY_TOL ) {
    //function and its first and second derivatives
    const double f_x = guess * ( 1. + guess * ( h2_term + guess * ( h3_term + guess * ( h4_term ) ) ) ) - Alpha_Ms ; 
    const double f_prime = 1. + guess * ( 2. * h2_term + guess * ( 3. * h3_term + guess * ( 4. * h4_term ) ) ) ; 
    const double f_pprime = 2. * h2_term + guess * ( 6. * h3_term + guess * ( 12. * h4_term ) ) ;
    //halley's method	  
    correction = -2. * f_x * f_prime / (2. * f_prime * f_prime - f_x * f_pprime) ;

    //printf( "CORRECTION %d %e \n" , counter , correction ) ;

    guess += correction ;

    // make sure it doesn't continue forever
    if( counter > 50 ) {
      printf( "Match up OS not converging ... Leaving \n" ) ;
      return -1 ;
    }
    counter ++ ;
  }
  return guess ;
} 

// matching down from nf to nf-1 flavours at scale Mh
static inline double
match_down_OS( const double Alpha_Ms , const double Mh , const int nf , const int loops )
{
  if( loops <= 2 )
    return Alpha_Ms ; // nice result, threshold effects only kick in at two loops
                      // but we always match the running to one lower order in the matching ...
  double h3_term , h4_term ;
  if( loops > 2 ) goto loop3 ;
  return Alpha_Ms ;
 loop3 :
  h3_term = -0.02955201190 ;
  if( loops > 3 ) goto loop4 ;
  return Alpha_Ms * ( 1.0 + Alpha_Ms * Alpha_Ms * h3_term ) ; 
 loop4 : 
  h4_term = -0.1717036285 + ( nf - 1 ) * 0.008465086429 ;
  return Alpha_Ms * ( 1.0 + Alpha_Ms * Alpha_Ms * ( h3_term + Alpha_Ms * h4_term ) ) ; 
} 

// Direct computation using the higher order terms
static double
lambda_MS( const double alpha , const double mu , const int nf , const int loops )
{
  init_betas( nf ) ;

  const double b1 = beta1 / beta0 ;
  const double b2 = beta2 / beta0 ;
  const double b3 = beta3 / beta0 ;

  // constant term ...
  double expr = b1 * log( beta0 ) ;

  if( loops > 1 ) { expr += ( 1.0 / alpha + b1 * log( alpha ) ) ; goto loop2 ; }
 loop2 :
  if( loops > 2 ) { expr += ( b1 - b2 * b2 ) * alpha ; goto loop3 ; }
 loop3 :
  if( loops > 3 ) { expr += ( b3 / 2.0 - b1 * b2 + b1 * b1 * b1 / 3.0 ) * alpha * alpha ; }

  return mu * exp( -expr / ( 2.0 * beta0 ) ) ;
}

// So, the point of this one is to run to two very large mu's
// and fit to 1/(their exponent) to get the asymptote at mu -> infinity
static double
lambda_MS_asymp( const double alpha , const double mu , const int nf , const int loops )
{
  init_betas( nf ) ;

  const double V1 = lambda_MS( RUN( mu , alpha , 1E10 , nf , loops ) , 1E10 , nf , loops ) ;
  const double V2 = lambda_MS( RUN( mu , alpha , 1E100 , nf , loops ) , 1E100 , nf , loops ) ;

  return ( V2 - ( V2 - V1 ) / ( 0.01 - 0.1 ) * 0.01 ) ;
}


// newton raphson ( actually Halley's method) solve for Lambda
// is gloriously fast, but doesn't give a particularly stable Lambda
// uses the 4 loop expression
static double
NR_lambda( const double alpha , const double mu , const int nf )
{
  // initialise number of flavours
  init_betas( nf ) ;

  // solve equation for L which is log( \mu^2 / \Lambda^2 )
  double guess = log( mu * mu / ( 0.2 * 0.2 ) ) , correction = 1.0 ;

  // uses chetyrkin's notation
  const double b0 = M_PI * beta0 ;
  const double b1 = M_PI * M_PI * beta1 / b0 ;
  const double b2 = M_PI * M_PI * M_PI * beta2 / b0 ;
  const double b3 = M_PI * M_PI * M_PI * M_PI * beta3 / b0 ;
  const double a = alpha / M_PI ;

  int counter = 0 ;
  while( fabs( correction ) > HALLEY_TOL ) {

    // cache in some stuff
    register const double b0L = b0 * guess ;
    register const double logL = log( guess ) ;

    // compute f, f' and f'' for Halley's method
    const double f = ( b3/2. -b1 * ( 3. * b2 * logL + b1 * b1 * ( 1.0 + logL * ( -4. + logL * ( -5. + logL * 2.0 ) ) ) / 2.0 )
		       +b0L * ( -(b2 + b1*b1*( -1 + logL * ( -1.0 + logL ))) + b0L * ( b1*logL + b0L * ( -1.0 + b0L * a ) ) ) ) ;

    const double fp = ( b1 * ( 2.*b1*b1 - 3.*b2 )
		       + b0L * ( 2.0 * b1 * b1 - b2 + b0L * ( b1 + b0L * ( -3. + 4.0 * b0L * a) ) ) 
		       + b1*logL*( 5.*b1*b1 + b0L * ( -b1 + 2. * b0L ) - b1*logL * (3*b1 + b0L) ) ) / guess ;

    const double fpp = ( +3.0 * b1 * ( b1*b1 * ( 1.0 + logL*logL ) + b2 ) 
			 +b0L * ( -b1*b1*( 1.0 + 2.0*logL ) 
				  + b0L * ( b1 * ( 3.0 + 2.0 * logL ) 
					    + b0L * ( -6.0 + 12. * b0L * a ) ) ) ) / ( guess * guess ) ;

    correction = -f*fp / (  fp*fp  - 0.5 * fpp*f ) ;

    guess += correction ;

    if( counter > 50 ) {
      printf( "Lambda not converging ... Leaving \n" ) ;
      return -1 ;
    }
    counter++ ;
  }

  return mu * exp( -guess * 0.5 ) ;
}

// running from nf=3 to MZ matching at the charm
// I want to automate this for many nfs
static double 
run_nf3_2MZ( const double alpha , const double mu , const int nf , const int loops )
{
  init_betas( 3 ) ;
  double alpha_mup = RUN( mu , alpha , MC , 3 , loops ) ;
  alpha_mup = match_up_OS( alpha_mup , MC , 4 , loops ) ;
  init_betas( 4 ) ;
  alpha_mup = RUN( MC , alpha_mup , MB , 4 , loops ) ;
  alpha_mup = match_up_OS( alpha_mup , MB , 5 , loops ) ;
  init_betas( 5 ) ;
  return RUN( MB , alpha_mup , MZ , 5 , loops ) ;
}

// run nf=3 to MZ
static double 
run_nf4_2MZ( const double alpha , const double mu , const int nf , const int loops )
{
  init_betas( 4 ) ;
  double alpha_mup = RUN( mu , alpha , MB , 4 , loops ) ;
  alpha_mup = match_up_OS( alpha_mup , MB , 5 , loops ) ;
  init_betas( 5 ) ;
  return RUN( MB , alpha_mup , MZ , 5 , loops ) ;
}

// running from nf=5 at MZ matching at the charm and back to some value in Nf=3
// I want to automate this for many nfs
static double 
run_MZ_2nf3( const double alpha , const double mu , const int nf , const int loops )
{
  init_betas( 5 ) ;
  double alpha_mup = RUN( MZ , alpha , MB , 5 , loops ) ;
  alpha_mup = match_down_OS( alpha_mup , MB , 5 , loops ) ;
  init_betas( 4 ) ;
  alpha_mup = RUN( MB , alpha_mup , MC , 4 , loops ) ;
  alpha_mup = match_down_OS( alpha_mup , MC , 4 , loops ) ;
  init_betas( 3 ) ;
  return RUN( MC , alpha_mup , mu , 3 , loops ) ;
}

// a little example to get us going
int main( void ) 
{
  const int loops = 5 ;
  const int nf = 3 ;
  const double mu = 2.0 ;

  printf( "%f \n" , run_nf3_2MZ( 0.2966 , mu , nf , loops ) ) ;

  return 0 ;
}
