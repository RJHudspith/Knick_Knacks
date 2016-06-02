#include "fitfunc.h"
#include "resampled_ops.h"

// choose the solver you want
#define RK4ADAPTIVE

// defined in the on shell scheme ...
#define MC 1.6 //1.29 // charm mass  OS Scheme :: 1.6
#define MB 4.7 //4.19 // bottom mass OS Scheme :: 4.7
#define MT 172.9 // Top mass
#define MZ 91.1876 // z-mass
#define MS_TOL 1E-10

#ifndef RK4ADAPTIVE
#define acc 1E-3 // stepsize of our integrator 1E-2 or 1E-3 suffice for RK4
                 // 1E-3 or 1E-4 for LEAPFROG and ~1E-6 for EULER 
#endif

#define _4Pi    0.079577471545947673
#define _16Pi2  0.006332573977646111
#define _64Pi3  0.00050393022551874201
#define _256Pi4 4.0101493182360688e-05

// cash-karp is my favourite
#define CASHKARP


// all hail glorious beta function!
static inline double
beta_function( const double alpha , const int nf , const int loops )
{
  register const double nf_3 = (double)nf / 3.0 ;
  double beta0 = 0. , beta1 = 0. , beta2 = 0. , beta3 = 0. ;
  // try with gotos
  beta0 = ( 11. - 2. * nf_3 ) * _4Pi ;
  if( loops > 1 ) goto loop2 ;
  return -alpha * alpha * beta0 ;
 loop2 :
  beta1 = ( 102.0 - 38. * nf_3 ) * _16Pi2 ;
  if( loops > 2 ) goto loop3 ;
  return -alpha * alpha * ( beta0 + alpha * beta1 ) ;
 loop3 : 
  beta2 = ( 1428.5 + ( nf_3 / 3.0 ) * ( -2516.5 + nf_3 * 162.5 ) ) * _64Pi3 ;
  if( loops > 3 ) goto loop4 ;
  return -alpha * alpha * ( beta0 + alpha * ( beta1 + alpha * beta2 ) ) ;
 loop4 :
  beta3 = ( 29242.96414 + nf * ( -6946.289617 + nf * ( +405.0890405 + nf * 1.499314129 ) ) ) * _256Pi4 ;
  return -alpha * alpha * ( beta0 + alpha * ( beta1 + alpha * ( beta2 + alpha * beta3 ) ) ) ;
}

/********************************** Integration Schemes live here ***********************************/

// run using the RK4
#define one_6 0.16666666666666666
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
inline double adaptfmax( const double a , const double b ) { return ( b < a ? a : b ) ; }
inline double adaptfmin( const double a , const double b ) { return ( a < b ? a : b ) ; }

/*
  These are the embedded cash-karp factors, their actual values
  can be found in NRC or wikipedia or something.
 */
const double b21 = 0.2 ;
const double b31 = 0.075 ;
const double b32 = 0.225 ;
const double b41 = 0.3 ;
const double b42 = -0.9 ;
const double b43 = 1.2 ;
const double b51 = -0.20370370370370369 ;
const double b52 = 2.5 ;
const double b53 = -2.5925925925925926 ;
const double b54 = 1.2962962962962963 ;
const double b61 = 0.029495804398148147 ;
const double b62 = 0.341796875 ;
const double b63 = 0.041594328703703706 ;
const double b64 = 0.40034541377314814 ;
const double b65 = 0.061767578125 ;
const double c1 = 0.097883597883597878 ;
const double c3 = 0.40257648953301128 ;
const double c4 = 0.21043771043771045 ;
const double c6 = 0.28910220214568039 ;
const double dc1 = -0.0042937748015873106 ;
const double dc3 = 0.018668586093857853 ;
const double dc4 = -0.034155026830808066 ;
const double dc5 = -0.019321986607142856 ;
const double dc6 = 0.039102202145680387 ;
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

/***************************************** Running codes live here ********************************/

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
  register const double murat = muprime/mu ;
  const double diff = 2.0 * log( murat ) ; 
  mu = RK4_step( mu , &alpha_mup , nf , loops , diff , murat ) ;

  return alpha_mup ;
}



/****************************************************************************

  Matching coefficients come from the chetyrkin, steinhauser and kniel paper

 ****************************************************************************/

// match a nf = N-1 to a nf = N theory
enum condition { NOT_MET , MET } ;
static inline double
match_up( const double Alpha_Ms , const int nf , const int loops )
{
  if( loops <= 2 )
    return Alpha_Ms ; // nice result, threshold effects only kick in at two loops
                      // but we always match the running to one lower order in the matching ...

  // have to use a Newton-raphson .. 
  const double h2_term = 0. ;
  const double h3_term = 0.152777777778 ;
  const double h4_term = 0.97205669 - 0.084651492 * ( nf - 1 ) ;

  //Loop size of Alpha
  const double hmgg = Alpha_Ms / M_PI ; 
  
  //set up guess and new_guess
  double guess = hmgg ;
  
  // don't need to do so many as halley's method converges like x^3
  int CONV_CRITERIA = NOT_MET , counter = 0 ;
  while( CONV_CRITERIA != MET ) {
    //function and its first and second derivatives
    const double f_x = guess * ( 1 + guess * ( h2_term + guess * ( h3_term + guess * ( h4_term ) ) ) ) - hmgg ; 
    const double f_prime = 1 + guess * ( 2 * h2_term + guess * ( 3 * h3_term + guess * ( 4 * h4_term ) ) ) ; 
    const double f_pprime = 2 * h2_term + guess * ( 6 * h3_term + guess * ( 12 * h4_term ) ) ;

    //halley's method	  
    const double update = guess - 2 * f_x * f_prime / (2 * f_prime * f_prime - f_x * f_pprime) ;

    if( ( fabs( update - guess ) < 1E-14 ) ) {
      CONV_CRITERIA = MET ;
    } else if( counter > 50 ) {
      #ifdef DEBUG_RUNNINGS
      printf( "Ill convergent series here (Match up) ! alpha :: %1.15f" , Alpha_Ms ) ;
      #endif
      CONV_CRITERIA = MET ;
    } else {
      guess = update ;
      counter ++ ;
    }
  }
  return  M_PI * guess ;
} 

// I invert this numerically opposed to using the inverse series so I have numerical forwards
//  -> backwards symmetry
static inline double
match_up_OS( const double Alpha_Ms , const double Mh , const int nf , const int loops )
{
  if( loops <= 2 )
    return Alpha_Ms ; // nice result, threshold effects only kick in at two loops
                      // but we always match the running to one lower order in the matching ...
  double h2_term = 0. , h3_term = 0. , h4_term = 0. ;
  // switch for the scheme we are matching to
  if( loops > 2 ) { h3_term = -0.02955201190 ; }
  if( loops > 3 ) { h4_term = -0.1717036285 + ( nf - 1 ) * 0.008465086429 ; }

  //Loop size of Alpha
  const double hmgg = Alpha_Ms ; 
  
  //set up guess and new_guess
  double guess = hmgg ;
  
  // don't need to do so many as halley's method converges like x^3
  int CONV_CRITERIA = NOT_MET , counter = 0 ;
  while( CONV_CRITERIA != MET ) {
    //function and its first and second derivatives
    const double f_x = guess * ( 1 + guess * ( h2_term + guess * ( h3_term + guess * ( h4_term ) ) ) ) - hmgg ; 
    const double f_prime = 1 + guess * ( 2 * h2_term + guess * ( 3 * h3_term + guess * ( 4 * h4_term ) ) ) ; 
    const double f_pprime = 2 * h2_term + guess * ( 6 * h3_term + guess * ( 12 * h4_term ) ) ;
    //halley's method	  
    const double update = guess - 2 * f_x * f_prime / (2 * f_prime * f_prime - f_x * f_pprime) ;

    if( fabs( update - guess ) < MS_TOL ) {
      CONV_CRITERIA = MET ;
    } else if ( counter > 50 ) {
      #ifdef DEBUG_RUNNINGS
      printf( "Ill convergent series here (Match up OS) ! alpha :: %1.15f" , Alpha_Ms ) ;
      #endif
      CONV_CRITERIA = MET ;
    } else {
      guess = update ;
      counter ++ ;
    }	
  }
  return ( guess) ;
} 

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

// run our lattice sim to the Z mass
static double
run_to_mz( const double Alpha_Ms , const double scale , 
	   const int nf , const int loops , const bool threshold )
{
  double alpha = Alpha_Ms ;
  if( threshold ) {
    // with matching ?
    /*
      If I understand this correctly, I have to run to the \mu scale mass in MS and then up
      to the On-shell mass, and this is where I can kick in the higher nf running
      It is only effective at the three loop and above level, because you use one fewer
      loop in the matching
    */
    alpha = RUN( scale , alpha , MC , 3 , loops ) ;
    alpha = match_up_OS( alpha , MC , 4 , loops ) ;
    alpha = RUN( MC , alpha , MB , 4 , loops ) ;
    alpha = match_up_OS( alpha , MB , 5 , loops ) ;
    alpha = RUN( MB , alpha , MZ , 5 , loops ) ;
  } else {
    alpha = RUN( scale , alpha , MC , 3 , loops ) ;
    alpha = RUN( MC , alpha , MB , 4 , loops ) ;
    alpha = RUN( MB , alpha , MZ , 5 , loops ) ;
  }
  return alpha ;
}

// Wraps the measurements back into Alpha_Ms at Mz , sped up for Jamie's sanity
struct resampled
boot_run_MS_quick( const double mu ,
		   const struct resampled Alpha_Ms ,
		   const int loops ,
		   const int nf ,
		   const bool threshold ) 
{
  struct resampled MZ_couple ;
  MZ_couple.resampled = malloc( Alpha_Ms.NSAMPLES * sizeof( double ) ) ;

  // equate
  equate( &MZ_couple , Alpha_Ms ) ;

  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < MZ_couple.NSAMPLES ; i++ ) {

    MZ_couple.resampled[i] = run_to_mz( MZ_couple.resampled[i] , 
					mu , nf , 
					loops , threshold ) ;
  }

  MZ_couple.avg = run_to_mz( MZ_couple.avg , 
			     mu , nf , 
			     loops , threshold ) ;
  compute_err( &MZ_couple ) ;

  return MZ_couple ;
}

// defined in the on shell scheme ...
#undef MC
#undef MB
#undef MT
#undef MZ
#undef MS_TOL

