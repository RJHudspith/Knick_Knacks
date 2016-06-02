/**
   A fit to the OPE for the Vacuum Polarisation Function to obtain alpha_s

   The form I use here is ::

   (form)

   so \alpha_s / Pi is in params[0]
 */

#include "fitfunc.h"
#include "Utils.h"
#include "coefficients.h"

// renormalisation scale
static int LOOP_ORDER = 5 ;

// diffs
static double q2[ 10 ] ;
static bool q2q3_set[ 10 ] ;

//#define NONMUL

// number of fit params
int D0_diff_multi_n( void )
{
  return 1 ;
}

//
void set_q2q3_D0_diff_multi( const double Q2 , 
			     const int NSIM ) 
{
  q2[ NSIM ] = sqrt( Q2 ) ; 
  
  q2q3_set[ NSIM ] = true ;

  printf( "Q2[%d] SET !! %f \n" , NSIM , q2[NSIM] ) ; 

  return ;
}

// set the constant mu
void set_loops_D0_diff_multi( const int loopsnew )
{
  LOOP_ORDER = loopsnew ;
  printf( "NUMBER OF LOOPS SET !!! %d \n" , LOOP_ORDER ) ;
  return ;
}

// provide a description of the function
void
D0_diff_multi_description( const char *message ,
		const double *__restrict params ,
		const struct x_descriptor X ,
		const int NPARAMS )
{
  printf( "ALPHA_GUESS :: %e + ( a^2 q^2 ) %f  \n" , 
	  params[0] , params[1] ) ;
  return ;
}

// I provide a guess too
void
D0_diff_multi_guess( double *__restrict params ,
	       void *data ,
	       const int NPARAMS ) // nparams is already in data no?
{
  const struct data* DATA = (const struct data*)data ;

  params[0] = 0.3 ; 
  int k , xposit = 0 ;

  // everything else we set to 0
  xposit = 0 ;
  for( k = 0 ; k < DATA->SIMS ; k++ ) {
    double p[ DATA-> NPARAMS ] ;
    p[0] = params[0] ;

    int j ;
    for( j = 1 ; j < DATA->NPARAMS ; j++ ) {
      p[ j ] = 0.0 ;
    }
    printf( "TEST2 :: %f %f \n" , DATA->quarks[0].ainverse , DATA->quarks[1].ainverse ) ;
    const struct mom_info dummy ;
    const struct x_descriptor XAXIS = { 0.0 , DATA->quarks[k] , dummy , 0 } ;
    #ifndef QUIET
    D0_diff_multi_description( "GUESS" , p , XAXIS , DATA-> NPARAMS ) ;
    #endif
    xposit += DATA->NDATA[k] ;
  }
  return ;
}

static inline double
num( const double *fitparams , const double XX , const double XXref , const double a2 , const double mu , const int der )
{
  const double a = fitparams[0] / M_PI ;
  return ( D0_OPE( log( ( XX ) / ( mu * mu ) ) , a ,  der , LOOP_ORDER )
	   - D0_OPE( log( ( XXref ) / ( mu * mu ) ) , a , der , LOOP_ORDER ) ) ;
}

// evaluate the function at the point "X" using the fitparams
inline double
D0_diff_multi_eval( const double *__restrict x ,
		    const struct x_descriptor X ,
		    const int NPARAMS )
{
  printf( "IN DIS :: %f \n" , X.quark.mu ) ;
  const double a2 = 1.0 / ( X.quark.ainverse * X.quark.ainverse ) ;
  return num( x , X.X * X.X , q2[X.LT] * q2[X.LT] , a2 , X.quark.mu , 0 ) ;
}

// evaluate the function and its difference from the data
int
D0_diff_multi_f( const gsl_vector *__restrict x , 
		 void *data , 
		 gsl_vector *__restrict f )
{
  const struct data* DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  int NSIM , count = 0 , xposit = 0 ;
  for( NSIM = 0 ; NSIM < DATA->SIMS ; NSIM++ ) {

    if( q2q3_set[ NSIM ] == false ) {
      printf( "Q2 and Q3 are not set!!! exiting \n" ) ;
      exit(1) ;
    }

    // local idx gives the unshared fit parameters
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON)+DATA->NCOMMON ;

    const double ai = ( DATA->quarks[NSIM].ainverse ) ;
    const int lopos = find_idx( DATA->FIT_LO*ai , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI*ai , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;

    printf( "WHAT %d %d %f %d \n" , lopos , hipos , ai , DATA->NDATA[NSIM] ) ;

    size_t i  ;
    for( i = lopos ; i <= hipos ; i++ ) {

      const struct x_descriptor XAXIS = { DATA->X[i] , DATA->quarks[NSIM] , DATA->mom[i] , NSIM } ;

      // local parameters ...
      double p[ DATA->NPARAMS ] ;

      p[0] = params[loc_idx] ;
      //p[1] = params[loc_idx] ;
      //p[2] = params[loc_idx+1] ;

      gsl_vector_set ( f , count , 
		       ( D0_diff_multi_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;

      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  free( (void*)params ) ;
  return GSL_SUCCESS;
}

// and the derivatives
int
D0_diff_multi_df( const gsl_vector *__restrict x, 
	    void *data, 
	    gsl_matrix *__restrict J )
{
  // have to check this
  if( q2q3_set == false ) {
    printf( "Q2 and Q3 are not set!!! exiting \n" ) ;
    exit(1) ;
  }

  const struct data* DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  int count = 0 , NSIM , xposit = 0 ; 
  for( NSIM = 0 ; NSIM < DATA-> SIMS ; NSIM++ ) {

    // make sure these are set
    if( q2q3_set[ NSIM ] == false ) {
      printf( "Q2 and Q3 are not set!!! exiting \n" ) ;
      exit(1) ;
    }

    // lattice spacings
    const double a = 1.0 / ( DATA->quarks[NSIM].ainverse ) ;
    const double a2 = a * a ;

    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON)+DATA->NCOMMON ;
    const int lopos = find_idx( DATA->FIT_LO/a , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI/a , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;

    size_t i ; 
    for( i = lopos ; i <= hipos ; i++ ) {

      const double _s = 1.0 / DATA->sigma[i] ;
 
      const double p2 = DATA->mom[i].p2 ;

      // initialise whole column to zero
      int k ;
      for( k = 0 ; k < DATA->LOGICAL_NPARS ; k++ ) {
	gsl_matrix_set( J , count , k ,  0.0 ) ;
      }

      // local parameters ...
      double p[ DATA->NPARAMS ] ;

      int j ;
      for( j = 0 ; j < DATA -> NPARAMS ; j++ ) {
	p[j] = params[j] ;
      }

      const double p2ref = q2[ NSIM ] * q2[ NSIM ] ;

      gsl_matrix_set ( J , count , loc_idx , _s * num( p , p2 , p2ref , a2 , DATA->quarks[NSIM].mu , 1 ) ) ;

      count ++ ;
    }
    xposit += DATA->NDATA[NSIM] ; 
  }
  free( (void*)params ) ;

  return GSL_SUCCESS;
}

int
D0_diff_multi_fdf (const gsl_vector *__restrict x , 
	     void *data ,
	     gsl_vector *__restrict f , 
	     gsl_matrix *__restrict J )
{
  D0_diff_multi_f( x , data , f ) ;
  D0_diff_multi_df( x , data , J ) ;
  return GSL_SUCCESS ;
}

 
