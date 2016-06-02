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
static double mu = 2.0 ;
static int LOOP_ORDER = 5 ;

// diffs
static double q2[ 10 ] ;
static double t2[ 10 ] ;
static bool q2q3_set[ 10 ] ;

// compute this
const static double pisq4 = 4.0 * M_PI * M_PI ;

//#define NONMUL

// number of fit params
int D0_diff_v1_n( void )
{
  return 2 ;
}

//
void set_q2q3_D0_diff_v1( const double Q2 , 
			  const int NSIM ) 
{
  q2[ NSIM ] = sqrt( Q2 ) ; 
  t2[ NSIM ] = log( ( Q2 ) / ( mu * mu ) ) ;
 
  q2q3_set[ NSIM ] = true ;

  printf( "Q2[%d] SET !! %f \n" , NSIM , q2[NSIM] ) ; 

  return ;
}

// set the constant mu
void set_mu_D0_diff_v1( const double munew )
{
  mu = munew ;
  printf( "RENORMALISATION SCALE SET !!! %f \n" , mu ) ;
  return ;
}

// set the constant mu
void set_loops_D0_diff_v1( const int loopsnew )
{
  LOOP_ORDER = loopsnew ;
  printf( "NUMBER OF LOOPS SET !!! %d \n" , LOOP_ORDER ) ;
  return ;
}

// provide a description of the function
void
D0_diff_v1_description( const char *message ,
		const double *__restrict params ,
		const struct x_descriptor X ,
		const int NPARAMS )
{
  printf( "ALPHA_GUESS :: %e + ( a^2 q^2 ) %f  \n" , 
	  params[0] , params[1] ) ;
  printf( "Renormalisation scale :: %f \n" , mu ) ;
  return ;
}

// I provide a guess too
void
D0_diff_v1_guess( double *__restrict params ,
	       void *data ,
	       const int NPARAMS ) // nparams is already in data no?
{
  //const struct data* DATA = (const struct data*)data ;
  /*
  params[0] = 0.3 ;
  params[1] = 0.0 ;
  params[2] = 0.3 ;
  params[3] = 0.0 ;
  params[4] = 0.3 ;
  params[5] = 0.0 ;
  */
  params[0] = 0.3 ;
  params[1] = 0.0 ;
  params[2] = 0.0 ;

  /*
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
      p[ j ] = 0.003 ;
    }
    printf( "TEST2 :: %f %f \n" , DATA->quarks[0].ainverse , DATA->quarks[1].ainverse ) ;
    const struct mom_info dummy ;
    const struct x_descriptor XAXIS = { 0.0 , DATA->quarks[k] , dummy , 0 } ;
    #ifndef QUIET
    D0_diff_v1_description( "GUESS" , p , XAXIS , DATA-> NPARAMS ) ;
    #endif
    xposit += DATA->NDATA[k] ;
  }
  */
  return ;
}

static inline double
num( const double *fitparams , const double XX , const double XXref , const double a2 , const int der )
{
  return ( D0_OPE( log( ( XX ) / ( mu * mu ) ) , 
		   fitparams[0] / M_PI + fitparams[2]*a2 , 
		   der , LOOP_ORDER )
	   - D0_OPE( log( ( XXref ) / ( mu * mu ) ) , 
		     fitparams[0] / M_PI + fitparams[2]*a2 , 
		     der , LOOP_ORDER ) ) ;
}

// evaluate the function at the point "X" using the fitparams
inline double
D0_diff_v1_eval( const double *__restrict x ,
	      const struct x_descriptor X ,
	      const int NPARAMS )
{
  const double a2 = 1.0 / ( X.quark.ainverse * X.quark.ainverse ) ;
  const double mulfact = pisq4 / log( X.X * X.X / ( q2[X.LT] * q2[X.LT] ) ) ;
  return 
    mulfact * 
    ( num( x , X.X * X.X , q2[X.LT] * q2[X.LT] , a2 , 0 ) +
      x[1] * ( X.X*X.X - q2[X.LT]*q2[X.LT] )
      ) ; //- 1.0 ;
}

// evaluate the function and its difference from the data
int
D0_diff_v1_f( const gsl_vector *__restrict x , 
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

    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;

    size_t i  ;
    for( i = lopos ; i <= hipos ; i++ ) {

      const struct x_descriptor XAXIS = { DATA->X[i] , DATA->quarks[NSIM] , DATA->mom[i] , NSIM } ;

      // local parameters ...
      double p[ DATA->NPARAMS ] ;
 
      if( DATA -> sim_params[0] == true ) {
	//
	if( DATA -> sim_params[1] == true ) { 
	  p[0] = params[0] ;
	  p[1] = params[1] ;
	} else {
	  p[0] = params[0] ;
	  p[1] = params[loc_idx] ;
	}
	//
      } else if( DATA -> sim_params[1] == true ) { 
	if( NSIM == 0 ) {
	  p[0] = params[0] ;
	  p[1] = params[1] ;
	} else {
	  p[0] = params[loc_idx] ;
	  p[1] = params[1] ;
	}
      } else {
	p[0] = params[loc_idx] ;
	p[1] = params[loc_idx+1] ;
      }

      gsl_vector_set ( f , count , 
		       ( D0_diff_v1_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;

      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  free( (void*)params ) ;
  return GSL_SUCCESS;
}

// and the derivatives
int
D0_diff_v1_df( const gsl_vector *__restrict x, 
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

    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON)+DATA->NCOMMON ;
    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;

    // lattice spacings
    const double a = 1.0 / ( DATA->quarks[NSIM].ainverse ) ;
    const double a2 = a * a ;

    size_t i ; 
    for( i = lopos ; i <= hipos ; i++ ) {

      const double _s = 1.0 / DATA->sigma[i] ;
 
      const double p2 = DATA->X[i] * DATA->X[i] ;

      // initialise whole column to zero
      int k ;
      for( k = 0 ; k < DATA->LOGICAL_NPARS ; k++ ) {
	gsl_matrix_set( J , count , k ,  0.0 ) ;
      }

      // local parameters ...
      double p[ DATA->NPARAMS ] ;

      const double p2ref = q2[ NSIM ] * q2[ NSIM ] ;
      const double t = pisq4 * log( p2 / p2ref ) ;

      if( DATA -> sim_params[0] == true ) {
	p[0] = params[0] ;
	if( DATA -> sim_params[1] == true ) {
	  p[1] = params[1] ;
	} else {
	  p[1] = params[loc_idx] ;
	}
      } else if( DATA -> sim_params[1] == true ) {
	p[1] = params[1] ;
	if( NSIM == 0 ) { 
	  p[0] = params[0] ;
	} else {
	  p[0] = params[loc_idx] ;
	}
      } else {
	p[0] = params[loc_idx] ;
	p[1] = params[loc_idx+1] ;
      }


      const double expr1 = _s * num( p , p2 , p2ref , a2 , 1 ) * t ;
      const double expr2 = _s * ( p2 - p2ref ) * t ;

      if( DATA -> sim_params[0] == true && DATA -> sim_params[1] == true ) {
	gsl_matrix_set ( J , count , 0         , expr1 ) ;
	gsl_matrix_set ( J , count , 1         , expr2 ) ;
      } else if( DATA -> sim_params[0] == true ) {
	gsl_matrix_set ( J , count , 0         , expr1 ) ;
	gsl_matrix_set ( J , count , loc_idx   , expr2 ) ;
      } else if( DATA -> sim_params[1] == true ) {
	if( NSIM == 0 ) { 
	  gsl_matrix_set ( J , count , 0 , expr1 ) ;
	} else {
	  gsl_matrix_set ( J , count , loc_idx , expr1 ) ;
	}
	gsl_matrix_set ( J , count , 1         , expr2 ) ;
      } else {
	gsl_matrix_set ( J , count , loc_idx   , expr1 ) ;
	gsl_matrix_set ( J , count , loc_idx+1 , expr2 ) ;
      }

      count ++ ;
    }
    xposit += DATA->NDATA[NSIM] ; 
  }
  free( (void*)params ) ;

  return GSL_SUCCESS;
}

int
D0_diff_v1_fdf (const gsl_vector *__restrict x , 
	     void *data ,
	     gsl_vector *__restrict f , 
	     gsl_matrix *__restrict J )
{
  D0_diff_v1_f( x , data , f ) ;
  D0_diff_v1_df( x , data , J ) ;
  return GSL_SUCCESS ;
}
