/**
   Two parameter fit to

   y = b / ( x + a )

   ders
   dy/da = -b / ( x + a )^2
   dy/db = 1 / ( x + a )
 */
#include "fitfunc.h"
#include "Utils.h"

// provide a description of the function
void
prop_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS )
{
  printf( "%s :: y = %f / ( x + %f ) \n" ,
	  message , params[1] , params[0] ) ;
  return ;
}

// I provide a guess too
void
prop_guess( double *__restrict params ,
	    void *data ,
	    const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;
  int k , xposit = 0 ;

  // params[0] is the mass
  double BEST = 1000. ;
  int BEST_IDX = 0 , end = DATA->n-1 ;
  for( k = 0 ; k < DATA->SIMS ; k++ ) {
    const int ZERO_IDX = find_idx( 0.0 , DATA->X , xposit + DATA->NDATA[k] , xposit ) ; // find the value closest to zero out of all the data
    if( fabs( DATA->X[ZERO_IDX] ) < BEST ) {
      BEST_IDX = ZERO_IDX ;
      BEST = fabs( DATA->X[ZERO_IDX] ) ;
      end = xposit + DATA->NDATA[k] - 1 ;
    }
    xposit += DATA->NDATA[k] ;
  }
  printf( "XBEST :: %f %f \n" , DATA->X[ BEST_IDX ] , DATA->X[ BEST_IDX + 1 ]  ) ;
  printf( "BEST_IDX %d END %d %f \n" , BEST_IDX , end , DATA->X[ end ] ) ;

  // test out new equation
  params[0] = ( DATA->y[ BEST_IDX ] * DATA->X[ BEST_IDX ] - DATA->y[ end ] * DATA->X[ end ] ) 
    / ( DATA->y[ end ] - DATA->y[ BEST_IDX ] ) ;

  // stuff for the print function
  double p[ 2 ] ;
  struct mom_info dummy ;
  struct x_descriptor XX = { 0.0 , DATA->quarks[0] , dummy , DATA->LT } ;

  xposit = 0 ;
  for( k = 1 ; k < DATA->LOGICAL_NPARS ; k++ ) {
    params[k] = ( params[0] + DATA -> X[ xposit ] ) * ( DATA -> y[ xposit ] ) ;
    p[0] = params[0] ; p[1] = params[k] ;
#ifndef QUIET
    prop_description( "GUESS" , p , XX , DATA->NPARAMS ) ;
#endif
    xposit += DATA -> NDATA[ k-1 ] ;
  }

  return ;
}

// evaluate the function at the point "X" using the fitparams
inline double
prop_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS )
{
  return x[1] / ( X.X + x[0] ) ;
}

// compute the difference of the function from the data
int
prop_f( const gsl_vector *__restrict x , 
	void *data , 
	gsl_vector *__restrict f )
{
  const struct data* DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  int NSIM , count = 0 , xposit = 0 ;
  for( NSIM = 0 ; NSIM < DATA->SIMS ; NSIM++ ) {

    // local idx gives the unshared fit parameters
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON)+DATA->NCOMMON ;

    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;

    // local parameters ...
    double p[ DATA->NPARAMS ] ;
    p[0] = params[0] ;
    p[1] = params[loc_idx] ;

    size_t i ;
    for( i = lopos ; i <= hipos ; i++ ) {
      const struct x_descriptor XAXIS = { DATA->X[i] } ;
      gsl_vector_set ( f , count , 
		       ( prop_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA -> NDATA[ NSIM ] ;
  }
  free( (void*)params ) ;

  return GSL_SUCCESS;
}

// the derivatives
int
prop_df( const gsl_vector *__restrict x, 
	 void *data, 
	 gsl_matrix *__restrict J )
{
  const struct data* DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  int count = 0 , NSIM , xposit = 0 ; 
  for( NSIM = 0 ; NSIM < DATA-> SIMS ; NSIM++ ) {
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON)+DATA->NCOMMON ;
    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    size_t i ; 
    for( i = lopos ; i <= hipos ; i++ ) {
      const double _s = 1.0 / DATA->sigma[i] ;
      const double denom = 1.0 / ( DATA->X[i] + params[0] ) ;

      // initialise whole column to zero
      int k ;
      for( k = 0 ; k < DATA->LOGICAL_NPARS ; k++ ) {
	gsl_matrix_set (J, count , k ,  0.0 ) ;
      }

      gsl_matrix_set ( J , count , 0 , -params[loc_idx] * denom * denom * _s ) ;
      gsl_matrix_set ( J , count , loc_idx , denom * _s ) ;

      count ++ ;
    }
    xposit += DATA -> NDATA[NSIM] ;
  }
  free( (void*)params ) ;
  return GSL_SUCCESS;
}

// wrapper for the function and the derivatives
int
prop_fdf( const gsl_vector *__restrict x , 
	  void *data ,
          gsl_vector *__restrict f , 
	  gsl_matrix *__restrict J )
{
  prop_f(x, data, f);
  prop_df(x, data, J);

  return GSL_SUCCESS;
}
