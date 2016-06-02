/* expfit.c -- model functions for exponential + background */

#include "fitfunc.h"
#include "Utils.h"

// provide a description of the function
void
staticpot_description( const char *message ,
		       const double *__restrict params ,
		       const struct x_descriptor X ,
		       const int NPARAMS )
{
  printf( "%s :: -%f/r + %f r + %f  \n" , message ,  
	  params[0] , params[1] , params[2] ) ;
  return ;
}

// I provide a guess too
void
staticpot_guess( double *__restrict params ,
		 void *data ,
		 const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;

  int xposit = 0 , NSIM = 0 ;
  for( NSIM = 0 ; NSIM < DATA->SIMS ; NSIM++ ) { // loop simultaneous datasets

    double loc_params[ 3 ] ;

    // get the low x evaluation first
    const double x = DATA->X[ xposit + 2 ] ;
    const double a = DATA->X[ xposit + 3 ] - x ;
    const double b = x - DATA->X[ xposit + 1 ] ;
    const double yxpa = DATA->y[ xposit + 3 ] ;
    const double yxmb = DATA->y[ xposit + 1 ] ;
    const double y = DATA->y[ xposit + 2 ] ;

    // approximate second derivative
    loc_params[ 0 ] = -x*x*x*( b * yxpa + a * yxmb - ( b + a ) * y ) / ( a*b*( a+b ) ) ;

    // compute the middle
    int mid , up ;
    mid = find_idx( 0.75 * ( DATA->X[ xposit + DATA->NDATA[NSIM] - 1 ] - DATA->X[ xposit ] ) , 
		    DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    up = find_idx( 0.95 * ( DATA->X[ xposit + DATA->NDATA[NSIM] - 1 ] - DATA->X[ xposit ] ) , 
		   DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;

    loc_params[ 1 ] = ( DATA->y[mid] - DATA->y[up] ) / ( DATA->X[mid] - DATA->X[up] ) ;
    loc_params[ 2 ] = ( DATA->X[up] * DATA->y[mid] - DATA->X[mid] * DATA->y[up] 
			+ loc_params[0] * ( DATA->X[up] / DATA->X[mid] - DATA->X[mid] / DATA->X[up]) ) \
      / ( DATA->X[up] - DATA->X[mid] ) ;

    params[ 0 + NSIM*3 ] = loc_params[0] ;
    params[ 1 + NSIM*3 ] = loc_params[1] ;
    params[ 2 + NSIM*3 ] = loc_params[2] ;

    const struct mom_info dummy ;
    struct x_descriptor XX = { 0.0 , DATA->quarks[NSIM] , dummy , DATA->LT } ;
    staticpot_description( "GUESS" , loc_params , XX , DATA->NPARAMS ) ;

    xposit += DATA->NDATA[NSIM] ;
  }

  return ;
}

// evaluate the function at the point "X" using the fitparams
double
staticpot_eval( const double *__restrict x ,
		const struct x_descriptor X ,
		const int NPARAMS )
{
  return -x[0]/X.X + x[1] * X.X + x[2] ;
}

// compute the difference of the function from the data
int
staticpot_f( const gsl_vector *__restrict x , 
	     void *data , 
	     gsl_vector *__restrict f )
{
  const struct data* DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  const struct mom_info dummy_mom ;

  int NSIM , count = 0 , xposit = 0 ;
  for( NSIM = 0 ; NSIM < DATA->SIMS ; NSIM++ ) { // loop simultaneous datasets

    // local idx gives the unshared fit parameters
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON) ;
    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;


    double p[ DATA->NPARAMS ] ;
    int j ;
    if( NSIM == 0 ) {
      for( j = 0 ; j < DATA -> NPARAMS ; j++ ) {
	p[j] = params[j] ;
      }
    } else {
      int check = DATA->NCOMMON ;
      for( j = 0 ; j < DATA -> NPARAMS ; j++ ) {
	if( DATA -> sim_params[j] == true ) {
	  p[j] = params[j] ;
	} else {
	  p[j] = params[ loc_idx + check ] ;
	  check++ ;
	}
      }
    }

    // loop each successive fit range of the flattened data
    size_t i ;
    for( i = lopos ; i < hipos ; i++ ) {

      // pack a local x-axis description
      const struct x_descriptor XAXIS = { DATA->X[i] , DATA->quarks[NSIM] , dummy_mom , 0 } ;

      // evaluate fit function(s)
      gsl_vector_set ( f , count , 
		       ( staticpot_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;
  return GSL_SUCCESS;
}

int
staticpot_df( const gsl_vector *__restrict x, 
	      void *data, 
	      gsl_matrix *__restrict J )
{
  const struct data* DATA = (const struct data*)data ;

  int count = 0 , NSIM , xposit = 0 ; 
  for( NSIM = 0 ; NSIM < DATA-> SIMS ; NSIM++ ) {

    // the next local idx
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON) ;

    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;

    size_t i ;
    for( i = lopos ; i < hipos ; i++ ) {

      // initial constant term
      const double _s = 1.0 / DATA->sigma[i] ;

      // initialise whole column to zero
      int k ;
      for( k = 0 ; k < DATA->LOGICAL_NPARS ; k++ ) {
	gsl_matrix_set (J, count , k ,  0.0 ) ;
      }

      // 
      gsl_matrix_set (J, count , loc_idx , -_s / DATA->X[i] ) ;
      gsl_matrix_set (J, count , loc_idx + 1 , _s * DATA->X[i] ) ;
      gsl_matrix_set (J, count , loc_idx + 2 , _s ) ;

      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  return GSL_SUCCESS;
}

// and the function evaluations ...
int
staticpot_fdf( const gsl_vector *__restrict x , 
	       void *data ,
	       gsl_vector *__restrict f , 
	       gsl_matrix *__restrict J )
{
  staticpot_f( x , data , f ) ;
  staticpot_df( x , data , J ) ;
  return GSL_SUCCESS;
}
