/* expfit.c -- model functions for exponential + background */

#include "fitfunc.h"

#include "svd.h"
#include "Utils.h"

// provide a description of the function
void
poly_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS )
{
  int i ;
  printf( "%s :: " , message ) ;
  for( i = 0 ; i < NPARAMS-1 ; i++ ) {
    printf( "%f x^%d + " , params[i] , i ) ;
  }
  printf( "%f x^%d \n" , params[i] , i ) ;
  return ;
}

// I provide a guess too
void
poly_guess( double *__restrict params ,
	    void *data ,
	    const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;

  printf( "NPARAMS !!!! %d \n" , NPARAMS ) ;

  int k , xposit = 0 ;
  for( k = 0 ; k < DATA->SIMS ; k++ ) {

    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[k] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[k] , xposit ) ;
    const int range = hipos - lopos + 1 ;

    double ydata[ range ] , xs[ range ] , coeffs[ NPARAMS ] , chisq , sigma[ range ] ;

    int i , count = 0 ;
    for( i = lopos ; i <= hipos ; i++ ) {
      ydata[ count ] = DATA->y[ i ] ;
      sigma[ count ] = DATA->sigma[ i ] ;
      xs[ count ] = DATA->X[ i ] ;
      count++ ;
    }
    xposit += ( DATA->NDATA[k] ) ;

    if( compute_coefficients( coeffs , &chisq , ydata , sigma , xs , 
			      range , NPARAMS ) )
      {
	// do something if it breaks ....
      }

    // tells us about coefficients
    struct x_descriptor X ;
    poly_description( "GUESS" , coeffs , X , NPARAMS ) ;
    printf( "CHISQ guess :: %f \n" , chisq / NPARAMS ) ;

    double p[ NPARAMS ] ;

    // set these
    for( i = 0 ; i < NPARAMS ; i++ ) {
      if( i < NPARAMS ) {
	p[i] = ( fabs( coeffs[i] ) > 100 ) ? 0.0 : coeffs[i] ;
      } else {
	p[i] = 0 ;
      }
      // set simultaneous parameters
      if( DATA -> sim_params[ i ] == true ) {
	params[ i ] += p[ i ] ;
      } else {
	params[ i + k * ( DATA->NPARAMS - DATA->NCOMMON ) ] = p[ i ] ;
      }
    }
  }

  // average shared params
  printf( "FIT coefficient guesses \n" ) ;
  for( k = 0 ; k < DATA->SIMS * ( NPARAMS - DATA->NCOMMON ) + DATA->NCOMMON ; k++ ) {
    if( DATA -> sim_params[ k ] == true ) {
      params[ k ] /= DATA->SIMS ;
    }
    printf( "params %d :: %f \n" , k , params[k] ) ;
  }
  printf( "\n" ) ; 

  return ;
}

// evaluate the function at the point "X" using the fitparams
double
poly_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS )
{
  int i ;
  if( NPARAMS < 2 ) return x[0] ;
  register double poly = X.X * x[ NPARAMS-1 ] ;
  for( i = NPARAMS-2 ; i > 0 ; i-- ) {
    poly = X.X * ( x[i] + poly ) ;
  }
  return poly + x[0] ;
}

// compute the difference of the function from the data
int
poly_f( const gsl_vector *__restrict x , 
	void *data , 
	gsl_vector *__restrict f )
{
  const struct data* DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

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
    for( i = lopos ; i <= hipos ; i++ ) {

      // pack a local x-axis description
      const struct x_descriptor XAXIS = { DATA->X[i] } ;

      // evaluate fit function(s)
      gsl_vector_set ( f , count , 
		       ( poly_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;
  return GSL_SUCCESS;
}

int
poly_df ( const gsl_vector *__restrict x, 
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
    for( i = lopos ; i <= hipos ; i++ ) {

      // initial constant term
      const double _s = 1.0 / DATA->sigma[i] ;

      // initialise whole column to zero
      int k ;
      for( k = 0 ; k < DATA->LOGICAL_NPARS ; k++ ) {
	gsl_matrix_set (J, count , k ,  0.0 ) ;
      }

      double xloc = 1.0 ;
      if( NSIM > 0 ) {
	int check = DATA->NCOMMON ;
	for( k = 0 ; k < DATA->NPARAMS ; k++ ) {
	  if( DATA -> sim_params[k] == true ) {
	    gsl_matrix_set( J , count , k , xloc * _s ) ;
	  } else {
	    gsl_matrix_set( J , count , check+loc_idx , xloc * _s ) ;
	    check++ ;
	  }
	  xloc *= DATA->X[i] ;
	}
      } else {
	for( k = 0 ; k < DATA->NPARAMS ; k++ ) {
	  gsl_matrix_set( J , count , k , xloc * _s ) ;
	  xloc *= DATA->X[i] ;
	}
      }

      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  return GSL_SUCCESS;
}

// and the function evaluations ...
int
poly_fdf (const gsl_vector *__restrict x , 
	   void *data ,
	   gsl_vector *__restrict f , 
	   gsl_matrix *__restrict J )
{
  poly_f( x , data , f ) ;
  poly_df( x , data , J ) ;
  return GSL_SUCCESS;
}


#if 0
// use vandermonde to solve linear regression
static void 
solvandermonde( double *f ,
		const double *x ,
		const int N ) 
{
  /// compute the "newton representation" of interpolating polynomial ///
  int i , k ;
  for( k = 0 ; k < N-1 ; k++ ) { 
    for( i = N-1 ; i > k ; i-- ) {
      f[i] = ( f[i] - f[i-1] ) / ( x[i] - x[i-k-1] ) ;
    }
  }

  /// plug in values for the interpolator ///
  for( k = N-1 ; k >= 0 ; k-- ) {
    for(i = k ; i < N-1 ; i++ ) {
      f[i] = f[i] - ( f[ i + 1 ] * x[k] ) ;
    }
  }
  return;
}
#endif
