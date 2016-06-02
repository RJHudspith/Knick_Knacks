/* expfit.c -- model functions for exponential + background */

#include "fitfunc.h"
#include "Utils.h"

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

// provide a description of the function
void
polystatalpha_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS )
{
  printf( "%s :: " , message ) ;
  printf( " ( a^2 %f ) + %f x + %f x^2 + %f x^3 \n" ,
	  params[0] , params[1] , params[2] , params[3] ) ;
  return ;
}

// I provide a guess too
void
polystatalpha_guess( double *__restrict params ,
	    void *data ,
	    const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;

  int k , xposit = 0 ;
  for( k = 0 ; k < DATA->SIMS ; k++ ) {

    const int lopos = find_idx( DATA->FIT_LO , DATA->X , 
				xposit + DATA->NDATA[k] , xposit ) ;

    const int NCOEFFS = 4 ;
    double coeffs[ NCOEFFS ] , xs[ NCOEFFS ] ;

    int i ;
    for( i = 0 ; i < NCOEFFS ; i++ ) {
      coeffs[i] = DATA->y[ xposit + lopos + i ] ;
      xs[i] = DATA->X[ xposit + lopos + i ] ;
    }

    solvandermonde( coeffs , xs , NCOEFFS ) ;
 
    double p[ NPARAMS ] ;

    // set these
    for( i = 0 ; i < NPARAMS ; i++ ) {
      if( i < NCOEFFS ) {
	p[i-1] = ( fabs( coeffs[i] ) > 100 ) ? 0.0 : coeffs[i] ;
      } else {
	p[i-1] = 0 ;
      }
      // set simultaneous parameters
      if( DATA -> sim_params[ i ] == true ) {
	params[ i ] += p[ i - 1 ] ;
      } else {
	params[ i + k * ( DATA->NPARAMS - DATA->NCOMMON ) ] = p[ i - 1 ] ;
      }
    }
  }

  // average shared params
  printf( "FIT coefficient guesses \n" ) ;
  for( k = 0 ; k < NPARAMS ; k++ ) {
    if( DATA -> sim_params[ k ] == true ) {
      params[ k ] /= DATA->SIMS ;
    }
    printf( "params %d :: %f \n" , k , params[k] ) ;
  }
  printf( "\n" ) ; 

  //exit(1) ;

  return ;
}

// evaluate the function at the point "X" using the fitparams
double
polystatalpha_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS )
{
  return ( x[0] / ( X.quark.ainverse * X.quark.ainverse ) ) + X.X * ( 
			    - ( -4.188790 * x[1] + x[1] * x[1] * ( -18.210642 + 9.424778 * log( 1.0 / ( 4. * X.X*X.X ) ) ) ) 
			    + X.X * ( x[2] + X.X * x[3] ) ) ;

}

// compute the difference of the function from the data
int
polystatalpha_f( const gsl_vector *__restrict x , 
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
    p[0] = params[0] ;
    p[1] = params[1] ;
    p[2] = params[loc_idx+1] ;
    p[3] = params[loc_idx+2] ;

    // loop each successive fit range of the flattened data
    size_t i ;
    for( i = lopos ; i < hipos ; i++ ) {

      // pack a local x-axis description
      const struct x_descriptor XAXIS = { DATA->X[i] , DATA->quarks[NSIM] } ;

      // evaluate fit function(s)
      gsl_vector_set ( f , count , 
		       ( polystatalpha_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;
  return GSL_SUCCESS;
}

int
polystatalpha_df ( const gsl_vector *__restrict x, 
	  void *data, 
          gsl_matrix *__restrict J )
{
  const struct data* DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

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

      gsl_matrix_set (J, count , 0 , _s / ( DATA->quarks[NSIM].ainverse * DATA->quarks[NSIM].ainverse ) ) ;
      gsl_matrix_set (J, count , 1 , -DATA->X[i] * ( -4.188790 + 2.0 * params[1] * ( -18.210642 + 
										     9.424778 * log( 1.0 / ( 4. * DATA->X[i]*DATA->X[i] ) ) ) ) * _s  ) ;
      gsl_matrix_set (J, count , loc_idx + 1 , pow( DATA->X[i] , 2 ) * _s ) ;
      gsl_matrix_set (J, count , loc_idx + 2 , pow( DATA->X[i] , 3 ) * _s ) ;

      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  return GSL_SUCCESS;
}

// and the function evaluations ...
int
polystatalpha_fdf (const gsl_vector *__restrict x , 
	   void *data ,
	   gsl_vector *__restrict f , 
	   gsl_matrix *__restrict J )
{
  polystatalpha_f( x , data , f ) ;
  polystatalpha_df( x , data , J ) ;
  return GSL_SUCCESS;
}
