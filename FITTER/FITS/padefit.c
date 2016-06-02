/**
   Pade coefficient fit
 */

#include "fitfunc.h"

#include "svd.h"
#include "Utils.h"

static int n = 1 , m = 1 ;

void 
pade_set_nm( const int new_n ,
	     const int new_m )
{
  if( new_n > 0 && new_m > 0 ) {
    n = new_n ;
    m = new_m ;
  }
}

void
pade_get_nm( int *new_n ,
	     int *new_m )
{
  *new_n = n ;
  *new_m = m ;
}

// provide a description of the function
void
pade_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS )
{
  // easiest one now
  printf( "%s :: \n" , message ) ;

  // numerator is first in our scheme
  printf( "                           ( " ) ;
  int i ;
  for( i = 0 ; i < n-1 ; i++ ) {
    printf( " (%1.2f) q^{%d} +" , params[ i + 1 ] , 2*(i+1) ) ;
  }
  printf( " (%1.2f) q^{%d} )\n" , params[ i + 1 ] , 2*(i+1) ) ;

  printf( "PADE( %d , %d ) :: %1.2f  +   " , n , m , params[ 0 ] ) ;
  for( i = 0 ; i < ( n>m?n:m ) ; i++ ) {
    printf( "-----------------" ) ;
  }
  printf( "---\n" ) ;

  // denominator is last
  printf( "                           ( 1 + " ) ;
  for( i = 0 ; i < m-1 ; i++ ) {
    printf( " (%1.2f) q^{%d} +" , params[ i + n + 1 ] , 2*(i+1) ) ;
  } 
  printf( " (%1.2f) q^{%d} ) \n" , params[ i + n + 1 ] , 2*(i+1) ) ;

  return ;
}

// I provide a guess too
void
pade_guess( double *__restrict params ,
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

    double ydata[ range ] , xs[ range ] , coeffs[ NPARAMS ] , sigma[ range ] , chisq ;

    int i , count = 0 ;
    for( i = lopos ; i <= hipos ; i++ ) {
      ydata[ count ] = DATA->y[ i ] ;
      sigma[ count ] = DATA->sigma[ i ] ;
      xs[ count ] = ( DATA->X[ i ] * DATA->X[ i ] ) ;
      count++ ;
    }
    xposit += ( DATA->NDATA[k] ) ;

    if( compute_coefficients( coeffs , &chisq , ydata , sigma , xs , 
			      range , NPARAMS ) )
      {
	// do something if it breaks ....
      }

    for( i = 0 ; i < NPARAMS ; i++ ) {
      printf( "[PADE] POLY-COEFFS %f \n" , coeffs[i] ) ;
    }

    // tells us about coefficients
    double p[ NPARAMS ] ;

    pades( p , coeffs , n , m ) ;

    struct x_descriptor X ;
    pade_description( "GUESS" , p , X , NPARAMS ) ;

    for( i = 0 ; i < NPARAMS ; i++ ) {
      params[ i ] = p[ i ] ;
    }
  }

  return ;
}

// evaluate the function at the point "X" using the fitparams
double
pade_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS )
{
  const double XX = X.X ; //* X.X ;
  register double numerator = XX * x[ n ] ;
  int i ;
  for( i = n-1 ; i > 0 ; i-- ) {
    numerator = XX * ( x[i] + numerator ) ;
  }
  register double denominator = XX * x[ n + m ] ;
  for( i = n + m - 1 ; i > n ; i-- ) {
    denominator = XX * ( x[i] + denominator ) ;
  }
  denominator += 1.0 ;
  return x[0] + ( numerator ) / ( denominator ) ;
}

// compute the difference of the function from the data
int
pade_f( const gsl_vector *__restrict x , 
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
		       ( pade_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;
  return GSL_SUCCESS;
}

int
pade_df( const gsl_vector *__restrict x, 
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
    for( i = lopos ; i <= hipos ; i++ ) {

      // initial constant term
      const double _s = 1.0 / DATA->sigma[i] ;

      // initialise whole column to zero
      int k ;
      for( k = 0 ; k < DATA->LOGICAL_NPARS ; k++ ) {
	gsl_matrix_set (J, count , k ,  0.0 ) ;
      }
      
      const double XX = DATA->X[i] ; //* DATA->X[i] ;

      // precompute numerator
      register double numerator = XX * params[ n ] ;
      for( k = n - 1 ; k > 0 ; k-- ) {
	numerator = XX * ( params[k] + numerator ) ;
      }

      // precompute denominator
      register double denominator = XX * params[ n + m ] ;
      for( k = n + m - 1 ; k > n ; k-- ) {
	denominator = XX * ( params[k] + denominator ) ;
      }
      denominator += 1.0 ;

      //printf( "Set?\n" ) ;
      const double factor = -( numerator ) / ( denominator * denominator ) ;

      gsl_matrix_set (J, count , loc_idx + 0 ,  1.0 * _s ) ;

      for( k = 1 ; k <= n ; k++ ) {
	gsl_matrix_set ( J , count , loc_idx + k , pow( XX , k ) / denominator  * _s ) ;
      }
      for( k = 0 ; k < m ; k++ ) {
	gsl_matrix_set ( J , count , loc_idx + n + k + 1 , pow( XX , k+1 ) * factor * _s ) ;
      }

      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;
  return GSL_SUCCESS ;
}

// and the function evaluations ...
int
pade_fdf (const gsl_vector *__restrict x , 
	  void *data ,
	  gsl_vector *__restrict f , 
	  gsl_matrix *__restrict J )
{
  pade_f( x , data , f ) ;
  pade_df( x , data , J ) ;
  return GSL_SUCCESS;
}


#if 0

      /*
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
      */

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
