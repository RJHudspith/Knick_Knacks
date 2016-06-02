/* expfit.c -- model functions for exponential + background */

#include "fitfunc.h"

#include "svd.h"
#include "Utils.h"

// provide a description of the function
void
chiral_cont_description( const char *message ,
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
chiral_cont_guess( double *__restrict params ,
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
    chiral_cont_description( "GUESS" , coeffs , X , NPARAMS ) ;

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
chiral_cont_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS )
{
  return x[0] + X.X * x[ 1 ] + x[ 2 ] / ( X.quark.ainverse * X.quark.ainverse ) ;
}

// compute the difference of the function from the data
int
chiral_cont_f( const gsl_vector *__restrict x , 
	void *data , 
	gsl_vector *__restrict f )
{
  const struct data* DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  int NSIM , count = 0 , xposit = 0 ;
  for( NSIM = 0 ; NSIM < DATA->SIMS ; NSIM++ ) { // loop simultaneous datasets

    // local idx gives the unshared fit parameters
    //const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON) ;
    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
 
    double p[ DATA->NPARAMS ] ;
    /*
    int j ;
    if( NSIM == 0 ) {
      for( j = 0 ; j < DATA -> NPARAMS-1 ; j++ ) {
	p[j] = params[j] ;
      }
    } else {
      int check = DATA->NCOMMON ;
      for( j = 0 ; j < DATA -> NPARAMS-1 ; j++ ) {
	if( DATA -> sim_params[j] == true ) {
	  p[j] = params[j] ;
	} else {
	  p[j] = params[ loc_idx + check ] ;
	  check++ ;
	}
      }
    }
    */
    p[ 0 ] = params[ 0 ] ;
    p[ 1 ] = params[ 1 ] ;
    p[ 2 ] = params[ 2 ] ;
    //p[ DATA -> NPARAMS-1 ] = params[ DATA -> NPARAMS-1 ] ;

    // loop each successive fit range of the flattened data
    size_t i ;
    for( i = lopos ; i <= hipos ; i++ ) {

      // pack a local x-axis description
      const struct x_descriptor XAXIS = { DATA->X[i] , DATA->quarks[ i ] } ;

      // evaluate fit function(s)
      gsl_vector_set ( f , count , 
		       ( chiral_cont_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;
  return GSL_SUCCESS;
}

int
chiral_cont_df ( const gsl_vector *__restrict x, 
	  void *data, 
          gsl_matrix *__restrict J )
{
  const struct data* DATA = (const struct data*)data ;

  int count = 0 , NSIM , xposit = 0 ; 
  for( NSIM = 0 ; NSIM < DATA-> SIMS ; NSIM++ ) {

    // the next local idx
    // const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON) ;

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

      /*
      double xloc = 1.0 ;
      if( NSIM > 0 ) {
	int check = DATA->NCOMMON ;
	for( k = 0 ; k < DATA->NPARAMS-1 ; k++ ) {
	  if( DATA -> sim_params[k] == true ) {
	    gsl_matrix_set( J , count , k , xloc * _s ) ;
	  } else {
	    gsl_matrix_set( J , count , check+loc_idx , xloc * _s ) ;
	    check++ ;
	  }
	  xloc *= DATA->X[i] ;
	}
      } else {
	for( k = 0 ; k < DATA->NPARAMS-1 ; k++ ) {
	  gsl_matrix_set( J , count , k , xloc * _s ) ;
	  xloc *= DATA->X[i] ;
	}
      }
     gsl_matrix_set( J , count , DATA->NPARAMS-1 , _s / ( DATA->quarks[i].ainverse * DATA->quarks[i].ainverse ) ) ;
      */

      gsl_matrix_set( J , count , 0 , _s  ) ;
      gsl_matrix_set( J , count , 1 , _s * DATA->X[i] ) ;
      gsl_matrix_set( J , count , 2 , _s / ( DATA->quarks[i].ainverse * DATA->quarks[i].ainverse )  ) ;

      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  return GSL_SUCCESS;
}

// and the function evaluations ...
int
chiral_cont_fdf (const gsl_vector *__restrict x , 
	   void *data ,
	   gsl_vector *__restrict f , 
	   gsl_matrix *__restrict J )
{
  chiral_cont_f( x , data , f ) ;
  chiral_cont_df( x , data , J ) ;
  return GSL_SUCCESS;
}



