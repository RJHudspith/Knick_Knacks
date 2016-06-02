/**
   @brief chiral fit forms for bags and ratios
 */

#include "fitfunc.h"

#include "svd.h"
#include "Utils.h"

// chiral scale in GeV
const static double chiral = 0.77 ;

static double fac = 3. / 2. ;

//const static double b02_f2 = 2 * 4.236 / ( 0.12229 * 0.12229 ) ;

//const static double _f2 = 1.0 / ( 0.12229 * 0.12229 ) ;

void
chiral_cont_R23_setfac( const double factor )
{
  fac = factor ;
  return ;
}

// provide a description of the function
void
chiral_cont_R23_description( const char *message ,
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
chiral_cont_R23_guess( double *__restrict params ,
		       void *data ,
		       const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;

  int k , xposit = 0 ;
  for( k = 0 ; k < DATA->SIMS ; k++ ) {

    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[k] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[k] , xposit ) ;
    const int range = hipos - lopos + 1 ;

    double ydata[ range ] , xs[ range ] , coeffs[ 2 ] , sigma[ range ] , chisq ;

    int i , count = 0 ;
    for( i = lopos ; i <= hipos ; i++ ) {
      ydata[ count ] = DATA->y[ i ] ;
      sigma[ count ] = DATA->sigma[ i ] ;
      xs[ count ] = DATA->X[ i ] ;
      count++ ;
    }
    xposit += ( DATA->NDATA[k] ) ;

    if( compute_coefficients( coeffs , &chisq , ydata , sigma , xs , 
			      range , 2 ) )
      {
	// do something if it breaks ....
      }

    // tells us about coefficients
    struct x_descriptor X ;
    chiral_cont_R23_description( "GUESS" , coeffs , X , 2 ) ;

    double p[ NPARAMS ] ;

    // set these
    for( i = 0 ; i < NPARAMS ; i++ ) {
      if( i < 2 ) {
	p[0] = ( fabs( coeffs[i] ) > 100 ) ? 0.0 : coeffs[0] ;
	p[1] = ( fabs( coeffs[i] ) > 100 ) ? 0.0 : coeffs[1] ;
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
chiral_cont_R23_eval( const double *__restrict x ,
		      const struct x_descriptor X ,
		      const int NPARAMS )
{
  return x[0] * ( 1 + 
		  X.X * ( x[ 1 ] + fac * ( log( X.X / ( chiral * chiral ) ) / ( 16 * M_PI * M_PI ) ) ) +  
		  x[ 2 ] / ( X.quark.ainverse * X.quark.ainverse ) ) ;
}

// compute the difference of the function from the data
int
chiral_cont_R23_f( const gsl_vector *__restrict x , 
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
    p[ 0 ] = params[ 0 ] ;
    p[ 1 ] = params[ 1 ] ;
    p[ 2 ] = params[ 2 ] ;

    // loop each successive fit range of the flattened data
    size_t i ;
    for( i = lopos ; i <= hipos ; i++ ) {

      // pack a local x-axis description
      const struct x_descriptor XAXIS = { DATA->X[i] , DATA->quarks[ i ] } ;

      // evaluate fit function(s)
      gsl_vector_set ( f , count , 
		       ( chiral_cont_R23_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;
  return GSL_SUCCESS;
}

int
chiral_cont_R23_df ( const gsl_vector *__restrict x, 
		     void *data, 
		     gsl_matrix *__restrict J )
{
  const struct data* DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  int count = 0 , NSIM , xposit = 0 ; 
  for( NSIM = 0 ; NSIM < DATA-> SIMS ; NSIM++ ) {

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

      gsl_matrix_set( J , count , 0 , _s * ( 1 + DATA->X[i] * ( params[ 1 ] + fac * ( log( DATA->X[i] / ( chiral * chiral ) ) / ( 16 * M_PI * M_PI ) ) ) + params[2] / ( DATA->quarks[i].ainverse * DATA->quarks[i].ainverse ) ) ) ;
      gsl_matrix_set( J , count , 1 , _s * DATA->X[i] * params[0] ) ;
      gsl_matrix_set( J , count , 2 , _s * params[0] / ( DATA->quarks[i].ainverse * DATA->quarks[i].ainverse )  ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  return GSL_SUCCESS;
}

// and the function evaluations ...
int
chiral_cont_R23_fdf (const gsl_vector *__restrict x , 
		     void *data ,
		     gsl_vector *__restrict f , 
		     gsl_matrix *__restrict J )
{
  chiral_cont_R23_f( x , data , f ) ;
  chiral_cont_R23_df( x , data , J ) ;
  return GSL_SUCCESS;
}



