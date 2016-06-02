/**
   @brief chiral fit forms for bags and ratios
 */

#include "fitfunc.h"

#include "svd.h"
#include "Utils.h"

// chiral scale in GeV
const static double chiral = 0.77 ;

static double fac = 3. / 2. ;
//const static double _f2 = 1.0 / ( 0.12229 * 0.12229 ) ;

void
chiral_R23_setfac( const double factor )
{
  fac = factor ;
  return ;
}

// provide a description of the function
void
chiral_R23_description( const char *message ,
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

static double *
compute_guess( const struct data* DATA ,
	       const int xposit ,
	       const int k ) 
{
  double *coeffs = malloc( 2 * sizeof( double ) ) ;
  const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[k] , xposit ) ;
  const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[k] , xposit ) ;
  const int range = hipos - lopos + 1 ;
  double ydata[ range ] , xs[ range ] , sigma[ range ] , chisq ;
  int i , count = 0 ;
  for( i = lopos ; i <= hipos ; i++ ) {
    ydata[ count ] = DATA->y[ i ] ;
    sigma[ count ] = DATA->sigma[ i ] ;
    xs[ count ] = DATA->X[ i ] ;
    count++ ;
  }
  compute_coefficients( coeffs , &chisq , ydata , sigma , xs , range , 2 ) ;
  printf( "%f %f \n" , coeffs[0] , coeffs[1] ) ;
  return coeffs ;
}

// I provide a guess too
void
chiral_R23_guess( double *__restrict params ,
		  void *data ,
		  const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;
  int xposit = 0 ;
  {
    double *coeffs1 = compute_guess( DATA , xposit , 0 ) ;
    xposit += ( DATA->NDATA[0] ) ;
    double *coeffs2 = compute_guess( DATA , xposit , 1 ) ;
    xposit += ( DATA->NDATA[1] ) ;
    params[0] = ( coeffs1[0] + coeffs2[0] ) / 2.0 ;
    params[1] = ( coeffs1[1] + coeffs2[ 1 ] ) / ( 2.0 * params[ 0 ] )  ;
    params[2] = -1.0 ;
    free( coeffs1 ) ; free( coeffs2 ) ;
  }
  {
    double *coeffs1 = compute_guess( DATA , xposit , 0 ) ;
    xposit += ( DATA->NDATA[0] ) ;
    double *coeffs2 = compute_guess( DATA , xposit , 1 ) ;
    xposit += ( DATA->NDATA[1] ) ;
    params[3] = ( coeffs1[0] + coeffs2[0] ) / 2.0 ;
    params[1] = ( params[1] + ( coeffs1[1] + coeffs2[ 1 ] ) / ( 2.0 * params[ 3 ] ) ) / 2.0 ;
    params[4] = -1.0 ;
    free( coeffs1 ) ; free( coeffs2 ) ;
  }
  printf( "FIT coefficient guesses \n" ) ;
  size_t k ;
  for( k = 0 ; k < DATA->SIMS * ( NPARAMS - DATA->NCOMMON ) + DATA->NCOMMON ; k++ ) {
    printf( "param %zu :: %f \n" , k , params[k] ) ;
  }

  return ;
}

// evaluate the function at the point "X" using the fitparams
double
chiral_R23_eval( const double *__restrict x ,
		 const struct x_descriptor X ,
		 const int NPARAMS )
{
  return x[0] * ( 1.0 + 
		  X.X * ( x[ 1 ] + fac * ( log( X.X / ( chiral * chiral ) ) / ( 16 * M_PI * M_PI ) ) ) + 
		  x[ 2 ] / ( X.quark.ainverse * X.quark.ainverse ) ) ;
}

// compute the difference of the function from the data
int
chiral_R23_f( const gsl_vector *__restrict x , 
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

    if( NSIM < DATA->SIMS/2 ) {
      p[ 0 ] = params[ 0 ] ;
      p[ 1 ] = params[ 1 ] ;
      p[ 2 ] = params[ 2 ] ;
    } else {
      p[ 0 ] = params[ 3 ] ;
      p[ 1 ] = params[ 1 ] ;
      p[ 2 ] = params[ 4 ] ;
    }

    // loop each successive fit range of the flattened data
    size_t i ;
    for( i = lopos ; i <= hipos ; i++ ) {

      // pack a local x-axis description
      const struct x_descriptor XAXIS = { DATA->X[i] , DATA->quarks[ i ] } ;

      // evaluate fit function(s)
      gsl_vector_set ( f , count , 
		       ( chiral_R23_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;
  return GSL_SUCCESS;
}

int
chiral_R23_df( const gsl_vector *__restrict x, 
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

      // fit R1
      if( NSIM < DATA-> SIMS/2 ) {
	gsl_matrix_set( J , count , 0 , _s * ( 1 
					       + DATA->X[i] * ( params[ 1 ] + fac * ( log( DATA->X[i] / ( chiral * chiral ) ) / ( 16 * M_PI * M_PI ) ) ) 
					       + params[2] / ( DATA->quarks[i].ainverse * DATA->quarks[i].ainverse )
					       ) ) ;
	gsl_matrix_set( J , count , 1 , _s * DATA->X[i] * params[0] ) ;
	gsl_matrix_set( J , count , 2 , _s * params[0] / ( DATA->quarks[i].ainverse * DATA->quarks[i].ainverse )  ) ;
	// fit R2
      } else if( NSIM/2 == 1 ) {
	gsl_matrix_set( J , count , 3 , _s * ( 1 
					       + DATA->X[i] * ( params[ 1 ] + fac * ( log( DATA->X[i] / ( chiral * chiral ) ) / ( 16 * M_PI * M_PI ) ) ) 
					       + params[4] / ( DATA->quarks[i].ainverse * DATA->quarks[i].ainverse )
					       ) ) ;
	gsl_matrix_set( J , count , 1 , _s * DATA->X[i] * params[3] ) ;
	gsl_matrix_set( J , count , 4 , _s * params[3] / ( DATA->quarks[i].ainverse * DATA->quarks[i].ainverse )  ) ;
	// fit R3
      } 
      
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  return GSL_SUCCESS;
}

// and the function evaluations ...
int
chiral_R23_fdf( const gsl_vector *__restrict x , 
		void *data ,
		gsl_vector *__restrict f , 
		gsl_matrix *__restrict J )
{
  chiral_R23_f( x , data , f ) ;
  chiral_R23_df( x , data , J ) ;
  return GSL_SUCCESS;
}
