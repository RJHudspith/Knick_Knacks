/**
   Gaussian distribution centered at x = 0

   y = params[1] * exp( -params[0] * x * x )

   Simultaneously fit for param[0] or param[1]
 */

#include "fitfunc.h"
#include "Utils.h"

// provide a description of the function
void
gaussian_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS )
{
  printf( "%s :: y = %f * ( exp( -%f * x^2 ) ) \n" ,
	  message , params[1] , params[0] ) ;
  return ;
}

// I provide a guess too
void
gaussian_guess( double *__restrict params ,
		void *data ,
		const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;
  int k = 0 , xposit = 0 ;
  double p[ 2 ] ;
  // simultaneous params are in here
  struct mom_info dummy ;
  const struct x_descriptor XX = { 0.0 , DATA->quarks[0] , dummy , DATA->LT } ;

  if( DATA -> sim_params[0] == true ) {
    // look for the best evaluation
    params[0] = 1.0 ; //set just in case the loop fails
    double BEST_EVAL = 1000. ;
    for( k = 1 ; k < DATA->LOGICAL_NPARS ; k++ ) {
      const int ZERO_IDX = find_idx( 0.0 , DATA->X , xposit + DATA->NDATA[k-1] , xposit ) ;
      const double diffy = ( DATA->y[ZERO_IDX] ) / ( DATA->y[ZERO_IDX + 1] ) ;
      const double diffx = ( DATA->X[ZERO_IDX] * DATA->X[ZERO_IDX] - DATA->X[ZERO_IDX + 1] * DATA->X[ZERO_IDX + 1] ) ;
      // search for the closest one to zero
      if( fabs( DATA-> X[ZERO_IDX] - 0.0 ) < BEST_EVAL 
	  && ( diffy ) > 0.0
	  && fabs( diffx ) > 0.0 ) 
	{
	  params[0] = -log( diffy ) / ( diffx ) ;
	}
      xposit += DATA->NDATA[k-1] ;
    }
    xposit = 0 ;
    for( k = 1 ; k < DATA->LOGICAL_NPARS ; k++ ) {
      const int ZERO_IDX = find_idx( 0.0 , DATA->X , xposit + DATA->NDATA[k-1] , xposit ) ;
      params[k] = DATA->y[ZERO_IDX] * exp( params[0] * DATA-> X[ZERO_IDX] * DATA-> X[ZERO_IDX] ) ;
      p[0] = params[0] ; p[1] = params[k] ;
      #ifndef QUIET
      gaussian_description( "GUESS" , p , XX , DATA->NPARAMS ) ;
      #endif
      xposit += DATA->NDATA[k-1] ;
    }
  } else {
    params[1] = 1000 ;
    {
      const int ZERO_IDX = find_idx( 0.0 , DATA->X , xposit + DATA->NDATA[k-1] , xposit ) ;
      const double diffy = ( DATA->y[ZERO_IDX] ) / ( DATA->y[ZERO_IDX + 1] ) ;
      const double diffx = ( DATA->X[ZERO_IDX] - DATA->X[ZERO_IDX + 1] ) ;
      params[0] = -log( diffy ) / ( diffx ) ;
      params[1] = DATA->y[xposit] * exp( params[0] * DATA-> X[xposit] * DATA-> X[xposit] ) ;
      p[0] = params[0] ; p[1] = params[1] ;
      #ifndef QUIET
      gaussian_description( "GUESS" , p , XX , DATA->NPARAMS ) ;
      #endif
      xposit += DATA->NDATA[0] ;
    }
    for( k = 2 ; k < DATA->LOGICAL_NPARS ; k++ ) {
      const int ZERO_IDX = find_idx( 0.0 , DATA->X , xposit + DATA->NDATA[k-1] , xposit ) ;
      const double diffy = ( DATA->y[ZERO_IDX] ) / ( DATA->y[ZERO_IDX + 1] ) ;
      const double diffx = ( DATA->X[ZERO_IDX] - DATA->X[ZERO_IDX + 1] ) ;
      params[k] = -log( diffy ) / ( diffx ) ;
      p[0] = params[k] ;
      p[1] = params[1] ;
      #ifndef QUIET
      gaussian_description( "GUESS" , p , XX , DATA->NPARAMS ) ;
      #endif
      xposit += DATA->NDATA[k-1] ;
    }
  }
  return ;
}

// evaluate the function at the point "X" using the fitparams
inline double
gaussian_eval( const double *__restrict x ,
	       const struct x_descriptor X ,
	       const int NPARAMS )
{
  return x[1] * exp( -x[0] * X.X * X.X ) ;
}

// compute the difference of the function from the data
int
gaussian_f( const gsl_vector *__restrict x , 
	    void *data , 
	    gsl_vector *__restrict f )
{
  const struct data *DATA = (const struct data *)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  int NSIM , count = 0 , xposit = 0 ;
  for( NSIM = 0 ; NSIM < DATA->SIMS ; NSIM++ ) {

    // local idx gives the unshared fit parameters
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON) ;

    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;

    // local parameters ...
    double p[ DATA->NPARAMS ] ;
    if( DATA->sim_params[0] == true ) {
      p[0] = params[0] ;
      p[1] = params[ loc_idx + 1 ] ;
    } else {
      p[ 1 ] = params[ 1 ] ;
      if( NSIM == 0 ) {
	p[ 0 ] = params[ 0 ] ;
      } else {
	p[ 0 ] = params[ loc_idx + 1 ] ;
      }
    }

    size_t i ;
    for( i = lopos ; i <= hipos ; i++ ) {
      const struct x_descriptor XAXIS = { DATA->X[i] } ;
      gsl_vector_set ( f , count , 
		       ( gaussian_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;

  return GSL_SUCCESS;
}

int
gaussian_df( const gsl_vector *__restrict x, 
	     void *data, 
	     gsl_matrix *__restrict J )
{
  const struct data *DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  int count = 0 , NSIM , xposit = 0 ; 
  for( NSIM = 0 ; NSIM < DATA-> SIMS ; NSIM++ ) {
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON) ;
    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    size_t i ; 
    for( i = lopos ; i <= hipos ; i++ ) {
      const double _s = 1.0 / DATA->sigma[i] ;
      //const double e = exp( -params[0] * DATA->X[i] * DATA->X[i] ) * _s ;

      // initialise whole column to zero
      int k ;
      for( k = 0 ; k < DATA->LOGICAL_NPARS ; k++ ) {
	gsl_matrix_set (J, count , k ,  0.0 ) ;
      }

      // pretty hideous extension here
      double e ;
      if( DATA->sim_params[0] == true ) {
	e = exp( -params[0] * DATA->X[i] * DATA->X[i] ) * _s ;
	gsl_matrix_set ( J , count , 0 ,  -DATA->X[i] * DATA->X[i] * params[loc_idx+1] * e ) ;
	gsl_matrix_set ( J , count , loc_idx+1 ,  e  ) ;
      } else {
	if( NSIM == 0 ) {
	  e = exp( -params[0] * DATA->X[i] * DATA->X[i] ) * _s ;
	  gsl_matrix_set ( J , count , 0 , -DATA->X[i] * DATA->X[i] * params[1] * e ) ;
	  gsl_matrix_set ( J , count , 1 , e ) ;
	} else {
	  e = exp( -params[loc_idx+1] * DATA->X[i] * DATA->X[i] ) * _s ;
	  gsl_matrix_set ( J , count , loc_idx+1 , -DATA->X[i] * DATA->X[i] * params[1] * e ) ;
	  gsl_matrix_set ( J , count , 1 , e ) ;
	}
      }
      count ++ ;
 
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  free( (void*)params ) ;

  return GSL_SUCCESS;
}

int
gaussian_fdf(const gsl_vector *__restrict x , 
	     void *data ,
	     gsl_vector *__restrict f , 
	     gsl_matrix *__restrict J )
{
  gaussian_f (x, data, f);
  gaussian_df (x, data, J);

  return GSL_SUCCESS;
}
