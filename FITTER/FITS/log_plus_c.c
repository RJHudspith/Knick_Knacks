/*
  Log plus constant

  y( x ) = p[0] + p[1] log( x )
 */

#include "fitfunc.h"
#include "Utils.h"

// provide a description of the function
void
log_plus_c_description( const char *message ,
			const double *__restrict params ,
			const struct x_descriptor X ,
			const int NPARAMS )
{
  printf( "%s :: y = %f + %f Log( x ) \n" ,
	  message , params[0] , params[1] ) ;
  return ;
}

// I provide a guess too
void
log_plus_c_guess( double *__restrict params ,
		  void *data ,
		  const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;
  double p[2] ;
  struct mom_info dummy ;
  struct x_descriptor XAXIS = { 0.0 , DATA->quarks[0] , dummy , 0 } ;
  int k , xposit = 0 ;

  if( DATA -> sim_params[0] == true ) {
    xposit = 0 ;
    for( k = 1 ; k < DATA->LOGICAL_NPARS ; k++ ) {
      const int loc_end = xposit + DATA->NDATA[k-1] - 2 ;
      params[k] = ( DATA->y[loc_end+1] - DATA->y[loc_end] ) / ( log( DATA->X[loc_end+1] ) - log( DATA->X[loc_end] ) ) ;
      p[0] = DATA->y[loc_end+1] - params[k] * log( DATA->X[loc_end+1] ) ; 
      p[1] = params[k] ;
      #ifndef QUIET
      log_plus_c_description( "GUESS" , p , XAXIS , DATA->NPARAMS ) ;
      #endif
      xposit += DATA->NDATA[k-1] ;
    }
  } else {
    // new version
    // (y2-y1)/(log(x2)-log(x1)) = p[1]
    xposit = 0 ;
    params[1] = ( DATA->y[2] - DATA->y[1] ) / ( log( DATA->X[2] ) - log( DATA->X[1] ) ) ;

    {
      p[0] = params[0] = DATA->y[2] - params[1] * log( DATA->X[2] ) ; 
      p[1] = params[1] ;
      #ifndef QUIET
      log_plus_c_description( "GUESS" , p , XAXIS , DATA->NPARAMS ) ;
      #endif
      xposit += DATA->NDATA[0] ;
    }

    for( k = 2 ; k < DATA->LOGICAL_NPARS ; k++ ) {
      const int loc_end = xposit + DATA->NDATA[k-1] - 2 ;
      p[0] = params[k] = DATA->y[loc_end+1] - params[1] * log( DATA->X[loc_end+1] ) ; 
      p[1] = params[1] ;
      #ifndef QUIET
      log_plus_c_description( "GUESS" , p , XAXIS , DATA->NPARAMS ) ;
      #endif
      xposit += DATA->NDATA[k-1] ;
    }
  }
  return ;
}

// evaluate the function at the point "X" using the fitparams
inline double
log_plus_c_eval( const double *__restrict x ,
		 const struct x_descriptor X ,
		 const int NPARAMS )
{
  return x[0] + x[1] * log( X.X )  ;
}

// compute the difference of the function from the data
int
log_plus_c_f( const gsl_vector *__restrict x , 
	      void *data , 
	      gsl_vector *__restrict f )
{
  const struct data* DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  int NSIM , count = 0 , xposit = 0 ;
  for( NSIM = 0 ; NSIM < DATA->SIMS ; NSIM++ ) {
    size_t i  ;
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON) ;

    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;

    // local parameters ...
    double p[ DATA->NPARAMS ] ;
    if( DATA->sim_params[0] == true ) {
      p[0] = params[0] ;
      p[1] = params[ loc_idx+1 ] ;
    } else {
      p[ 1 ] = params[ 1 ] ;
      if( NSIM == 0 ) {
	p[ 0 ] = params[ 0 ] ;
      } else {
	p[ 0 ] = params[ loc_idx + 1 ] ;
      }
    }
 
    for( i = lopos ; i <= hipos ; i++ ) {
      const struct x_descriptor XAXIS = { DATA->X[i] } ;
      gsl_vector_set ( f , count , 
		       ( log_plus_c_eval( p , XAXIS , 2 ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  free( (void*)params ) ;

  return GSL_SUCCESS;
}

// derivatives
int
log_plus_c_df( const gsl_vector *__restrict x, 
	       void *data, 
	       gsl_matrix *__restrict J )
{
  const struct data* DATA = (const struct data*)data ;

  int count = 0 , NSIM , xposit = 0 ; 
  for( NSIM = 0 ; NSIM < DATA-> SIMS ; NSIM++ ) {
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON) ;
    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    size_t i ;
    for( i = lopos ; i <= hipos ; i++ ) {
      const double _s = 1.0 / DATA->sigma[i] ;

      int k ;
      for( k = 0 ; k < DATA->LOGICAL_NPARS ; k++ ) {
	gsl_matrix_set (J, count , k ,  0.0 ) ;
      }

      if( DATA->sim_params[0] == true ) {
	gsl_matrix_set ( J , count , 0 ,  _s ) ;
	gsl_matrix_set ( J , count , loc_idx+1 , log( DATA->X[i] ) * _s ) ;
      } else {
	if( NSIM == 0 ) {
	  gsl_matrix_set ( J , count , 0 , _s ) ;
	  gsl_matrix_set ( J , count , 1 , log( DATA->X[i] ) * _s ) ;
	} else {
	  gsl_matrix_set ( J , count , loc_idx+1 , _s ) ;
	  gsl_matrix_set ( J , count , 1 , log( DATA->X[i] ) * _s ) ;
	}
      }
      count ++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  return GSL_SUCCESS;
}

int
log_plus_c_fdf (const gsl_vector *__restrict x , 
		void *data ,
		gsl_vector *__restrict f , 
		gsl_matrix *__restrict J )
{
  log_plus_c_f( x , data , f ) ;
  log_plus_c_df( x , data , J ) ;
  return GSL_SUCCESS;
}
