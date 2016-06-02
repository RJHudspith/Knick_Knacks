/**
   cosh fit

   y = params[1] * ( exp( -params[0] * t ) - exp( -params[0] * ( T - t ) ) ) 

   Simultaneously fit for params[0]

   Hmmm, need to provide "T" in data and X-descriptor, yuck
 */

#include "fitfunc.h"
#include "Utils.h"

// provide a description of the function
void
sinh_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS )
{
  printf( "%s :: y = %f * ( exp( -%f * x ) - exp( -%f * ( %d - x ) ) ) \n" ,
	  message , params[1] , params[0] , params[0] , X.LT ) ;
  return ;
}

// regression for the effective mass guesses
/*
  TODO :: run this backwards in the hope that we find the lightest meff first, subtract it
  and then the second and so on
 */
static bool
effmass_regress( double *mass , 
		 double *amp ,
		 const double *y , 
		 const double *X ,
		 const int NDATA ,
		 const double MASS_TOL )
{
  // effective mass solution
  int i = 0 , N = 0 ;
  double effmass_prev = -log( fabs( y[i+1] / y[i] ) ) / ( X[i+1] - X[i] ) ;
  double effmass_this = 0.0 ;

  int negs = 0 ;
  double sumx = 0.0 , sumy = 0.0 , sumxy = 0.0 , sumxx = 0.0 ;
  for( i = 1 ; i < NDATA/2 ; i++ ) {
    if( y[i] < 0.0 ) negs++ ;
    effmass_this = -log( fabs( y[i+1] / y[i] ) ) / ( X[i+1] - X[i] ) ;
    if( effmass_this < 0 && i > 2 ) break ;
    if( fabs( effmass_this - effmass_prev ) < MASS_TOL ) {
      // regression coefficients
      sumx += X[i] ;
      sumy += log( fabs( y[i] ) ) ;
      sumxx += X[i] * X[i] ;
      sumxy += X[i] * log( fabs( y[i] ) ) ;
      N++ ;
    }
    effmass_prev = effmass_this ;
  }

  // just in case
  if( N < 2 ) {
    printf( "None accepted ... Using previous parameters !! \n" ) ;
    return false ;
  }

  const double slope = ( sumx * sumy - N * sumxy ) / ( sumx * sumx - N * sumxx ) ;
  const double y_intercept = ( sumy - slope * sumx ) / ( N ) ;

  #if DEBUG
  printf( "sumx %f \n" , sumx ) ;
  printf( "sumy %f \n" , sumy ) ;
  printf( "sumxy %f \n" , sumxy ) ;
  printf( "sumxx %f \n" , sumxx ) ;
  printf( "%d \n" , N ) ;
  printf( "slope %f \n" , slope ) ;
  printf( "y-intercept %f \n" , y_intercept ) ;
  #endif

  *mass = -slope ;
  *amp = exp( y_intercept ) ;
  if( negs > NDATA/4 ) *amp *= -1 ;
  return true ;
}

// I provide a guess too
void
sinh_guess( double *__restrict params ,
	    void *data ,
	    const int NPARAMS )
{
  const struct data *DATA = (const struct data *)data ;

  //double summass = 0.0 ;
  params[0] = 1.0 ;
  int k , xposit = 0 ;
  for( k = 0 ; k < DATA->SIMS ; k++ ) {

    double ydata[ DATA -> NDATA[k] ] ;
    double xdata[ DATA -> NDATA[k] ] ;

    int i ;
    for( i = 0 ; i < DATA->NDATA[k] ; i++ ) {
      ydata[i] = DATA->y[ xposit + i ] ;
      xdata[i] = DATA->X[ xposit + i ] ;
    }
    xposit += DATA->NDATA[k] ;

    double p[ 2 ] = { 0 , 0 } ;
    if( effmass_regress( &p[0] , &p[1] ,
		     ydata , 
		     xdata ,
		     DATA->NDATA[k] ,
		     0.1 // this is a fudge parameter :: estimates how constant our meff is
			 ) == false ) {
      if( DATA->sim_params[0] == true ) {
	p[ 0 ] = params[0] ;
	p[ 1 ] = ( ydata[ 1 ] + ydata[2] ) / 
	  ( exp( -params[0] * xdata[1] ) + exp( -params[0] * xdata[2] )) ;
      }
    }
    /*
    summass += params[0] ;
    params[0] = summass / (double)( DATA->SIMS ) ;
    */
    if( DATA->sim_params[0] == true ) {
      params[0]   = p[0] ;
      params[k+1] = p[1] ;
    } else if( DATA->sim_params[1] == true )  {
      params[1] = p[1] ;
      params[k==0?0:k+1] = p[0] ;
    } else {
      params[2*k]   = p[0] ;
      params[2*k+1] = p[1] ;
    }

    struct mom_info momdummy ;
    struct chiral qdummy ;
    struct x_descriptor X = { 0 , qdummy , momdummy , DATA->LT } ;
    sinh_description( "Guess" , p , X , NPARAMS ) ;
  }
  return ;
}

// evaluate the function at the point "X" using the fitparams
inline double
sinh_eval( const double *__restrict x ,
	   const struct x_descriptor X ,
	   const int NPARAMS )
{
  return x[1] * ( exp( -x[0] * X.X ) - exp( -x[0] * ( X.LT - X.X ) ) ) ;
}

// compute the difference of the function from the data
int
sinh_f( const gsl_vector *__restrict x , 
	void *data , 
	gsl_vector *__restrict f )
{
  const struct data *DATA = (const struct data *)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  // dummy struct fillers
  struct mom_info dummy_mom ;
  struct chiral dummy_chiral ;

  int NSIM , count = 0 , xposit = 0 ;
  for( NSIM = 0 ; NSIM < DATA->SIMS ; NSIM++ ) {

    // local idx gives the unshared fit parameters
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON) ;

    // lo and hi for this x-range
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

    size_t i ;
    for( i = lopos ; i <= hipos ; i++ ) {
      const struct x_descriptor XAXIS = { DATA->X[i] , dummy_chiral , dummy_mom , DATA -> LT } ;
      gsl_vector_set ( f , count , 
		       ( sinh_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;

  return GSL_SUCCESS;
}

int
sinh_df( const gsl_vector *__restrict x, 
	 void *data, 
	 gsl_matrix *__restrict J )
{
  const struct data *DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  //printf( "%f \n" , LT_2 ) ;

  int count = 0 , NSIM , xposit = 0 ; 
  for( NSIM = 0 ; NSIM < DATA-> SIMS ; NSIM++ ) {
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON) ;
    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;

    size_t i ; 
    for( i = lopos ; i <= hipos ; i++ ) {
      const double _s = 1.0 / DATA->sigma[i] ;

      // initialise whole column to zero
      int k ;
      for( k = 0 ; k < DATA->LOGICAL_NPARS ; k++ ) {
	gsl_matrix_set (J, count , k ,  0.0 ) ;
      }

     // common factor cache
      const double forward = DATA->X[ i ] ;
      const double backward = (double)DATA->LT - DATA->X[ i ] - 1 ;

      // pretty hideous extension here
      double e1 , e2 ;
      if( DATA->sim_params[0] == true ) {
	e1 = exp( -params[0] * forward ) * _s ;
	e2 = exp( -params[0] * backward ) * _s ;
	gsl_matrix_set ( J , count , 0 ,  -params[loc_idx] * ( forward * e1 - backward * e2 ) ) ;
	gsl_matrix_set ( J , count , loc_idx+1 ,  e1 - e2  ) ;
      } else {
	if( NSIM == 0 ) {	
	  e1 = exp( -params[0] * forward ) * _s ;
	  e2 = exp( -params[0] * backward ) * _s ;
	  gsl_matrix_set ( J , count , 0 , -DATA->X[i] * params[1] * ( forward * e1 - backward * e2 ) ) ;
	  gsl_matrix_set ( J , count , 1 ,  e1 - e2  ) ;
	} else {
	  e1 = exp( -params[loc_idx+1] * forward ) * _s ;
	  e2 = exp( -params[loc_idx+1] * backward ) * _s ;
	  gsl_matrix_set ( J , count , loc_idx+1 , -DATA->X[i] * params[1] * ( forward * e1 - backward * e2 ) ) ;
	  gsl_matrix_set ( J , count , 1 ,  e1 - e2  ) ;
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
sinh_fdf(const gsl_vector *__restrict x , 
	 void *data ,
	 gsl_vector *__restrict f , 
	 gsl_matrix *__restrict J )
{
  sinh_f (x, data, f);
  sinh_df (x, data, J);

  return GSL_SUCCESS;
}
