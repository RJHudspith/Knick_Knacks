/**
   Single exponential fit

   y = params[1] * exp( -params[0] * x )

   Simultaneously fit for params[0] or params[1], if that
   is what you're into

   And compare to UKhadron ... Eventually
 */

#include "fitfunc.h"
#include "Utils.h"

// provide a description of the function
void
exp_description( const char *message ,
		 const double *__restrict params ,
		 const struct x_descriptor X ,
		 const int NPARAMS )
{
  printf( "%s :: y = %f * ( exp( -%f * x ) ) \n" ,
	  message , params[1] , params[0] ) ;
  return ;
}

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
  double effmass_prev = -log( y[i+1] / y[i] ) / ( X[i+1] - X[i] ) ;
  double effmass_this = 0.0 ;

  double sumx = 0.0 , sumy = 0.0 , sumxy = 0.0 , sumxx = 0.0 ;
  for( i = 1 ; i < NDATA ; i++ ) {
    if( y[i] < 0.0 ) break ;
    effmass_this = -log( fabs( y[i+1] / y[i] ) ) / ( X[i+1] - X[i] ) ;
    if( effmass_this < 0 && i > 2 ) break ;
    if( fabs( effmass_this - effmass_prev ) < MASS_TOL ) {
      printf( "ACCEPTED :: %f \n" , effmass_this ) ;
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
  const double y_intercept = ( sumy - slope * sumx) / ( N ) ;

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
  return true ;
}

// I provide a guess too
void
exp_guess( double *__restrict params ,
	   void *data ,
	   const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;
  double tempdata[ DATA -> NDATA[0] ] ;

  int i ;
  for( i = 0 ; i < DATA -> NDATA[0] ; i++ ) {
    tempdata[ i ] = DATA -> y[ i ] ;
  }

  // todo loop k
  int k ;
  for( k = 0 ; k < DATA->SIMS ; k++ ) {
    double p[ 2 ] = { 0 , 0 } ;
    effmass_regress( &p[0] , &p[1] ,
		     tempdata , 
		     DATA->X ,
		     DATA->NDATA[k] ,
		     0.075 + 0.66 * 1 ) ;
    if( DATA->sim_params[0] == true ) {
      params[0] += p[0] ;
      params[1+k] = p[1] ;
    } else if( DATA->sim_params[1] == true ) {
      params[0] += p[1] ;
      params[k+1] = p[0] ;
    } else {
      params[k+0] = p[0] ;
      params[k+1] = p[1] ;
    }

    struct x_descriptor X ;
    exp_description( "Guess" , p , X , NPARAMS ) ;
  }
  if( DATA->sim_params[0] == true || DATA->sim_params[1] == true ) {
    params[0] /= ( DATA -> SIMS ) ;
  }

  return ;
}

// evaluate the function at the point "X" using the fitparams
inline double
exp_eval( const double *__restrict x ,
	  const struct x_descriptor X ,
	  const int NPARAMS )
{
  return x[1] * exp( -x[0] * X.X ) ;
}

// compute the difference of the function from the data
int
exp_f( const gsl_vector *__restrict x , 
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
      p[0] = params[ 0 ] ;
      p[1] = params[ loc_idx+1 ] ;
    } else if ( DATA->sim_params[1] == true ) {
      p[0] = params[ loc_idx+1 ] ;
      p[1] = params[ 0 ] ;
    } else {
      p[0] = params[ loc_idx+0 ] ;
      p[1] = params[ loc_idx+1 ] ;
    }

    size_t i ;
    for( i = lopos ; i <= hipos ; i++ ) {
      const struct x_descriptor XAXIS = { DATA->X[i] } ;
      gsl_vector_set( f , count , ( exp_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;

  return GSL_SUCCESS;
}

int
exp_df( const gsl_vector *__restrict x, 
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

      int k ;
      // initialise whole column to zero
      for( k = 0 ; k < DATA->LOGICAL_NPARS ; k++ ) {
	gsl_matrix_set ( J , count , k ,  0.0 ) ;
      }

      // pretty hideous extension here
      double e ;
      if( DATA->sim_params[0] == true ) {
	const double _s = 1.0 / DATA->sigma[i] ;
	e = exp( -params[0] * DATA->X[i] ) * _s ;
	gsl_matrix_set ( J , count , 0 , -DATA->X[i] * params[loc_idx+1] * e ) ;
	gsl_matrix_set ( J , count , loc_idx+1 ,  e ) ;
	//////////////////////////////////////////////////////////////////////////
      } else if( DATA->sim_params[1] == true ) {
	// this is wrong -> FIX!!!!!
	if( NSIM == 0 ) {
	  const double _s = 1.0 / DATA->sigma[i] ;
	  e = exp( -params[0] * DATA->X[i] ) * _s ;
	  gsl_matrix_set ( J , count , 0 , -DATA->X[i] * params[1] * e ) ;
	  gsl_matrix_set ( J , count , 1 ,  e  ) ;
	} else {
	  const double _s = 1.0 / DATA->sigma[i] ;
	  e = exp( -params[loc_idx+1] * DATA->X[i] ) * _s ;
	  //////////////////////////////////////////////////////////////////////////
	  gsl_matrix_set ( J , count , loc_idx+1 , -DATA->X[i] * params[1] * e ) ;
	  gsl_matrix_set ( J , count , 1 ,  e  ) ;
	}

      } else {

	const double _s = 1.0 / DATA->sigma[i] ;
	e = exp( -params[loc_idx] * DATA->X[i] ) * _s ;
	gsl_matrix_set ( J , count , loc_idx+0 , -DATA->X[i] * params[loc_idx+1] * e ) ;
	gsl_matrix_set ( J , count , loc_idx+1 ,  e  ) ;
      }
      
      count ++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  free( (void*)params ) ;

  return GSL_SUCCESS;
}

int
exp_fdf(const gsl_vector *__restrict x , 
	void *data ,
	gsl_vector *__restrict f , 
	gsl_matrix *__restrict J )
{
  exp_f( x , data , f ) ;
  exp_df( x , data , J ) ;

  return GSL_SUCCESS;
}
