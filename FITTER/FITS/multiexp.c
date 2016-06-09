/**
   Multiple exponential fit, solves the system

   | 1    1    1    .. 1     | | A_0 |   | Y(0)  |
   | a_1  a_2  a_3  .. a_m   | |  .  | = | Y(1)  |
   | a_1^2 ............a_m^2 | |  .  | = | .     |
   | a_1^n .............     | | A_m |   | Y(LT) |

   For the parameters a_i and A_i s.t.

   Y(t) = \sum_{i} A_i * a_i^t

   For exponentially-decaying correlation functions
   
   a_i = exp( -M_i )

   Maybe this should be done explicitly?

   Also, M < N and can be anything we want

   I pack all of these into the fit parameters, there will be 2M of them

   the first M will be the a_i's and the second M will be the A_i's

   1)
      I force the parameters to be positive so you will need to take their square root
      to obtain the masses

   2)
      I "order" the masses, so the first is e^{-M_0^2 t} the second is e^{-( M_0^2 + M_1^2 )*t } and so on
      this will need to be addressed too
 */

#include "fitfunc.h"
#include "Utils.h"
#include "blackbox.h"

const static int M = 2 ;

int getm( void )
{
  return 2*M ;
}

// provide a description of the function
void
mexp_description( const char *message ,
		 const double *__restrict params ,
		 const struct x_descriptor X ,
		 const int NPARAMS )
{
  int m ;
  for( m = 0 ; m < M ; m++ ) {
    printf( "%s :: y = %f * ( mexp( -%f * x ) ) \n" ,
	    message , params[1 + 2*m] , params[0+2*m] ) ;
  }
  return ;
}

// I provide a guess too
void
mexp_guess( double *__restrict params ,
	   void *data ,
	   const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;

  int k , xposit = 0 ;
  for( k = 0 ; k < DATA->SIMS ; k++ ) {
    double tempdata[ DATA -> NDATA[k] ] ;

    int i ;
    for( i = 0 ; i < DATA -> NDATA[k] ; i++ ) {
      tempdata[ i ] = DATA -> y[ xposit + i ] ;
    }

    double masses[ M ][ DATA -> NDATA[k] ] ;
    blackbox( tempdata , DATA -> NDATA[k] , M , masses ) ;

    int m ;
    for( m = 0 ; m < M ; m++ ) {
      params[0+2*m] = masses[m][0] ;
      params[1+2*m] = (tempdata[1]+tempdata[2]) /			\
	(exp( -masses[m][0]*DATA->X[xposit+1] )+exp( -masses[m][0]*DATA->X[xposit+2] )) ;
    }
    struct x_descriptor X ;
    mexp_description( "Multiexp" , params , X , NPARAMS ) ;
    xposit += DATA->NDATA[k] ;
  }

  return ;
}

// evaluate the function at the point "X" using the fitparams
inline double
mexp_eval( const double *__restrict x ,
	  const struct x_descriptor X ,
	  const int NPARAMS )
{
  int i ;
  register double sum = 0.0 ;
  for( i = 0 ; i < M ; i++ ) {
    sum += x[1+2*i] * exp( -x[2*i] * X.X ) ;
    if( fabs( x[2*i] ) > 10.0 ) {
      sum *= 100 ;
    }
  }
  if( x[0] > 0 ) {
    sum += 100 ;
  }
  return sum ;
}

// compute the difference of the function from the data
int
mexp_f( const gsl_vector *__restrict x , 
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
    size_t i ;
    for( i = 0 ; i < M ; i++ ) {
      if( DATA->sim_params[ i ] == true ) {
	p[ 2*i ] = params[ 2*i ] ;
	p[ 2*i + 1 ] = params[ loc_idx + 2*i + 1 ] ;
      } else {
	p[ 2*i + 0 ] = params[ loc_idx + 2*i + 0 ] ;
	p[ 2*i + 1 ] = params[ loc_idx + 2*i + 1 ] ;
      }
    }

    for( i = lopos ; i <= hipos ; i++ ) {
      const struct x_descriptor XAXIS = { DATA->X[i] } ;
      gsl_vector_set( f , count , ( mexp_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;

  return GSL_SUCCESS;
}

int
mexp_df( const gsl_vector *__restrict x, 
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
      const double _s = 1.0 / DATA->sigma[i] ;
      for( k = 0 ; k < M ; k++ ) {
	if( DATA->sim_params[ k ] == true ) {
	  double e = exp( -params[2*k] * DATA->X[i] ) * _s ;
	  gsl_matrix_set ( J , count , 2*k , -DATA->X[i] * params[loc_idx+2*k+1] * e ) ;
	  gsl_matrix_set ( J , count , loc_idx+1+2*k ,  e  ) ;
	} else {
	  double e = exp( -params[loc_idx+2*k] * DATA->X[i] ) * _s ;
	  gsl_matrix_set ( J , count , loc_idx+2*k , -DATA->X[i] * params[loc_idx+2*k+1] * e ) ;
	  gsl_matrix_set ( J , count , loc_idx+1+2*k ,  e  ) ;
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
mexp_fdf(const gsl_vector *__restrict x , 
	void *data ,
	gsl_vector *__restrict f , 
	gsl_matrix *__restrict J )
{
  mexp_f( x , data , f ) ;
  mexp_df( x , data , J ) ;

  return GSL_SUCCESS;
}
