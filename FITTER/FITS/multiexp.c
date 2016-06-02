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

const static int M = 2 ;

//#define STANDARD
//#define SQUARED
#define SUM
//#define SQUARED_SUM

// get the number of states
int
getm( )
{
  return 2 * M ;
}

// provide a description of the function
void
mexp_description( const char *message ,
		 const double *__restrict params ,
		 const struct x_descriptor X ,
		 const int NPARAMS )
{
  printf( "%s multi-exponential fit\n" , message ) ;
  char str[ 256 ] = {} ;
  int i ;
  for( i = 0 ; i < M ; i++ ) {
    sprintf( str , "%s%f^2" , str , params[i] ) ;
    printf( "(%f)^2 * exp{ -( %s ) } \n" , params[M+i] , str ) ;
    sprintf( str , "%s + " , str ) ;
  }

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
  double effmass_prev = -log( y[i+1] / y[i] ) / ( X[i+1] - X[i] ) ;
  double effmass_this = 0.0 ;

  double sumx = 0.0 , sumy = 0.0 , sumxy = 0.0 , sumxx = 0.0 ;
  for( i = 1 ; i < NDATA ; i++ ) {
    if( y[i] < 0.0 ) break ;
    effmass_this = -log( fabs( y[i+1] / y[i] ) ) / ( X[i+1] - X[i] ) ;
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
mexp_guess( double *__restrict params ,
	    void *data ,
	    const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;
  double tempdata[ DATA -> NDATA[0] ] ;

  int i ;
  for( i = 0 ; i < DATA -> NDATA[0] ; i++ ) {
    tempdata[ i ] = DATA -> y[ i ] ;
  }

  double summass = 0.0 , summasssq = 0.0 , sumamp = 0.0 ;
  int k ;
  for( k = 0 ; k < M ; k++ ) {

    // regression using logs to get the mass and amplitude
    if( effmass_regress( &params[k] , &params[k+M] ,
			 tempdata , 
			 DATA->X ,
			 DATA->NDATA[0] ,
			 0.075 + 0.66 * k // this is a fudge parameter :: estimates how constant our meff is
			 ) == false ) {
      if( k > 0 ) {
	params[k] = params[k-1] + 0.5 ;
	params[k+M] = params[k+M-1] + 100 ;
      } else {
	params[k] = 1.0 ;
	params[k+M] = 100. ;
      }
    } else {
      // correct the params
      #if (defined STANDARD)
      params[ k ] = fabs( params[k] ) ;
      #elif (defined SQUARED)
      params[ k ] = sqrt( fabs( params[k] ) ) ;
      params[k+M] = sqrt( fabs( params[k+M] ) ) ;
      #elif (defined SUM)
      params[ k ] = fabs( params[k] ) - summass ;
      #elif (defined SQUARED_SUM)
      params[ k ] = sqrt( fabs( fabs( params[k] ) - summasssq ) ) ;
      params[k+M] = sqrt( fabs( params[k+M] ) ) ;
      #endif
    }

    summass += params[k] ;
    summasssq += ( params[k] * params[k] ) ;
    #if (defined SQUARED) || (defined SQUARED_SUM)
    sumamp += params[k+M] * params[k+M] ;
    #else
    sumamp += params[k+M] ;
    #endif

    for( i = 0 ; i < DATA -> NDATA[0] ; i++ ) {
      #if defined STANDARD
      tempdata[ i ] -= ( params[k+M] * exp( -( params[k] ) * DATA->X[i] ) ) ;
      #elif defined SQUARED
      tempdata[ i ] -= ( ( params[k+M] * params[k+M] ) * exp( -( params[k] * params[k] ) * DATA->X[i] ) ) ;
      #elif defined SUM
      tempdata[ i ] -= ( params[k+M] * exp( -( summass ) * DATA->X[i] ) ) ;
      #elif defined SQUARED_SUM
      tempdata[ i ] -= ( ( params[k+M] * params[k+M] ) * exp( -( summasssq ) * DATA->X[i] ) ) ;
      #endif
    }
  }

  // sort the guesses?

  struct x_descriptor X ;
  mexp_description( "Guess" , params , X , NPARAMS ) ;

  return ;
}

// evaluate the function
inline double
mexp_eval( const double *__restrict p ,
	  const struct x_descriptor X ,
	  const int NPARAMS )
{
  // is a big matrix multiply
  register double sum = 0.0 , mass_sum = 0.0 ;
  int i ;
  for( i = 0 ; i < M ; i++ ) {
    #ifdef STANDARD
    sum += exp( -p[ i ] * X.X ) * p[ i + M ] ;
    #elif defined SQUARED
    sum += exp( -( p[ i ] * p[ i ] ) * X.X ) * p[ i + M ] * p[ i + M ] ;
    #elif defined SUM
    mass_sum += ( p[ i ] ) ;
    sum += exp( -mass_sum * X.X ) * ( p[ i + M ] ) ;
    #elif defined SQUARED_SUM
    mass_sum += ( p[ i ] * p[ i ] ) ;
    sum += exp( -mass_sum * X.X ) * ( p[ i + M ] * p[ i + M ] ) ;
    #endif
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
  const double *params = set_params( x , 2*M ) ;
  const struct mom_info dummy_mom ;

  int NSIM , count = 0 , xposit = 0 ;
  for( NSIM = 0 ; NSIM < DATA->SIMS ; NSIM++ ) {

    // local idx gives the unshared fit parameters
    const int loc_idx = NSIM*(DATA->NPARAMS-DATA->NCOMMON) ;

    const int lopos = find_idx( DATA->FIT_LO , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , xposit + DATA->NDATA[NSIM] , xposit ) ;

    // local parameters ...
    double p[ 2*M ] ;
    size_t i ;
    for( i = 0 ; i < 2*M ; i++ ) {
      p[ i ] = params[ loc_idx + i ] ;
    }

    for( i = lopos ; i <= hipos ; i++ ) {
      const struct x_descriptor XAXIS = { DATA->X[i] , DATA->quarks[0] , dummy_mom , DATA->NDATA[NSIM] } ;
      gsl_vector_set ( f , count , 
		       ( mexp_eval( p , XAXIS , DATA->NPARAMS ) - DATA->y[i] ) / DATA->sigma[i] ) ;
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
  const double *params = set_params( x , 2*M ) ;

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
      for( k = 0 ; k < 2*M ; k++ ) {
	gsl_matrix_set ( J , count , k+loc_idx , 0.0 ) ;
      }


#if (defined SUM) || (defined SQUARED_SUM)
      double mass_sum = 0.0 , sumexp = 0.0 ;
      for( k = 0 ; k < M ; k++ ) {
	const double mass = params[k+loc_idx] ;
	const double amp = params[k+M+loc_idx] ;
	#if defined SUM
	mass_sum += mass ;
	sumexp += amp * exp( -( mass_sum ) * DATA->X[i] ) * _s ;
	#elif defined SQUARED_SUM
	mass_sum += mass * mass ;
	sumexp += ( amp * amp ) * exp( -( mass_sum ) * DATA->X[i] ) * _s ;
	#endif
      }
      mass_sum = 0.0 ;
      for( k = 0 ; k < M ; k++ ) {
	const double mass = params[k+loc_idx] ;
	const double amp = params[k+M+loc_idx] ;
	#if defined SUM
	mass_sum += mass ;
	const double e = exp( -( mass_sum ) * DATA->X[i] ) * _s ;
	gsl_matrix_set ( J , count , k+loc_idx , -DATA->X[i] * sumexp ) ;
	gsl_matrix_set ( J , count , k+M+loc_idx , e ) ;
	sumexp -= ( amp * e ) ;
        #elif defined SQUARED_SUM
	mass_sum += mass * mass ;
	const double e = exp( -( mass_sum ) * DATA->X[i] ) * _s ;
	gsl_matrix_set ( J , count , k+loc_idx , -2.0 * DATA->X[i] * mass * sumexp ) ;
	gsl_matrix_set ( J , count , k+M+loc_idx , 2.0 * amp * e ) ;
	sumexp -= ( amp * amp * e ) ;
	#endif
      }
      //printf( "%e \n" , sumexp ) ;
#else
      for( k = 0 ; k < M ; k++ ) {
	const double mass = params[k+loc_idx] ;
	const double amp = params[k+M+loc_idx] ;
	#ifdef STANDARD
	const double e = exp( -params[k+loc_idx] * DATA->X[i] ) * _s ;
	gsl_matrix_set ( J , count , k+loc_idx , -DATA->X[i] * params[k+M+loc_idx] * e ) ;
	gsl_matrix_set ( J , count , k+M+loc_idx , e ) ;
	#elif defined SQUARED
	// squared parameters
	const double e = exp( -( mass * mass ) * DATA->X[i] ) * _s ;
	gsl_matrix_set ( J , count , k+loc_idx , -DATA->X[i] * 2.0 * amp * amp * mass * e ) ;
	gsl_matrix_set ( J , count , k+M+loc_idx , 2.0 * amp * e ) ;
	#endif
      }
#endif

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
