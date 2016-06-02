/**
   @file pp_ap_pa_aa.c
   @brief simultaneous fit over pseudoscalar channels
 */
#include "fitfunc.h"
#include "Utils.h"

// cute little enum
enum { VWVL , VWVW , VTWVTL , VTWVTW } ;

// provide a description of the function
void
VV_VTVT_WW_description( const char *message ,
		  const double *__restrict params ,
		  const struct x_descriptor X ,
		  const int NPARAMS )
{
  printf( "In guess \n" ) ;
  printf( "%s :: y = %f * ( exp( -%f * x ) + exp( -%f * ( %d - x - 1 ) ) ) \n" ,
	  message , params[1] , params[0] , params[0] , X.LT ) ;
  printf( "In guess \n" ) ;
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
  int i = 0 , N = 0 , negs = 0 ;
  double effmass_prev = -log( fabs( y[i+1] / y[i] ) ) / ( X[i+1] - X[i] ) ;
  double effmass_this = 0.0 ;
  
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
  if( negs > NDATA/4 ) *amp *= -1 ;
  return true ;
}

// I provide a guess too
void
VV_VTVT_WW_guess( double *__restrict params ,
	       void *data ,
	       const int NPARAMS )
{
  const struct data* DATA = (const struct data*)data ;

  //double summass = 0.0 ;
  int k = 0 , xposit = 0 ;
  for( k = 0 ; k < 1 ; k++ ) {

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
		     0.075 // this is a fudge parameter :: estimates how constant our meff is
			 ) == false ) {
      p[ 0 ] = params[0] ;
      p[ 1 ] = ( ydata[1] + ydata[2] ) / 
	( exp( -params[0] * xdata[1] ) + exp( -params[0] * xdata[2] ) ) ;
    }

    params[0] = p[0] ;

    struct mom_info momdummy ;
    struct chiral qdummy ;
    struct x_descriptor X = { 0 , qdummy , momdummy , DATA->LT } ;
    VV_VTVT_WW_description( "VV_VTVT_WW guess" , p , X , 2 ) ;
  }

  params[1] = 15 ;
  params[2] = 2400 ;
  params[3] = 15 ;
  params[4] = 2400 ;
  return ;
}

// evaluate the function at the point "X" using the fitparams
inline double
VV_VTVT_WW_eval( const double *__restrict x ,
	       const struct x_descriptor X ,
	       const size_t chan )
{
  switch( chan ) {
    // WALL-LOCAL
  case VWVL :
    return x[2] * x[1] * ( exp( -x[0] * X.X ) + exp( -x[0] * ( X.LT - X.X ) ) ) ;
  case VWVW :
    return x[2] * x[2] * ( exp( -x[0] * X.X ) + exp( -x[0] * ( X.LT - X.X ) ) ) ;
  case VTWVTL :
    return x[4] * x[3] * ( exp( -x[0] * X.X ) + exp( -x[0] * ( X.LT - X.X ) ) ) ;
  case VTWVTW :
    return x[4] * x[4] * ( exp( -x[0] * X.X ) + exp( -x[0] * ( X.LT - X.X ) ) ) ;
  }
  return 0 ;
}

// compute the difference of the function from the data
int
VV_VTVT_WW_f( const gsl_vector *__restrict x , 
	void *data , 
	gsl_vector *__restrict f )
{
  const struct data *DATA = (const struct data *)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  struct chiral dummy_chiral ;
  struct mom_info dummy_mominfo ;

  int NSIM , count = 0 , xposit = 0 ;
  for( NSIM = 0 ; NSIM < DATA->SIMS ; NSIM++ ) {

    // lo and hi for this x-range
    const int lopos = find_idx( DATA->FIT_LO , DATA->X , 
				xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , 
				xposit + DATA->NDATA[NSIM] , xposit ) ;

    size_t i ;
    for( i = lopos ; i <= hipos ; i++ ) {
      const struct x_descriptor XAXIS = { DATA->X[i] , dummy_chiral , dummy_mominfo , DATA -> LT } ;
      gsl_vector_set ( f , count , 
		       ( VV_VTVT_WW_eval( params , XAXIS , NSIM ) - DATA->y[i] ) / DATA->sigma[i] ) ;
      count++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }
  free( (void*)params ) ;

  return GSL_SUCCESS;
}

int
VV_VTVT_WW_df( const gsl_vector *__restrict x, 
	 void *data, 
	 gsl_matrix *__restrict J )
{
  const struct data *DATA = (const struct data*)data ;
  const double *params = set_params( x , DATA->LOGICAL_NPARS ) ;

  int count = 0 , NSIM , xposit = 0 ; 
  for( NSIM = 0 ; NSIM < DATA-> SIMS ; NSIM++ ) {

    const int lopos = find_idx( DATA->FIT_LO , DATA->X , 
				xposit + DATA->NDATA[NSIM] , xposit ) ;
    const int hipos = find_idx( DATA->FIT_HI , DATA->X , 
				xposit + DATA->NDATA[NSIM] , xposit ) ;

    size_t i ; 
    for( i = lopos ; i <= hipos ; i++ ) {
      const double _s = 1.0 / DATA->sigma[i] ;

      // initialise whole column to zero
      int k ;
      for( k = 0 ; k < DATA->LOGICAL_NPARS ; k++ ) {
	gsl_matrix_set (J, count , k ,  0.0 ) ;
      }

      // precompute forward and backward props
      const double fw = DATA->X[i] ;
      const double bw = DATA->LT - DATA->X[i] ;
      double e1 = 0.0 , e2 = 0.0 ;

      e1 = exp( -params[0] * fw ) * _s ;
      e2 = exp( -params[0] * bw ) * _s ;

      // and the shared ones
      switch( NSIM%4 ) {
      case VWVL :
	gsl_matrix_set ( J , count , 0 , -params[2] * params[1] * ( fw * e1 + bw * e2 ) ) ;
	gsl_matrix_set ( J , count , 1 ,  params[2] * ( e1 + e2 ) ) ;
	gsl_matrix_set ( J , count , 2 ,  params[1] * ( e1 + e2 ) ) ;
	break ;
       case VWVW :
	gsl_matrix_set ( J , count , 0 , -params[2] * params[2] * ( fw * e1 + bw * e2 ) ) ;
	gsl_matrix_set ( J , count , 2 ,  2 * params[2] * ( e1 + e2 ) ) ;
	break ;
      case VTWVTL :
	gsl_matrix_set ( J , count , 0 , -params[4] * params[3] * ( fw * e1 + bw * e2 ) ) ;
	gsl_matrix_set ( J , count , 3 ,  params[4] * ( e1 + e2 ) ) ;
	gsl_matrix_set ( J , count , 4 ,  params[3] * ( e1 + e2 ) ) ;
	break ;
       case VTWVTW :
	gsl_matrix_set ( J , count , 0 , -params[4] * params[4] * ( fw * e1 + bw * e2 ) ) ;
	gsl_matrix_set ( J , count , 4 ,  2 * params[4] * ( e1 + e2 ) ) ;
	break ;
      }
      count ++ ;
    }
    xposit += DATA->NDATA[NSIM] ;
  }

  free( (void*)params ) ;

  return GSL_SUCCESS;
}

int
VV_VTVT_WW_fdf(const gsl_vector *__restrict x , 
	 void *data ,
	 gsl_vector *__restrict f , 
	 gsl_matrix *__restrict J )
{
  VV_VTVT_WW_f (x, data, f);
  VV_VTVT_WW_df (x, data, J);
  return GSL_SUCCESS;
}
