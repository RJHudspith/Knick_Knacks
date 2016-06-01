/**
   @file adaptive_simpson.c
   @brief adaptive simpson's integration
 */
#include <stdio.h>
#include <math.h>

// simpson evaluation
static inline double
simp_eval( double (*f)( const double ) , const double x , const double h )
{
  return h * ( f( x ) + 4. * f( x + h/2. ) + f( x + h ) ) ;
}

// classic simpson's
static double
simpsons( double (*f)( const double ) ,
	  const double low , 
	  const double upp , 
	  const int nsteps )
{
  register double sum = 0.0 ;
  const double h = ( upp - low ) / (double)nsteps ;
  int i ;
  for( i = 0 ; i < nsteps ; i++ ) {
    sum += simp_eval( f , low + i * h , h ) ;
  }
  return sum / 6. ;
}

// driver for the adaptive simpson's
// usual two half-step approach
static double
adaptive_simpsons( double (*f)( const double ) ,
 		   const double low , 
		   const double upp ,
		   const double eps )
{
  const double inveps = 1.0 / eps ; 
  double fx , fxph_4 , fxph_2 , fxp3h_4 , fxph ; // some evaluations I use a lot
  double h = ( upp - low ) / 100. ; // initial guess
  register double sum = 0.0 ;
  double x = low ;
  while( x < upp ) {
    double diff , step = 0.0 , step2 = 0.0 ;
    int nsteps = 0 ;
    while( nsteps < 50 ) {

      fx      = f( x ) ;
      fxph_4  = f( x + h/4. ) ;
      fxph_2  = f( x + h/2. ) ;
      fxp3h_4 = f( x + 3.*h/4. ) ;
      fxph    = f( x + h ) ;

      step  = h * ( fx + 4. * fxph_2 + fxph ) ;
      step2 = h * ( fx + 4. * ( fxph_4 + fxp3h_4 ) + 2. * fxph_2 + fxph ) / 2.0 ;
 
      diff = fabs( step - step2 ) * inveps ;
 
      if( diff < 1.0 ) 
	break ;
      if( nsteps == 49 ) {
	printf( "Adaptive not converging !!! Leaving \n" ) ;
	exit(1) ;
      }
      h *= 0.9 * pow( diff , -0.25 ) ; 
      nsteps++ ;
    }
    sum = sum + ( step2 + ( step2 - step ) / 15. ) ;
    x += h ;
    h *= 0.9 * pow( diff , -0.2 ) ;
  }
  // small correction exactly to "x"
  return ( sum + simp_eval( f , x , upp-x ) ) / 6. ;
}

// mainfunc
int main( void )
{
  printf( "%e \n" , 1.0 - adaptive_simpsons( sin , 0.0 , M_PI/2. , 1E-14 ) ) ;
  printf( "%e \n" , 1.0 - simpsons( sin , 0.0 , M_PI/2. , 1000 ) ) ;
  return 0 ;
}
