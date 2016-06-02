/**
   @file integrators.c
   @brief numerical integrators
 */
#include <stdio.h>
#include <math.h>

// simpson evaluation
static inline double
simp_eval( double (*f)( const double ) , const double x , const double h )
{
  return h * ( f( x ) + 4. * f( x + h/2. ) + f( x + h ) ) ;
}

// driver for the adaptive simpson's
// usual two half-step approach to estimate the error
double
adaptive_simpsons( double (*f)( const double ) ,
 		   const double low , 
		   const double upp ,
		   const double eps )
{
  const double inveps = 1.0 / eps ; 
  register double fx , fxph_4 , fxph_2 , fxp3h_4 , fxph ; // some evaluations I use a lot
  register double sum = 0.0 , diff ;
  double h = ( upp - low ) / 100. ; // initial guess
  double x = low ;
  const size_t max_steps = 150 ;
  while( x < upp ) {
    register double step = 0.0 , step2 = 0.0 ;
    int nsteps = 0 ;
    while( nsteps < max_steps ) {

      fx      = f( x ) ;
      fxph_4  = f( x + h/4. ) ;
      fxph_2  = f( x + h/2. ) ;
      fxp3h_4 = f( x + 3 * h/4. ) ;
      fxph    = f( x + h ) ;

      step  = h * ( fx + 4. * fxph_2 + fxph ) ;
      step2 = h * ( fx + 4. * ( fxph_4 + fxp3h_4 ) + 2. * fxph_2 + fxph ) / 2.0 ;
 
      diff = fabs( step - step2 ) * inveps ;
 
      // break if the diff is within tolerance
      if( diff < 1.0 ) 
	break ;

      // if we run into too much trouble we return a NAN
      if( nsteps == max_steps-1 ) {
	printf( "[ADS] ill convergence met\n" ) ;
	printf( "[ADS] %e \n" , diff ) ;
	return ( sum + simp_eval( f , x , upp-x ) ) / 6.0 ;
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

static double
general_simpsons_array_O3( const double *x , 
			   const double *y , 
			   const int i )
{
  const register double a = x[ i ] ;
  const register double b = x[ i + 1 ] ;
  const register double c = x[ i + 2 ] ;
  const register double fa = y[ i ] / 6.0 ;
  const register double fb = y[ i + 1 ] / 6.0 ;
  const register double fc = y[ i + 2 ] / 6.0 ;

  // up to some tolerance we just return the trapezoid evaluation
  if( ( b - a ) < 1E-20 ) {
    return ( fc - fb ) * ( c - b ) * 0.5 ;
  } else if( ( c - b ) < 1E-20 ) {
    return ( fb - fa ) * ( b - a ) * 0.5 ;
  } else {
    register double loc_sum = a * a * ( fb - fc ) ;
    loc_sum += c * c * ( fb - fa ) ;
    loc_sum += -3 * b * b * ( fa + fc ) ;
    loc_sum += 2 * b * c * ( 2 * fa + fc ) ;
    loc_sum += 2 * a * b * ( fa + 2 * fc ) ;
    loc_sum += -2 * a * c * ( fa + fb + fc ) ;
    loc_sum *= ( c - a ) / ( ( a - b ) * ( b - c ) ) ;
    return loc_sum ;
  }
}

static double
general_simpsons_array_O4( const double *x , 
			   const double *y , 
			   const int i )
{
  const double a = x[ i ] ;
  const double b = x[ i + 1 ] ;
  const double c = x[ i + 2 ] ;
  const double d = x[ i + 3 ] ;

  const double fa = y[ i ] ;
  const double fb = y[ i + 1 ] ;
  const double fc = y[ i + 2 ] ;
  const double fd = y[ i + 3 ] ;

  long double loc_sumfa , loc_sumfb , loc_sumfc , loc_sumfd ;

  // f(a) multiplier
  loc_sumfa  = 3 * a * a + 6 * b * c - 2 * ( b + c ) * d ;
  loc_sumfa += d * d + 2 * a * ( -2 * ( b + c ) + d ) ;
  loc_sumfa *= -fa / ( ( a - b ) * ( a - c ) ) ;

  // f(b) multiplier
  loc_sumfb  = ( a - d ) * ( a - d ) * ( a - 2. * c + d ) ;
  loc_sumfb *= fb / ( ( b - a ) * ( b - c ) * ( b - d ) ) ;
  
  // f(c) multiplier
  loc_sumfc  = ( a - d ) * ( a - d ) * ( a - 2. * b + d ) ;
  loc_sumfc *= fc / ( ( c - a ) * ( c - b ) * ( c - d ) ) ;

  // f(d) multiplier
  loc_sumfd  = a * a + 6 * b * c - 2 * a * ( b + c - d ) ;
  loc_sumfd += -4. * b * d - 4. * c * d + 3 * d * d ;
  loc_sumfd *= fd / ( ( b - d ) * ( d - c ) ) ;

  return ( loc_sumfa + loc_sumfb + loc_sumfc + loc_sumfd ) * ( a - d ) / 12. ;
}

// numerical integrator
double
simpsons_arr5( const double *y , 
	       const double *x ,
	       const int N )
{
  register double sum = 0.0 ;

  const double lo = general_simpsons_array_O4( x , y , 0 )
    - general_simpsons_array_O3( x , y , 1 ) ;

  int i ;
  for( i = 0 ; i < N - 2 ; i++ ) {
    sum += general_simpsons_array_O3( x , y , i ) ;
  }

  const double hi = general_simpsons_array_O4( x , y , N-4 )
    - general_simpsons_array_O3( x , y , N - 4 ) ;

  return 0.5 * ( lo + sum + hi ) ;
}
