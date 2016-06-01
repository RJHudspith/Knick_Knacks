/**
   @file generic_spacing.c
   @brief integrators for data evenly and unevenly spaced
 */
#include <stdio.h>

// simple simpson's rule integrator
// (upp-low)/( 6 * N ) * sum_x ( f( x ) + 4*f( x + h/2 ) + f( x + h ) )
static double
simpsons2( double (*f)( const double ) ,
	   const double low , 
	   const double upp , 
	   const int nsteps )
{
  register double sum = 0.0 ;
  const double h = ( upp - low ) / (double)nsteps ;
  int i ;
  for( i = 0 ; i < nsteps ; i++ ) {
    const double x = low + i * h ;
    sum += ( f( x ) + 2. * f( x + h/2. ) ) ;
  }
  sum += 0.5 * ( f( upp ) - f( low ) ) ;
  return h * sum / 3. ;
}

// trapezoid rule integration
// 0.5 * sum_{i} ( x[i+1] - x[i] ) * [ f( x[i] ) + f( x[i+1] ) ]
// corresponds to y[i] = f( x[i] ) 
// where x must be sorted in ascending order
static double
nintegrate_arr( const double *y , 
		const double *x ,
		const int LEN )
{
  register double sum = 0.0 ;
  int i ;
  for( i = 0 ; i < LEN-1 ; i++ ) {
    const double h = ( x[i+1] - x[i] ) ;
    sum += h * ( y[i] + y[i+1] ) ;
  }
  return 0.5 * sum ;
}

// general simpson's rule integration
// corresponds to y[i] = f( x[i] ) 
// where x must be sorted in ascending order
static double
simpsons_arr( const double *y , 
	      const double *x ,
	      const int N )
{
  register double sum = 0.0 ;
  register double loc_sum ;
  int i ;
  for( i = 0 ; i < N - 2 ; i+=2 ) {

    const double a = x[ i ] ;
    const double b = x[ i + 1 ] ;
    const double c = x[ i + 2 ] ;

    const double fa = y[ i ] ;
    const double fb = y[ i + 1 ] ;
    const double fc = y[ i + 2 ] ;

    loc_sum = 0.0 ;

    loc_sum += a * a * ( fb - fc ) ;
    loc_sum += c * c * ( fb - fa ) ;
    loc_sum += -3 * b * b * ( fa + fc ) ;
    loc_sum += 2 * b * c * ( 2 * fa + fc ) ;
    loc_sum += 2 * a * b * ( fa + 2 * fc ) ;
    loc_sum += -2 * a * c * ( fa + fb + fc ) ;
    loc_sum *= ( c - a ) / ( ( a - b ) * ( b - c ) ) ;

    sum += loc_sum ;
  }
  return ( N&1 != 0 )? sum/6. : sum/6. + 0.5 * ( x[i+1] - x[i] ) * ( y[i] + y[i+1] ) ;
}

// generic (separation agnostic) trapezoid method
static inline double
general_trapezoid( const double *x , const double *y , const int i ) 
{
  return -0.5 * ( x[i] - x[i+1]) * ( y[i] + y[i+1] ) ;
}

// third order spacing agnostic simpsons method
static inline double
general_simpsons_array_O3( const double *x , const double *y , const int i )
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

// fourth order generic spacing simpson's method
static inline double
general_simpsons_array_O4( const double *x , const double *y , const int i )
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


// general simpson's rule integration
// corresponds to y[i] = f( x[i] ) 
// where x must be sorted in ascending order
static double
simpsons_arr2( const double *y , 
	       const double *x ,
	       const int N )
{
  register double sum = 0.0 ;
  const double lo = general_trapezoid( x , y , 0 ) ;

  int i ;
  for( i = 0 ; i < N - 2 ; i++ ) {
    sum += general_simpsons_array_O3( x , y , i ) ;
  }
  const double hi = general_trapezoid( x , y , N-2 ) ;

  return 0.5 * ( lo + sum + hi ) ;
}

// fourth order with trapezoids on either end
static double
simpsons_arr3( const double *y , 
	       const double *x ,
	       const int N )
{
  register double sum = 0.0 ;

  int i ;
  for( i = 0 ; i < N - 3 ; i++ ) {
    sum += general_simpsons_array_O4( x , y , i ) ;
  }
    
  return ( + 2.0 * general_trapezoid( x , y , 0 ) 
	   + general_trapezoid( x , y , 1 ) 
	   + sum 
	   + general_trapezoid( x , y , N - 3 ) 
	   + 2.0 * general_trapezoid( x , y , N - 2 ) ) / 3.0 ;

}

// trapezoids fill in the spaces for the simpson's on this one
static double
simpsons_arr4( const double *y , 
	       const double *x ,
	       const int N )
{
  register double sum = 0.0 ;

  int i ;
  for( i = 0 ; i < N - 3 ; i+=3 ) {
    sum += general_simpsons_array_O4( x , y , i ) ;
  }
  
  //printf( "%f \n" , sum ) ;
  switch( (N-1)%3 ) {
  case 0 :
    return sum ;
  case 1 :
    return sum + general_trapezoid( x , y , N - 2 ) ;
  case 2 :
    return sum + general_trapezoid( x , y , N - 3 ) + general_trapezoid( x , y , N - 2 ) ;
  }
}

// computes a correction at either end and then O(3) simpson's in the middle
static double
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

// mainfunc
int main( void )
{
  return 0 ;
}
