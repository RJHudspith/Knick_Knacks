/**
   Reads in data of length "N" and after interpolation to a point
   packs the data into the resampled form ...
 */
#include "fitfunc.h"

static double
sum_residuals( const double *x , 
	       const double *y , 
	       const int NDATA ,
	       const double slope ,
	       const double intercept )
{
  int i ;
  double sumres = 0. ;
  for( i = 0 ; i < NDATA ; i++ ) {
    const double res = y[i] - ( slope * x[i] + intercept ) ;
    sumres += res * res ;      
  }
  return sumres ;
}

static double
simple_least_squares( const double *x ,
		      const double *y ,
		      const int NDATA ,
		      const double extrap ,
		      double *x_extrap ) 
{
  double sumx = 0. , sumy = 0. , sumxy = 0. , sumxx = 0. ;
  int i ;
  for( i = 0 ; i < NDATA ; i++ ) {
    sumx += x[i] ;
    sumy += y[i] ;
    sumxy += x[i] * y[i] ;
    sumxx += x[i] * x[i] ;
  }

  // compute the slope and the intercept
  const double slope = ( sumx * sumy - NDATA * sumxy ) / ( sumx * sumx - NDATA * sumxx ) ;
  const double y_intercept = ( sumy - slope * sumx) / ( NDATA ) ;

  // pass by reference the x-value
  *x_extrap = ( extrap - y_intercept ) / slope ;

  // tell us what the equation is
  printf( "Equation is :: Y = %f + x * %f \n" , y_intercept , slope ) ;

  // extrapolation point and residuals
  printf( "Extrap %f :: Residuals :: %e \n" , *x_extrap , 
	  sum_residuals( x , y , NDATA , slope , y_intercept ) ) ;

  return slope ;
}

// computes a distribution of t_0, which can be resampled and fit
struct resampled
obtain_t0( const char *filename ,  // the filename, expects a %d number at the end of it
	   const double VAL ,      // the point to measure t_0 (0.3)
	   const int TRAJ_START , 
	   const int TRAJ_END ,
	   const int TRAJ_INC ,    
	   const int NDATA ,       // number of data points
	   const bool DER           // are we taking the derivative?
	   )
{
  const int N = ( TRAJ_END - TRAJ_START ) / TRAJ_INC ;
  struct resampled T0 ;
  T0.resampled = malloc( N * sizeof( struct resampled ) ) ;
  T0.NSAMPLES = N ;

  int i ;
#pragma omp parallel for private(i)
  for( i = TRAJ_START ; i < TRAJ_END ; i+= TRAJ_INC ) {
    
    const int meas = ( i - TRAJ_START ) / TRAJ_INC ;

    // open files
    char str[ 128 ] ;
    sprintf( str , "%s.%d" , filename , i ) ;
    FILE *file = fopen( str , "r" ) ;
    if( file == NULL ) { 
      printf( "FILE %s not found !!! \n" , str ) ; 
      exit(1) ;
    }

    int j ;
    double X[ NDATA ] , Y[ NDATA ] ;
    printf( "\n" ) ;
    for( j = 0 ; j < NDATA ; j++ ) {
      double x , y ;
      if( !fscanf( file , "%lf %lf" , &x , &y ) ) {
	exit(1) ;
      }
      printf( "%f %f \n" , x , y ) ;
      Y[j] = y ;
      X[j] = x ;
    }
    double T0_extrap ;
    const double slope = simple_least_squares( X , Y , NDATA , VAL , &T0_extrap ) ;

    if( DER == false ) {
      T0.resampled[ meas ] = T0_extrap ;
    } else {
      printf( "EXTRAP :: %f \n" , 0.3 / slope ) ;
      T0.resampled[ meas ] = 0.3 / slope ;
    }
    fclose( file ) ;
  }

  T0.restype = RAWDATA ;

  compute_err( &T0 ) ;
  printf( "%f %f \n" , T0.avg , T0.err ) ;

  return T0 ;
}
