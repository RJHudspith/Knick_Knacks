////////////////// PRONY's method guess, courtesy of Julien Frison ///////////////////
#include "fitfunc.h"

#include <gsl/gsl_poly.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <complex.h>

// creates a polynomial from a row of the cofactor matrix
static void 
blackbox_det( double *coefs , 
	      const double *y , 
	      const int NSTATES )
{
  gsl_matrix *A = gsl_matrix_alloc( NSTATES , NSTATES ) ;
  gsl_permutation *p = gsl_permutation_alloc( NSTATES ) ;
  double det ;
  int signum , i , j , n ;
    
  for( n = 0 ; n < NSTATES+1 ; n++ ) {
    // we fill the A matrix with a reduction of the det related to the X^n term 
    for( i = 0 ; i < NSTATES ; i++ ) {
      for( j = 0 ; j < NSTATES ; j++ ) {
	gsl_matrix_set( A , i , j , y[ ((i<n)?i:(i+1))+j*(NSTATES+1) ] ) ;
      }
    }
    // we compute the det and multiply by the right sign
    gsl_linalg_LU_decomp( A , p , &signum ) ;
    det = gsl_linalg_LU_det( A , signum ) ;
    coefs[n] = ( n&1 ) ? ( det ) :( -det ) ;
  }
  gsl_matrix_free(A);
  gsl_permutation_free(p);
  return ;
}

// qsort comparison
static int 
comp( const void *elem1 , 
      const void *elem2 ) 
{
  const double f = *( (double*)elem1 ) ;
  const double s = *( (double*)elem2 ) ;
  if (f > s) return  1 ;
  if (f < s) return -1 ;
  return 0 ;
}

#if 0
static void
write_matrix( const double *A[] , const int size )
{
  int i , j ;
  printf( "\n" ) ;
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      printf( " %f " , A[i][j] ) ;
    }
    printf("\n" ) ;
  }
  printf( "\n" ) ;
  return ;
}
#endif

// compute the blackbox effective mass for a correlator of y-data
void
blackbox( const double *data ,
	  const int NDATA ,
	  const int NSTATES ,
	  double masses[ NSTATES ][ NDATA ] )
{
  const int smallt = NDATA ; // some cut off t

  // Fleming's notation 
  const int nt = NDATA ;

  size_t t , i , j ;
  for( t = 0 ; t < smallt ; t++ ) {
    // set up the matrix
    double y[ NSTATES*(NSTATES+1) ] ;
    double complex x[ NSTATES ] ;
    double e[ NSTATES ] ;
    for( i = 0 ; i < NSTATES ; i++ ) {
      for( j = 0 ; j < NSTATES+1 ; j++ ) {
	y[ j + i*(NSTATES+1) ] = data[ ( t + j + i )%nt ] ;
      }
    }
    // usual effective mass solution
    if( NSTATES == 1 ) {
      x[0] = y[1] / y[0];
    } else if( NSTATES == 2 ) {
      // we form the polynom ax2+bx+c = det(y0,y1,y2;y3,y4,y5;1,x,x2)
      const double a = y[0] * y[4] - y[1] * y[3] ;
      const double b = y[3] * y[2] - y[0] * y[5] ;
      const double c = y[1] * y[5] - y[2] * y[4] ;
      // we solve it for x1,x2
      const double complex delta = fabs(b) * csqrt( 1.0 - 4.0 * a * c / ( b * b ) ) ;
      x[0] = -( delta - b ) / ( 2.0 * a ) ;
      x[1] =  ( delta + b ) / ( 2.0 * a ) ;
      // otherwise use linear pred
    } else {

      gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc( NSTATES + 1 ) ;
      double *z = (double*)malloc( 2 * NSTATES * sizeof(double*) ) ;
      double coefs[ NSTATES + 1 ] ;

      // compute a row of the cofactor matrix
      blackbox_det( coefs , y , NSTATES ) ;
      gsl_poly_complex_solve( coefs , NSTATES+1 , w , z ) ;

      // 
      for( i = 0 ; i < NSTATES ; i++ ) {
	x[ i ] = z[ 2*i ] + I * z[ 2*i+1 ] ;
      }
      gsl_poly_complex_workspace_free( w ) ;
      free( z ) ;
    }

    // loggy log log .... x[i] = exp( -e[i] )
    for( i = 0 ; i < NSTATES ; i++ ) {
      e[ i ] = creal( -clog( x[ i ] ) ) ;
    }

    // sort the e's ?
    qsort( e , NSTATES , sizeof(double) , comp ) ;
    for( i = 0 ; i < NSTATES ; i++ ) {
      masses[i][t] = e[i] ;
    }
  }

  return ;
}

// driving function, creates a distribution of effective mass
struct resampled **
prony_effmass( const struct resampled **bootavg ,
	       const int *NDATA ,
	       const int NSLICES ,
	       const int NSTATES ,
	       const int STATE )
{
  struct resampled **effmass = malloc( NSLICES * sizeof( struct resampled* ) ) ;
  size_t i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    size_t k , j ;
    effmass[ i ] = malloc( NDATA[i] * sizeof( struct resampled ) ) ;
    for( j = 0 ; j < NDATA[i] ; j++ ) {
      effmass[i][j].resampled = malloc( bootavg[i][j].NSAMPLES * sizeof( double ) ) ;
      effmass[i][j].restype   = bootavg[i][j].restype ;
      effmass[i][j].NSAMPLES  = bootavg[i][j].NSAMPLES ;
    }
    // loop samples
    #pragma omp parallel for private(k)
    for( k = 0 ; k < bootavg[i][0].NSAMPLES ; k++ ) {
      size_t l ;
      // copy the data for this number of samples
      double y[ NDATA[i] ] , masses[ NSTATES ][ NDATA[i] ] ;
      for( l = 0 ; l < NDATA[i] ; l++ ) {
	y[l] = bootavg[i][l].resampled[k] ;
      }
      // compute the black-box solution
      blackbox( y , NDATA[i] , NSTATES , masses ) ;
      // copy first state back into effective mass
      for( l = 0 ; l < NDATA[i] ; l++ ) {
	effmass[i][l].resampled[k] = masses[STATE][l] ;
      }
      //
    }
    // set up the error distribution
    for( j = 0 ; j < NDATA[i] ; j++ ) {
      compute_err( &effmass[i][j] ) ;
    }
  }
  //
  return effmass ;
}

