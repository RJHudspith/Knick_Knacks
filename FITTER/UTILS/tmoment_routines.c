#include "fitfunc.h"

#include "tmoment_routines.h"

#include "svd.h" // my svd wrapper
#include <gsl/gsl_linalg.h> // LU decomp
#include "Utils.h"

//#define WRITE_MATRIX

#define CYCLE_ERR (-3)

//#define LU_INVERT

// factorials
static long double factorial[ 30 ] ;

// initialise?
void
init_factorial( void )
{
  factorial[0] = factorial[1] = 1 ;
  size_t i ;
  for( i = 2 ; i < 30 ; i++ ) {
    factorial[i] = factorial[i-1] * (long double)i ;
  }
  return ;
}

// Computes
//
//   | m |
//   | k |
//
static double
binomial( const int n ,
	  const int k )
{
  if( k < 0 || k > n ) { 
    return 0 ; 
  } else if( k == 0 || k == n ) {
    return 1 ;
  } else if( k == 1 ) {
    return n ;
  } else {
    const int kp = k < n - k ? k : n - k ;
    register long double c = 1. ;
    int i ;
    for( i = 0 ; i < kp ; i++ ) {
      c *= ( n - i ) / (long double)( i + 1 ) ;
    }
    return c ;
  }
  // should probably put a check here
  return -1 ;
}

// cycle of derivatives
static double
cycle_comp( const double q ,
	    const double mhalf ,
	    const int n ,
	    const int k ,
	    const bool tilde ) 
{
  // q->0 limit only these survive
  if( fabs( q ) < 1E-12 ) {
    switch( n&3 ) {
    case 0 : return -1.0 ;
    case 1 : return  0.0 ; 
    case 2 : return  1.0 ; 
    case 3 : return  0.0 ; }
  }
  // generic
  if( tilde ) {
    switch( n&3 ) {
    case 0 : return -cos( 2 * ( mhalf - k )*q ) ;
    case 1 : return  sin( 2 * ( mhalf - k )*q ) ; 
    case 2 : return  cos( 2 * ( mhalf - k )*q ) ; 
    case 3 : return -sin( 2 * ( mhalf - k )*q ) ; }
  } else {
    switch( n&3 ) {
    case 0 : return -cos( ( mhalf - k )*q ) ;
    case 1 : return  sin( ( mhalf - k )*q ) ; 
    case 2 : return  cos( ( mhalf - k )*q ) ; 
    case 3 : return -sin( ( mhalf - k )*q ) ; }
  }
  return CYCLE_ERR ;
}

// uses my general formula for the moments :: computes \Delta^{(n)} \hat{q}^{m} with respect
// to the momenta q, for our use we should specify q = 0.0.
static double
genfac( const int n , 
	const int m , 
	const double q , 
	const int LT ,
	const moment_type moments ,
	const int MomDef )
{
  if( MomDef == PSQ_MOM ) {
    // Matrix is purely diagonal for these
    if( moments == HPQCD_MOMENTS ) {
      return n != m ? 0.0 : 1.0 ;
    } else {
      // something going awry with this!
      if( n > m ) return 0.0 ;
      double twiddle ;
      if( moments == PORTELLI_MOMENTS ) {
	twiddle = M_PI / (double)LT ;
      } else {
	twiddle = 2.0 * M_PI / (double)LT ;
      }
      //twiddle =1 ;
      const double mhalf = (double)n/2. ;
      register double sum = 0.0 ;
      int k ;
      for( k = 0 ; k <= n ; k++ ) {
	sum += k&1 ?\
	  -binomial( n , k ) * pow( q + 2 * twiddle * ( mhalf - k ) , m ) :	\
	   binomial( n , k ) * pow( q + 2 * twiddle * ( mhalf - k ) , m ) ;
      }
      // it took like 2 days to find this factor of 2
      return sum / pow( 2 , n ) ;
    }
  }
  // twosin mom definition is cyclic and that is a bit of an issue
  if( MomDef == TWOSIN_MOM ) {
    // HPQCD MOMENTS
    if( moments == HPQCD_MOMENTS ) {
      // this is really just the LT -> Infinity limit of the one below
      if( m > n ) return 0.0 ;
      // otherwise we go on
      const double mhalf = (double)m/2. ;
      double sum = 0.0 ;
      int k ;
      for( k = 0 ; k <= ( mhalf - 1 ) ; k++ ) {
	sum += binomial( m , k ) *					\
	  pow( -1 , mhalf - k - 1 ) *					\
	  pow( mhalf - k , n ) *					\
	  cycle_comp( q , mhalf , n , k , false )  / factorial[ n ] ;
      }
      return 2.0 * sum ; 
      // DISCRETE MOMENTS
    } else {
      double twiddle ;
      if( moments == PORTELLI_MOMENTS ) {
	twiddle = M_PI / (double)LT ;
      } else {
	twiddle = 2.0 * M_PI / (double)LT ;
      }
      const double mhalf = (double)m/2. ;
      double sum = 0.0 ;
      int k ;
      for( k = 0 ; k <= ( mhalf - 1 ) ; k++ ) {
	sum += binomial( m , k ) *					\
	  pow( -1 , mhalf - k - 1 ) *					\
	  pow( sin( (mhalf-k) * twiddle ) , n ) *			\
	  cycle_comp( q , mhalf , n , k , false ) ;
      }
      return 2.0 * sum ; 
    }
  }
  // twosin mom definition is cyclic and that is a bit of an issue
  if( MomDef == SIN_MOM ) {
    // HPQCD MOMENTS
    if( moments == HPQCD_MOMENTS ) {
      // this is really just the LT -> Infinity limit of the one below
      if( m > n ) return 0.0 ;
      // otherwise we go on
      const double mhalf = (double)m/2. ;
      double sum = 0.0 ;
      int k ;
      for( k = 0 ; k <= ( mhalf - 1 ) ; k++ ) {
	sum += binomial( m , k ) *					\
	  pow( -1 , mhalf - k - 1 ) *					\
	  pow( mhalf - k , n ) *
	  cycle_comp( q , mhalf , n , k , true ) / factorial[ n ] ;
      }
      return pow( 2.0 , n - m + 1 ) * sum ; 
      // DISCRETE MOMENTS
    } else {
      double twiddle ;
      if( moments == PORTELLI_MOMENTS ) {
	twiddle = M_PI / (double)LT ;
      } else {
	twiddle = 2.0 * M_PI / (double)LT ;
      }
      const double mhalf = (double)m/2. ;
      double sum = 0.0 ;
      int k ;
      for( k = 0 ; k <= ( mhalf - 1 ) ; k++ ) {
	sum += binomial( m , k ) *					\
	  pow( -1 , mhalf - k - 1 ) *					\
	  pow( sin( (mhalf-k) * 2 * twiddle ) , n ) * \
	  cycle_comp( q , mhalf , n , k , true ) ;
      }
      return pow( 2 , 1 - m ) * sum ; 
    }
  }
  printf( "WHAT? Should not get here ... genfac \n" ) ;
  exit(1) ;
}

// set matrix ( NxN square matrix )
void
init_moments_matrix( double *mat ,
		     const int N ,
		     const int LT ,
		     const moment_type moments ,
		     const int MomDef )
{
  int i , j ;
  #ifdef WRITE_MATRIX
  printf( "\n* Derivative Matrix *\n" ) ;
  printf( "\n" ) ;
  #endif
  for( i = 0 ; i < N ; i++ ) {
    for( j = 0 ; j < N ; j++ ) {
      mat[ j + N*i ] = genfac( 2*(i+1) , 2*(j+1) , 0.0 , LT , moments , MomDef ) ;
      #ifdef WRITE_MATRIX
      printf( "%e " , mat[ j + N*i ] ) ;
      #endif
    }
    #ifdef WRITE_MATRIX
    printf( "\n" ) ;
    #endif
  }
  #ifdef WRITE_MATRIX
  printf( "\n" ) ;
  #endif
  return ;
}

#ifndef LU_INVERT
// and the actual solving goes here
static void
LU_solver( gsl_matrix_view LUdecomp ,
	   gsl_permutation *p ,
	   double *solutions , 
	   const int N )
{
  // there we go
  gsl_matrix_view mp = LUdecomp ;

  // solution vector
  gsl_vector *x = gsl_vector_alloc( N ) ;

  gsl_vector_view b = gsl_vector_view_array( solutions , N ) ;
  gsl_linalg_LU_solve( &mp.matrix , p , &b.vector , x ) ;

  int i ;
  for( i = 0 ; i < N ; i++ ) {
    solutions[ i ] = gsl_vector_get( x , i ) ;
  }

  gsl_vector_free( x ) ;
  return ;
}
#else
// and the actual solving goes here
static void
LU_invert( gsl_matrix_view LUdecomp ,
	   gsl_permutation *p ,
	   double *minv ,
	   const int N )
{
  // there we go
  gsl_matrix_view mp = LUdecomp ;

  gsl_matrix *inverse = gsl_matrix_alloc( N , N ) ;
  gsl_linalg_LU_invert( &mp.matrix , p , inverse ) ;

  printf( " *Inverted* \n" ) ;

  size_t i , j ;
  for( i = 0 ; i < N ; i++ ) {
    for( j = 0 ; j < N ; j++ ) {
      minv[ j + i*N ] = gsl_matrix_get( inverse , i , j ) ;
      printf( "%e " , minv[ j + i*N ] ) ;
    }
    printf( "\n" ) ;
  }
  printf( "\n" ) ;

  gsl_matrix_free( inverse ) ;
  return ;
}
#endif

// LU decomposition solve for the PIs
//
// schematically ( mats ).( PI ) = ( MOMENTS ), solves for PIS
//
// LU decomp for the matrix mat
// solve by back substitution against MOMENTS to get PIS
// mat is (NxN) matrix and MOMENTS is a N-length vector of bootstraps
void
gsl_inversion( struct resampled *PIS ,
	       double *mat ,
	       const struct resampled *MOMENTS ,
	       const int N ,
	       const moment_type moments ,
	       const int MomDef )
{
  // HPQCD's is the most trivial for this choice
  if( ( moments == HPQCD_MOMENTS ) && ( MomDef == PSQ_MOM ) ) {
    size_t i ; 
    for( i = 0 ; i < N ; i++ ) {
      equate( &PIS[i] , MOMENTS[i] ) ;
    }
    return ;
    // everything else I use an inversion
  } else {

    gsl_matrix_view m = gsl_matrix_view_array( mat , N , N ) ;

    gsl_permutation *p = gsl_permutation_alloc( N ) ;
    
    int s ;
    const int check = gsl_linalg_LU_decomp( &m.matrix , p , &s ) ;

    if( check != 0 ) {
      printf( "[TMOMENTS] LU decomp failed \n" ) ;
      exit(1) ;
    }

    #ifdef LU_INVERT
    double *minv = (double*)malloc( N*N*sizeof( double ) ) ;
    LU_invert( m , p , minv , N ) ;
    #endif

    // this loop allows for computation over all the boots
    int stress ;
    #pragma omp parallel for private(stress)
    for( stress = 0 ; stress < MOMENTS[0].NSAMPLES ; stress++ ) {

      // convert boot-moments to a double vector
      size_t j ;
      double moms[ N ] ;
      for( j = 0 ; j < N ; j++ ) {
	moms[j] = MOMENTS[j].resampled[stress] ;
      }

      #ifdef LU_INVERT
      // apply the inverse to the moms
      for( j = 0 ; j < N ; j++ ) {
	size_t i ;
	register double sum = 0.0 ;
	for( i = 0 ; i < N ; i++ ) {
	  sum += minv[ i + j*N ] * moms[i] ;
	}
	PIS[j].resampled[stress] = sum ;
      }
      #else
      // write solution out into PI
      LU_solver( m , p , moms , N ) ;
      for( j = 0 ; j < N ; j++ ) {
	PIS[j].resampled[stress] = moms[j] ;
      }
      #endif
    }

    // convert boot-moments to a double vector
    size_t j ;
    double moms[N] ;
    for( j = 0 ; j < N ; j++ ) {
      moms[j] = MOMENTS[j].avg ;
    }

    #ifdef LU_INVERT
    for( j = 0 ; j < N ; j++ ) {
      size_t i ;
      register double sum = 0.0 ;
      for( i = 0 ; i < N ; i++ ) {
	sum += minv[ i + j*N ] * moms[i] ;
      }
      printf( "sum %f \n" , sum ) ;
      PIS[j].avg = sum ;
    }
    free( minv ) ;
    #else
    // write solution out into PI
    LU_solver( m , p , moms , N ) ;
    for( j = 0 ; j < N ; j++ ) {
      PIS[j].avg = moms[j] ;
    }
    #endif

    gsl_permutation_free( p ) ;

    // compute err
    for( j = 0 ; j < N ; j++ ) {
      compute_err( &PIS[ j ] ) ;
    }

    return ;
  } 

  // should not get here
  return ;
}

// compute the n^th moment
static void
compute_nth_moment( struct resampled *MOMENTS , 
		    const struct resampled *tcorr ,
		    const int n ,
		    const int LT , 
		    const moment_type moments )
{
  // set the n-1th moment to zero
  equate_constant( &MOMENTS[ n - 1 ] , 0.0 , tcorr[0].NSAMPLES ,
		   tcorr[0].restype ) ;

  // temporary distribution
  struct resampled temp ;
  temp.resampled = malloc( tcorr[0].NSAMPLES * sizeof( double ) ) ;

  // logic for the moment types
  if( moments == HPQCD_MOMENTS ) {
    init_factorial() ;
    // and compute the moments ?
    int t ;
    for( t = -LT/2+1 ; t < LT/2 ; t++ ) {
      // compute the position in the array
      const int posit = ( t + LT )%LT ; 
      //
      equate( &temp , tcorr[ posit ] ) ;
      mult_constant( &temp , (double)pow( t , 2*n ) ) ;
      // if it is in the fit range we use the fit
      add( &MOMENTS[ n - 1 ] , temp ) ;
    }
    // multiply by (-1)^n
    if( n & 1 ) {
      mult_constant( &MOMENTS[ n - 1 ] , -1.0 ) ;
    }
    mult_constant( &MOMENTS[ n - 1 ] , 1.0/factorial[ 2*n ] ) ;
    // all other moment definitions need this
  } else {
    double MULFACT ;
    if( moments == PORTELLI_MOMENTS ) {
      MULFACT = M_PI / (double)LT ;
    } else {
      MULFACT = 2.0 * M_PI / (double)LT ;
    }
    int t ;
    for( t = 0 ; t < LT ; t++ ) {
      // compute the position in the array
      const int posit = t ;
      // set the temp array
      equate( &temp , tcorr[ posit ] ) ;
      mult_constant( &temp , (double)pow( sin( t*MULFACT )  , 2*n )  ) ;
      // add to moments array
      add( &MOMENTS[ n - 1 ] , temp ) ;
    }
    // multiply by (-1)^n
    if( n & 1 ) {
      mult_constant( &MOMENTS[ n - 1 ] , -1.0 ) ;
    }
  }
  free( temp.resampled ) ;
  return ;
}

void
compute_moments( struct resampled *MOMENTS , 
		 const struct resampled *tcorr ,
		 const int NMAX ,
		 const int LT , 
		 const moment_type moments )
{
  int n ;
  for( n = 0 ; n < NMAX ; n++ ) {
    compute_nth_moment( MOMENTS , tcorr , n+1 , LT , moments ) ;
    compute_err( &MOMENTS[ n ] ) ;
    //printf( "[MOMENTS] %e %e \n" , MOMENTS[ n ].avg , MOMENTS[ n ].err ) ;
  }
  return ;
}

// computest the moments and solves for the PIs
struct resampled*
compute_PIS( const struct resampled *tcorr ,
	     const int NMAX ,
	     const int LT , 
	     const moment_type moments ,
	     const int momtype )
{
  // allocate moments
  struct resampled *MOMENTS = malloc( NMAX * sizeof( struct resampled ) ) ;
  struct resampled *PIS = malloc( NMAX * sizeof( struct resampled ) ) ;
  size_t n ;
  for( n = 0 ; n < NMAX ; n++ ) {
    MOMENTS[ n ].resampled = malloc( tcorr[0].NSAMPLES * sizeof( double ) ) ;
    PIS[ n ].resampled = malloc( tcorr[0].NSAMPLES * sizeof( double ) ) ;
  }
  
  // compute moments
  compute_moments( MOMENTS , tcorr , NMAX , LT , moments ) ;

  // do the inversion and solve for the PIS
  double *mat = (double*)malloc( NMAX*NMAX*sizeof( double ) ) ;
  
  init_moments_matrix( mat , NMAX , LT , moments , momtype ) ;

  // set PIS to zero
  for( n = 0 ; n < NMAX ; n++ ) {
    equate_constant( &PIS[ n ] , 0.0 , tcorr[0].NSAMPLES ,
		     tcorr[0].restype ) ;

  }

  // and compute
  gsl_inversion( PIS , mat , MOMENTS , NMAX , moments , momtype ) ;
  
  // free the memories of moments
  free( mat ) ;
  free_resampled( MOMENTS , NMAX ) ;

  return PIS ;
}

// compute the pades
struct resampled**
compute_PADE( struct resampled **PIS ,
	      const int NSLICES ,
	      const int n ,
	      const int m )
{
  struct resampled **PADES = malloc( NSLICES * sizeof( struct resampled* ) ) ;
  size_t i ;
  // compute pade coeffs using the SVD method
  for( i = 0 ; i < NSLICES ; i++ ) {
    PADES[i] = malloc( (n+m+1) * sizeof( struct resampled ) ) ; 
    size_t j ;
    for( j = 0 ; j < (n+m+1) ; j++ ) {
      PADES[i][j].resampled = malloc( PIS[i][0].NSAMPLES * sizeof( double ) ) ;
      // zero
      equate_constant( &PADES[i][j] , 0.0 , 
		       PIS[i][0].NSAMPLES , PIS[i][0].restype ) ;
    }
    // loop boots
    #pragma omp parallel for private(j)
    for( j = 0 ; j < PIS[i][0].NSAMPLES ; j++ ) {
      double p[ n+m+1 ] , coeffs[ n+m+1 ] ;
      // set coeffs
      size_t k ;
      for( k = 0 ; k < ( n + m + 1 ) ; k++ ) {
	coeffs[ k ] = PIS[ i ][ k ].resampled[ j ] ;
      }
      // svd computation of pade coefficients
      pades( p , coeffs , n , m ) ;
      // set pades
      for( k = 0 ; k < ( n + m + 1 ) ; k++ ) {
	PADES[ i ][ k ].resampled[ j ] = p[ k ] ;
      }
    }
    // and the average...
    double p[ n+m+1 ] , coeffs[ n+m+1 ] ;
    // set coeffs
    size_t k ;
    for( k = 0 ; k < ( n + m + 1 ) ; k++ ) {
      coeffs[ k ] = PIS[ i ][ k ].avg ;
    }
    // svd computation of pade coefficients
    pades( p , coeffs , n , m ) ;
    // set the average
    for( k = 0 ; k < ( n + m + 1 ) ; k++ ) {
      PADES[ i ][ k ].avg = p[ k ] ;
      compute_err( &PADES[ i ][ k ] ) ;
    }
    //
  }
  return PADES ;
}
