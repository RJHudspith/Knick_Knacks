/**
   @file pade_from_data.c
   @brief computes directly a pade representation from some data

   Direct Pade evaluation using the formula

   Z^-1(y*Z) b = a

   With Z being the usual Vandermonde Z-matrix

   | 1 q_1 q_1^2 ... |
   | 1 q_2 q_2^2 ... | 
   | 1 ............. |
   | 1 ............. |
 */
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define SVD_COL_BALANCE

#define FAILURE 0
#define SUCCESS !FAILURE

static gsl_rng *r ;
static long unsigned int SEED = -1 ;

double
rng_double( void ) { return gsl_rng_uniform( r ) ; }

// between 0 and some number "max_idx"
int
rng_int( const int max_idx ) { 
  return gsl_rng_uniform_int( r , max_idx ) ;
}

double
rng_gaussian( const double sigma ) { 
  return gsl_ran_gaussian( r , sigma ) ; 
}

void
init_rng( long unsigned int Seed )
{
  gsl_rng_env_setup( ) ;
  
  const gsl_rng_type *type = gsl_rng_default ;
  r = gsl_rng_alloc (type) ;

  if( Seed == 0 ) {
    FILE *urandom = fopen( "/dev/urandom" , "r" ) ;
    if( urandom == NULL ) exit( 1 ) ;
    if( fread( &Seed , sizeof( Seed ) , 1 , urandom ) != 1 ) exit(1) ;
    fclose( urandom ) ;
  }

  SEED = Seed ;
  printf( "USING GSL seed %lu \n" , SEED ) ;

  gsl_rng_set( r , Seed ) ;

  return ;
}

void
rng_reseed( void ) 
{
  if( SEED != -1 ) {
    gsl_rng_set( r , SEED ) ;
  } else {
    printf( "GSL rng not seeded properly!\n" );
    exit(1) ;
  }
  return ;
}

void
free_rng( void )
{
  gsl_rng_free( r );
  return ;
}

// little print utility
#ifdef VERBOSE
static void
printmatrix( const double **mat ,
	     const int NROWS , 
	     const int NCOLS )
{
  printf( "\n" ) ;
  int i , j ;
  for( i = 0 ; i < NROWS ; i++ ) {
    printf( "|" ) ;
    for( j = 0 ; j < NCOLS ; j++ ) {
      printf( " %f " , mat[i][j] ) ;
    }
    printf( "|\n" ) ;
  }
  printf( "\n" ) ;
  return ;
}
#endif

/**
   Computes A = U S V^T for a general M * N matrix where M > N
 */
static int 
svd_inverse( double *__restrict *__restrict Ainv , 
	     const double *__restrict *__restrict A ,
	     const int NCOLS ,
	     const int NROWS )
{
  // initial case dies
  if( NROWS < NCOLS ) {
    printf( "WHAT M :: %d vs %d \n" , NROWS , NCOLS ) ;
    return FAILURE ;
  }

  // allocations
  gsl_matrix *m  = gsl_matrix_alloc( NROWS , NCOLS ) ;
  gsl_matrix *Q  = gsl_matrix_alloc( NCOLS , NCOLS ) ;
  gsl_vector *S  = gsl_vector_alloc( NCOLS ) ;
  gsl_vector *WORK  = gsl_vector_alloc( NCOLS ) ;
  double Diag[ NCOLS ] , InvDiag[ NCOLS ] ;
  double diff = 0.0 ;
  register double sum ;
  double tmp ;
  int i , j , k , FLAG = SUCCESS ;

  for( i = 0 ; i < NROWS ; i++ ){
   for( j = 0 ; j < NCOLS ; j++ ) {
     gsl_matrix_set( m , i , j , A[i][j] ) ;
   }
  }

  // do the decomposition
#ifdef SVD_COL_BALANCE
  gsl_vector *D = gsl_vector_alloc( NCOLS ) ;
  gsl_linalg_balance_columns( m , D ) ;
#endif

  if( gsl_linalg_SV_decomp( m , Q , S , WORK ) ) {
    printf( "GSL SVD comp failure \n" ) ;
    FLAG = FAILURE ;
    goto FREE ;
  }

  // loop these set 1.0 / Diag[i] to 0.0 if Diag[i] is crazy small
  for( i = 0 ; i < NCOLS ; i++ ) {
    Diag[ i ] = gsl_vector_get( S , i ) ; 
    tmp = 1.0 / Diag[i] ;
    InvDiag[ i ] = fabs( Diag[i] ) < 1.E-32 ? 0.0 : tmp ;
    #ifdef VERBOSE
    printf( "SVD %d %le\n" , i , Diag[i] ) ;
    printf( "SVD %d %le\n" , i , InvDiag[i] ) ;
    #endif
  }

  // test the solution to make sure it isn't too bad
  for( i = 0 ; i < NROWS ; i++ ) {
    for( k = 0 ; k < NCOLS ; k++ ) {
      sum = 0.0 ;
      for( j = 0 ; j < NCOLS ;j++ ) {
	sum += gsl_matrix_get( m , i , j ) * Diag[j] * gsl_matrix_get( Q , k , j ) ;
      }
      #ifdef SVD_COL_BALANCE
      diff += fabs( sum - A[i][k] / gsl_vector_get( D , k ) ) ;
      #else
      diff += fabs( sum - A[i][k] ) ;
      #endif
    }
  }

  // tell us how good the solution is and cry if it is too bad
  printf( "Decomposition accuracy :: %le \n" , diff ) ;
  diff /= (double)( NCOLS * NROWS ) ;
  if( diff > 1E-8 ) {
    printf( "SVD accuracy considered too low %e \n" , diff ) ;
    FLAG = FAILURE ;
    goto FREE ;
  } 

  // compute the product Ainv = V * ( 1.0 / Diag ) * U^T
  for( i = 0 ; i < NCOLS ; i++ ) {
    for( k = 0 ; k < NROWS ; k++ ) {
      sum = 0.0 ;
      for( j = 0 ; j < NCOLS ;j++ ) {
	if( InvDiag[ j ] != 0.0 ) {
	  sum += gsl_matrix_get( Q , i , j ) * InvDiag[j] * gsl_matrix_get( m , k , j ) ;
	}
      }
      #ifdef SVD_COL_BALANCE
      Ainv[i][k] = sum / gsl_vector_get( D , i ) ;
      #else
      Ainv[i][k] = sum ;
      #endif
    }
  }

  #ifdef VERBOSE
  // test our solution, should equal the NxN identity
  double **sol = malloc( NCOLS * sizeof( double* ) ) ;
  for( i = 0 ; i < NCOLS ; i++ ) {
    sol[i] = (double*)malloc( NCOLS * sizeof( double ) ) ;
    for( j = 0 ; j < NCOLS ; j++ ) {
      sum = 0.0 ;
      for( k = 0 ; k < NROWS ; k++ ) {
	sum += Ainv[i][k] * A[k][j] ;
      }
      sol[i][j] = sum ;
    }
  }
  printmatrix( (const double**)A , NROWS , NCOLS ) ;
  printmatrix( (const double**)Ainv , NCOLS , NROWS ) ;
  printmatrix( (const double**)sol , NCOLS , NCOLS ) ;

  for( i = 0 ; i < NCOLS ; i++ ) {
    free( sol[i] ) ;
  }
  free( sol ) ;
  #endif

 FREE :
  gsl_matrix_free( m ) ;
  gsl_matrix_free( Q ) ;
  gsl_vector_free( S ) ;
  gsl_vector_free( WORK ) ;

#ifdef SVD_COL_BALANCE
  gsl_vector_free( D ) ;
#endif

  return FLAG ;
}

// pade coefficient function
static int
pades( double *pade_coeffs ,
       const double **AinvYA ,
       const int n ,
       const int m )
{
  // the matrix order equation we must solve
  const int order = n + m + 1 ;

  // initialise solve matrix and its inverse
  double **sub_solve = malloc( m * sizeof( double* ) ) ;
  double **Ainv = malloc( m * sizeof( double* ) ) ;

  // vectors we use
  double *solves = malloc( m * sizeof( double ) ) ;

  // flag for whether it worked
  int FLAG = FAILURE ;

  // set up the submatrix we need to solve for the
  // left hand side
  int i , j ;
  for( i = 0 ; i < m ; i++ ) {
    sub_solve[ i ] = malloc( m * sizeof( double) ) ;
    Ainv[ i ] = malloc( m * sizeof( double) ) ;
    for( j = 0 ; j < m ; j++ ) {
      sub_solve[ i ][ j ] = AinvYA[ i + n ][ j + 1 ] ;
    }
    solves[i] = -AinvYA[ n + i ][ 0 ] ;
    #ifdef VERBOSE
    printf( "sols [ %d ] = %d :: %f \n" , i , m + i , solves[i] ) ;
    #endif
  }

  // printhimout
#ifdef VERBOSE
  printmatrix( (const double**)sub_solve , m , m ) ;
#endif

  // perform an inverse
  if( svd_inverse( Ainv , (const double**)sub_solve , m , m ) == FAILURE ) {
    FLAG = FAILURE ;
    goto FREE ;
  }

  // multiply ainverse by sols to get left hand side coefficients
  pade_coeffs[ n ] = 1.0 ;
  for( i = 1 ; i < order - 1 ; i++ ) {
    if( i <= m ) {
      for( j = 0 ; j < m ; j++ ) {
	pade_coeffs[ n + i ] += Ainv[i-1][j] * solves[j] ;
      }
    }
  }

  // multiply coefficient matrix by left hand side
  for( i = 0 ; i < n ; i++ ) {
    register double sum = 0.0 ;
    for( j = 0 ; j < m+1 ; j++ ) {
      sum += AinvYA[ i ][ j ] * pade_coeffs[ n + j ] ;
    }
    pade_coeffs[ i ] = sum ;
    #ifdef VERBOSE
    printf( "\n" ) ;
    printf( "right[ %d ] = %f \n" , i , pade_coeffs[ i ] ) ;
    #endif
  }

  // numerator is first in our scheme
  printf( "( (%1.2f) +" , pade_coeffs[ 0 ] ) ;
  for( i = 1 ; i < n-1 ; i++ ) {
    printf( " (%1.2f) q^{%d} +" , pade_coeffs[ i ] , 2*(i) ) ;
  }
  printf( " (%1.2f) q^{%d} ) / ( " , pade_coeffs[ i ] , 2*(i+1) ) ;

  // denominator is last
  printf( "( (%1.2f) +" , pade_coeffs[ n ] ) ;
  for( i = 1 ; i < m ; i++ ) {
    printf( " (%1.2f) q^{%d} +" , pade_coeffs[ i + n ] , 2*(i+1) ) ;
  } 
  printf( " (%1.2f) q^{%d} ) \n" , pade_coeffs[ i + n ] , 2*(i+1) ) ;

  // free workspaces
 FREE :
  for( i = 0 ; i < m ; i++ ) {
    free( sub_solve[i] ) ;
    free( Ainv[i] ) ;
  }
  free( sub_solve ) ;
  free( Ainv ) ;
  free( solves ) ;

  return FLAG ;
}

// compute ( Z^{-1} ( y*Z ) )
double **
compute_Zinv_y_Z( const double *xdata ,
		  const double *ydata ,
		  const size_t NROWS , // this is the number of data
		  const size_t NCOLS ) // this is the order of the pade
{
  // set up a temporary y
  double yt[ NROWS ] ;
  size_t i , j , k ;
  for( i = 0 ; i < NROWS ; i++ ) {
    yt[ i ] = ydata[ i ] ;
  }
  // set Z
  double **Z = malloc( NROWS * sizeof( double* ) ) ;
  for( i = 0 ; i < NROWS ; i++ ) {
    Z[i] = malloc( NCOLS * sizeof( double ) ) ;
  }

  for( i = 0 ; i < NROWS ; i++ ) {
    register double xs = 1 , xx = xdata[i] ;
    for( j = 0 ; j < NCOLS ; j++ ) {
      Z[ i ][ j ] = xs ;
      xs *= ( xx ) ;
    }
  }

  // set up and compute the inverse
  double **Zinv = malloc( NCOLS * sizeof( double* ) ) ;
  for( i = 0 ; i < NCOLS ; i++ ) {
    Zinv[i] = malloc( NROWS * sizeof( double ) ) ;
  }

  if( svd_inverse( Zinv , (const double**)Z , NCOLS , NROWS ) == FAILURE ) {
    return NULL ;
  }

#ifdef VERBOSE
  printf( "Z and Zinv\n" ) ;
  printmatrix( (const double**)Z , NROWS , NCOLS ) ;
  printmatrix( (const double**)Zinv , NCOLS , NROWS ) ;
#endif

  // form Zinv * ( y * Z ) 
  for( i = 0 ; i < NROWS ; i++ ) {
    for( j = 0 ; j < NCOLS ; j++ ) {
      Z[ i ][ j ] *= yt[ i ] ;
    }
  }

#ifdef VERBOSE
  printmatrix( (const double**)Z , NROWS , NCOLS ) ;
#endif

  // matrix product ( Z^{-1}( y * Z ) )
  for( i = 0 ; i < NCOLS ; i++ ) {
    size_t j , k ;
    for( j = 0 ; j < NCOLS ; j++ ) {
      register double sum = 0.0 ;
      for( k = 0 ; k < NROWS ; k++ ) {
	sum += Zinv[i][k] * Z[k][j] ;
      }
      yt[ j ] = sum ;
    }
    for( j = 0 ; j < NCOLS ; j++ ) {
      Zinv[i][j] = yt[j] ;
    }
  }

  return Zinv ;
}

// generate some fake data
static void
fakedata( double *x , 
	  double *y ,
	  const size_t NDATA ,
	  const size_t n ,
	  const size_t m )
{
  // compute Z.a and Z.b
  double b[ m + 1 ] , a[ n ] ;

  // init rng
  init_rng( 1234 ) ;

  size_t i , j ;
  for( i = 0 ; i < n ; i++ ) {
    a[ i ] = rng_double() ;
    printf( "a[ %zu ] :: %f \n" , i , a[i] ) ;
  }
  for( i = 0 ; i < m+1 ; i++ ) {
    b[ i ] = i == 0 ? 1 : rng_double() ;
    printf( "b[ %zu ] :: %f \n" , i , b[i] ) ;
  }

  for( i = 0 ; i < NDATA ; i++ ) {
    double num = 0.0 , denom = 0.0 ;
    for( j = 0 ; j < n ; j++ ) {
      num += a[j] * pow( i+1 , j ) ;
    }
    for( j = 0 ; j < m+1 ; j++ ) {
      denom += b[j] * pow( i+1 , j ) ;
    }
    y[i] = num / denom ;
    x[i] = i + 1 ;
    printf( "DATA :: %f %f \n" , x[i] , y[i] ) ;
  }
  free_rng( ) ;

  return ;
}

//
int main( void )
{
  // general n, m pade
  const int n = 1 ;
  const int m = 1 ;
  const int order = n + m + 1; // number of coefficients

  // create some fake data
  const size_t NDATA = 15 ;

  double x[ NDATA ] , y[ NDATA ] ;

  // create some fake data
  fakedata( x , y , NDATA , n , m ) ;

  // compute the important matrix
  double **Zinv = compute_Zinv_y_Z( x , y , NDATA , order ) ;

  double pade_coeffs[ order ] ; // outputs
  pades( pade_coeffs , (const double**)Zinv , n , m ) ;

  return 0 ;
}
