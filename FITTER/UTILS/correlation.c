/**
   @file correlation.c
   @brief computes the correlation matrix of a dataset
 */
#include "fitfunc.h"

#include "svd.h"

//#define DIAGONAL

void
write_corrmatrix_to_file( FILE *outfile , 
			  const double **correlation ,
			  const int NCUT )
{
  // print out correlations
  size_t i , j ;
  for( i = 0 ; i < NCUT ; i++ ) {
    for( j = 0 ; j < NCUT ; j++ ) {
      fprintf( outfile , "%1.3E " , correlation[i][j] ) ;
    }
    fprintf( outfile , "\n" ) ;
  }
  return ;
}

void
write_corrmatrix_mathematica( double **correlation ,
			      const int NCUT )
{
  // print out correlations
  size_t i , j ;
  printf( "\n{{" ) ;
  for( i = 0 ; i < NCUT-1 ; i++ ) {
    for( j = 0 ; j < NCUT-1 ; j++ ) {
      printf( "%1.10f," , correlation[i][j] ) ;
    }     
    printf( "%1.10f},{\n" , correlation[i][j] ) ;
  }
  for( j = 0 ; j < NCUT-1 ; j++ ) {
    printf( "%1.10f," , correlation[i][j] ) ;
  }
  printf( "%1.10f}}\n" , correlation[i][j] ) ;
  return ;
}

// computes the upper section
static void
compute_upper_correlation( double **correlation , 
			   const struct resampled *data ,
			   const int NDATA ,
			   const bool diagonal )
{
  const int NSAMPLES = data[0].NSAMPLES ;

  // compute the correct normalisation
  double NORM = 1.0 ;
  switch( data[0].restype ) {
  case RAWDATA :
    NORM = 1.0 / (double)( NSAMPLES * ( NSAMPLES - 1.0 ) ) ;
    break ;
  case JACKDATA :
    NORM = ( NSAMPLES - 1.0 ) / (double)NSAMPLES ;
    break ;
  case BOOTDATA :
    NORM = 1.0 / (double)NSAMPLES ;
    break ;
  }

  // set to zero
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < NDATA ; i++ ) {
    size_t j ;
    for( j = 0 ; j < NDATA ; j++ ) {
      correlation[i][j] = 0.0 ;
    }
  }

  // diagonal one is pretty straightforward
  if( diagonal == true ) {
#pragma omp parallel for private(i)
    for( i = 0 ; i < NDATA ; i++ ) {     
      const register double ave = data[i].avg ;
      register double sum = 0.0 ;
      size_t k ;
      for( k = 0 ; k < NSAMPLES ; k++ ) {
	sum += 
	  ( data[i].resampled[k] - ave ) *  
	  ( data[i].resampled[k] - ave ) ;
      }
      correlation[i][i] = sum * NORM ;
    }
    // off diagonal is only a bit more complicated
  } else {
#pragma omp parallel for private(i)
    for( i = 0 ; i < NDATA ; i++ ) {
      const double avei = data[i].avg ;
      size_t j ;
      for( j = i ; j < NDATA ; j++ ) {
	const double avej = data[j].avg ;
	register double sum = 0.0 ;
	size_t k ;
	for( k = 0 ; k < NSAMPLES ; k++ ) {
	  sum += 
	    ( data[i].resampled[k] - avei ) *  
	    ( data[j].resampled[k] - avej ) ;
	}
	correlation[i][j] = sum * NORM ;
      }
    }
    // done
  }

  return ;
}

// divide by sigma^2
static void
modified_covariance( double **correlation ,
		     const int NDATA )
{
  // compute sigma
  double sigma[ NDATA ] ;
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < NDATA ; i++ ) {
    sigma[ i ] = sqrt( correlation[i][i] ) ;
  }
  
  // rescale correlation matrix by variance
#pragma omp parallel for private(i)
  for( i = 0 ; i < NDATA ; i++ ) {
    size_t j ;
    for( j = i ; j < NDATA ; j++ ) {
      // edge case that I quite often hit
      if( sigma[i] == 0.0 || sigma[j] == 0.0 ) {
	if( j == i ) {
	  correlation[i][i] = 1.0 ;
	} else {
	  correlation[i][i] = 0.0 ;
	}
      // otherwise rescale
      } else {
	correlation[i][j] /= ( sigma[i] * sigma[j] ) ;
      }
      //
    }
  }
  return ;
}

// fill in the lower half by symmetry
static void
fill_lower_triangular( double **correlation ,
		       const int NDATA )
{
  // fill in the rest by symmetry
  int i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < NDATA ; i++ ) {
    int j ;
    for( j = 0 ; j < i ; j++ ) {
      correlation[i][j] = correlation[j][i] ;
    }
  }
  return ;
}

// compute the covariance matrix, normalised by sigma_i sigma_j
void
correlations( double **correlation , 
	      const struct resampled *data ,
	      const int NDATA )
{
  // compute this correlation matrix
  compute_upper_correlation( correlation , data , NDATA , false ) ;

  // divide by sigma^2
  modified_covariance( correlation , NDATA ) ;

  // correlation matrix is symmetric
  fill_lower_triangular( correlation , NDATA ) ;

  return ;
}

// compute the covariance matrix, normalised by sigma_i sigma_j
double **
correlations_inv( const struct resampled *data ,
		  const int NDATA )
{
  double **Cinv = malloc( NDATA * sizeof( double* ) ) ;
  double **C = malloc( NDATA * sizeof( double* ) ) ;

  int i ;
  for( i = 0 ; i < NDATA ; i++ ) {
    Cinv[ i ] = malloc( NDATA * sizeof( double ) ) ;
    C[ i ] = malloc( NDATA * sizeof( double ) ) ;
  }

  // compute this correlation matrix (false means not diagonal)
  compute_upper_correlation( C , data , NDATA , false ) ;

  // divide by sigma^2, this ruins the chisq!!!
  //modified_covariance( C , NDATA ) ;

  // correlation matrix is symmetric
  fill_lower_triangular( C , NDATA ) ;

#ifdef DIAGONAL
  for( i = 0 ; i < NDATA ; i++ ) {
    size_t j ;
    for( j = 0 ; j < NDATA ; j++ ) {
      if( i == j ) {
	Cinv[i][i] = 1.0 / ( C[i][i] ) ;
      } else {
	Cinv[i][j] = 0.0 ;
      }
      //
    }
  }
#else
  // compute the inverse
  svd_inverse( Cinv , (const double**)C , NDATA , NDATA ) ;
#endif

  for( i = 0 ; i < NDATA ; i++ ) {
    free( C[i] ) ;
  }
  free( C ) ;

  // print him out
  //write_corrmatrix( Cinv , NDATA ) ;

  return Cinv ;
}
