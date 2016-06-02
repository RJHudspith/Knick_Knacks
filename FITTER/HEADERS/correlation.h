#ifndef GLUFIT_CORRELATION_H
#define GLUFIT_CORRELATION_H

void
write_corrmatrix_to_file( FILE *outfile , 
			  double **correlation ,
			  const int NCUT ) ;

void
write_corrmatrix_mathematica( double **correlation ,
			      const int NCUT ) ;

void
correlations( double **correlation , 
	      const struct resampled *data ,
	      const int NDATA ) ;

// inverse correlation matrix
double **
correlations_inv( const struct resampled *data ,
		  const int NDATA ) ;

#endif
