#ifndef WRITE_DISTRIBUTION_H
#define WRITE_DISTRIBUTION_H

void
write_singledist( struct resampled data ,
		  FILE *file ) ;

void
write_distribution_arr( struct resampled *data ,
			const char *filename ,
			const int NRES ) ;

void
write_distribution_arr2( const double *X , 
			 struct resampled *data ,
			 const char *filename ,
			 const int NRES ) ;
#endif
