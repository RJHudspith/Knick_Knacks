#ifndef UTILS_H
#define UTILS_H

int 
find_idx( const double target , 
	  const double *X , 
	  int pivot , 
	  int diff ) ;

double *
set_params( const gsl_vector *x , 
	    const int NPARAMS ) ;

void
free_resampled( struct resampled *data , 
		const int NDATA ) ;

void
free_resampled_dp( struct resampled **data , 
		   const int idx1 , 
		   const int *idx2 ) ;

void
free_double_dp( double **data , 
		const int idx1 ) ;

void
free_mominfo( struct mom_info **data , 
	      const int idx1 ) ;

void
printboots( const struct resampled R ) ;

struct resampled*
flatten_resampled_array( const struct resampled **R ,
			 const int idx1 ,
			 const int *idx2 ) ;

double*
flatten_double_array( const double **R ,
		      const int idx1 ,
		      const int *idx2 ) ;

struct mom_info*
flatten_mom_array( const struct mom_info **R ,
		   const int idx1 ,
		   const int *idx2 ) ;

double*
resampled_to_double( const struct resampled *R ,
		     const int size , 
		     const int which ) ;

#endif
