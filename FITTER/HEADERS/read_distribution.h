#ifndef READ_DISTRIBUTION_H
#define READ_DISTRIBUTION_H

struct resampled *
read_distribution( int *NDATA , 
		   FILE *file ) ;

int
read_single_distribution( FILE *file , 
			  struct resampled *RES ,
			  struct mom_info *moms ,
			  const int NSAMPLES ,
			  const resample_type res_check  ,
			  const int *dimensions ,
			  const momtype type ) ;

struct resampled **
read_dispersions( double ***X ,
		  struct input_params *INPARAMS ,
		  struct mom_info ***mominfo ,
		  int *NSLICES ) ;

int
can_flatten( double *array , 
	     struct resampled *res , 
	     const int N ) ;

struct resampled **
read_UKhadron( double ***X ,
	       struct input_params *INPARAMS ,
	       int *NSLICES ) ;

struct resampled **
read_SUSY( double ***X ,
	   struct input_params *INPARAMS ,
	   int *NSLICES ) ;

#endif
