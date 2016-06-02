#ifndef READ_CORR_H
#define READ_CORR_H

void
copy_inparams( struct input_params *INPARAMS ,
	       const int idx1 ,
	       const int idx2 ) ;

void
expand_inparams( struct input_params *INPARAMS ,
		 const int NCHANNELS ) ;

struct resampled *
read_rawcorr( struct input_params *INPARAMS ,
	      struct mom_info **mom ,
	      const char *filename ,
	      const int nfile ,
	      const foldselection fold ,
	      const int src , 
	      const int snk ) ;

struct resampled *
corrs( double ***X , 
       struct input_params *INPARAMS ,
       struct mom_info **mom ,
       const char *filename ,
       const foldselection fold ,
       const int src , 
       const int snk ,
       const int idx ) ;

#endif
