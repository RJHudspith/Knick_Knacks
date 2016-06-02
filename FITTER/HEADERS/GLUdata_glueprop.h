#ifndef GLUDATA_GLUEPROP_H
#define GLUDATA_GLUEPROP_H

struct mom_info
fill_mominfo( const int ND ,
	      const int mom[ ND ] ,
	      const int *dimensions ,
	      const momtype type ) ;

struct mom_info *
read_momlist( int *NDATA ,
	      int *ND ,
	      FILE *file ,
	      const int *dimensions ,
	      const momtype mom_type , 
	      const bool dummy ) ;

struct resampled*
read_rawGLU( struct mom_info *mominfo ,
	     struct input_params *INPARAMS ,
	     int nfile ,
	     const momtype mom_type ) ;

// read GLU data file
struct resampled**
read_GLUprop( struct mom_info ***mominfo ,
	      double ***X ,
	      struct input_params *INPARAMS ,
	      int *NSLICES ,
	      const momtype mom_type ) ;

// super distribution
struct resampled**
read_superdist( struct mom_info ***mominfo ,
		double ***X , // is the momentum^2
		struct input_params *INPARAMS ,
		int *NSLICES ,
		const momtype mom_type ) ;

#endif
