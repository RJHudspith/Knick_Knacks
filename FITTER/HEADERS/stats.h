#ifndef STATS_H
#define STATS_H

void
init_boot_order( const int NBOOTS , 
		 const int NRAW ) ;

void
free_boot_order( void ) ;

// simple numerically-stable average and variance
void
average( double *ave , double *err , 
	 const double *meas , const int N ) ;

void
raw_err( struct resampled *replicas ) ;

// calculator of the error from the jackknifed distribution
void
jackknife_error( struct resampled *replicas ) ;

// Jackknife resample the raw data
void
jackknife( struct resampled *Jackknife ,
	   const struct resampled Raw ) ;

// Compute the bootstrapped data's error ...
void
bootstrap_error( struct resampled *replicas ) ;

// Bootstrap the raw data
void
bootstrap( struct resampled *Bootstrap ,
	   const struct resampled Raw ) ;

void
compute_err( struct resampled *replicas ) ;

struct resampled*
resample_data( const struct resampled *RAW ,
	       const int N ,
	       const resample_type restype ,
	       const int NBOOTS ) ;

void
resample( struct resampled *BOOT ,
	  const int N ) ;

struct resampled **
resample_array( const struct resampled **RAW ,
		const int NSLICES ,
		const int *NDATA ,
		const int resample , 
		const int NBOOTS ) ;

struct resampled
bin_the_data( const struct resampled RAW ,
	      const int binning ) ;

struct resampled
extrapolate( const struct resampled *fitparams ,
	     const double extrap_point ,
	     const int NPARAMS ,
	     const struct chiral quarks ,
	     const struct mom_info mom_info ,
	     const fitfunc fit ,
	     const int LT ,
	     const int SLICE ) ;

#endif
