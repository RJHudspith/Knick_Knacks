#ifndef STATS_H
#define STATS_H

void
average( double *ave , double *err , 
	 const double *meas , const int N ) ;

void
jackknife_error( struct resampled *replicas ) ;

void
jackknife( struct resampled *Jackknife ,
	   const struct resampled Raw ) ;

void
equate( struct resampled *a , 
	const struct resampled b ) ;

void
subtract( struct resampled *a , 
	  const struct resampled b ) ;

void
init_stats( struct resampled *raw , struct resampled *rawsq ,
	    struct resampled *jack , struct resampled *jacksq ,
	    struct resampled *susc , const int MEASUREMENTS ) ;

#endif
