#ifndef AUTOCORR_H
#define AUTOCORR_H

/**
   @fn int autocorrelation( const struct resampled RAW , const int NSEP , const char *output )
   @brief computes the integrated autocorrelation time of the raw sample
   @param RAW :: Raw data distribution
   @param NSEP :: measurement separation of each RAW data point
   @param output :: output file name for the autocorrelations
   
   @warning writes out a file and needs fftw at the moment
 */
int
autocorrelation( const struct resampled RAW ,
		 const int NSEP ,
		 const char *output ) ;

#endif
