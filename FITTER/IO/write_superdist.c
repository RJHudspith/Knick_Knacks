/**
   @file write_superdist.c
   @brief super distribution writer
 */

#include "fitfunc.h"

// write out the data of a super-distribution?
void
super_distribution( struct resampled **RAW ,
		    const struct mom_info **mominfo ,
		    const struct input_params INPARAMS ,
		    const int NSLICES ,
		    const int ND ,
		    const bool renormalise )
{
  // loop different samples
  int j , i , nmom = 0 , nsamp = 0 ;
  for( j = 0 ; j < NSLICES ; j++ ) {

    // check for rawness
    /*
    if( RAW[j][0].restype != RAWDATA ) {
      printf( "Data must be raw  \n" ) ;
      return ;
    }
    */

    if( renormalise == true ) {
      // multiply by Z_V
      printf( "\n--> Renormalising <--\n" ) ;

      // renormalise
#pragma omp parallel for private(i)
      for( i = 0 ; i < INPARAMS.NDATA[j] ; i++ ) {
	mult_constant( &RAW[j][i] , INPARAMS.quarks[j].ZV ) ;
      }
    }
    nmom += INPARAMS.NDATA[j] ;

    // compute total number of samples
    nsamp += RAW[j][0].NSAMPLES ;
  }

  printf( "IN DIS \n" ) ;
  
  // check mom lists are all the same size
  if( nmom != NSLICES * INPARAMS.NDATA[0] ) {
    printf( "Error, all momlists should be equal %d vs. %d \n" ,
	    nmom , NSLICES * INPARAMS.NDATA[0] ) ;
    exit(1) ;
  }

  FILE *outfile = fopen( INPARAMS.graph_name , "wb" ) ;

  printf( "--> Writing a SUPER-distribution to %s <--\n" ,  
	  INPARAMS.graph_name ) ;

  // writing out, magic num_mom and momlist is first
  const uint32_t magic[1] = { 717685 } ;
  const uint32_t num_mom[1] = { INPARAMS.NDATA[0] } ;

  fwrite( magic , sizeof(uint32_t) , 1 , outfile ) ;
  fwrite( num_mom , sizeof(uint32_t) , 1 , outfile ) ;

  // write the momentum list
  for( i = 0 ; i < INPARAMS.NDATA[0] ; i++ ) {
    const uint32_t Nd[1] = { ND } ;
    fwrite( Nd , sizeof(uint32_t) , 1 , outfile ) ;
    fwrite( mominfo[0][i].n , sizeof(uint32_t) , ND , outfile ) ;
    // all zero?
  }

  // write out the length of the distribution again
  const uint32_t nsamples[1] = { nsamp } ;
  for( i = 0 ; i < INPARAMS.NDATA[0] ; i++ ) {
    fwrite( nsamples , sizeof(uint32_t) , 1 , outfile ) ;
    for( j = 0 ; j < NSLICES ; j++ ) {
      fwrite( RAW[j][i].resampled , sizeof(double) , RAW[j][i].NSAMPLES , outfile ) ;
    }
  }

  // write out what resampling we used
  const uint32_t resample[1] = { (uint32_t)RAW[0][0].restype } ;
  fwrite( resample , sizeof(uint32_t) , 1 , outfile ) ;

  fclose( outfile ) ;

  return ;
}
