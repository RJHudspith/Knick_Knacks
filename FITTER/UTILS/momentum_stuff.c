/**
   @file momentum_stuff.c
   @brief momentum routines I use quite often
 */
#include "fitfunc.h"

// convert to physical momentum
void
convert_to_physmom( struct mom_info **mavg , 
		    double **x ,
		    const struct input_params *INPARAMS , 
		    const size_t NSLICES ) 
{
  size_t i , j ;
  // OK, momentum first
  for( j = 0 ; j < NSLICES ; j++ ) {
    // ainverse multipliers
    const double aI  = INPARAMS -> quarks[j].ainverse ;
    const double aI2 = pow( INPARAMS -> quarks[j].ainverse , 2 ) ;
    const double aI4 = pow( INPARAMS -> quarks[j].ainverse , 4 ) ;
    const double aI6 = pow( INPARAMS -> quarks[j].ainverse , 6 ) ;

    #pragma omp parallel for private(i)
    for( i = 0 ; i < INPARAMS -> NDATA[j] ; i++ ) {

      // set to physical lattice spacing xavg[j][i] is now |p| in GeV!!
      x[j][i] = aI * sqrt( x[j][i] ) ;
      // and the moms
      mavg[j][i].p2 *= aI2 ;
      mavg[j][i].p4 *= aI4 ;
      mavg[j][i].p6 *= aI6 ;

      // check for some consistency
      if( fabs( mavg[j][i].p2 - x[j][i]*x[j][i] ) > 1E-12 ) {
	printf( "P2 Broken !! %e \n" , mavg[j][i].p2 - x[j][i]*x[j][i] ) ;
	exit(1) ;
      } 
    }
  }
  return ;
}
