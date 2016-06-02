/**
   @file mass_splittings.h
   @brief compute various subtractions and ratios
 */
#ifndef MASS_SPLITTINGS_H
#define MASS_SPLITTINGS_H

void
mass_splittings_eval( double **xavg ,
		      struct resampled **bootavg ,
		      struct mom_info **mominfo ,
		      struct input_params *INPARAMS ,
		      const int NSLICES ,
		      const int LT ) ;

#endif
