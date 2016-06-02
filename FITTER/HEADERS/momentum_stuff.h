/**
   @file momentum_stuff.h
   @brief commmon momentum-based routines
 */
#ifndef MOMENTUM_STUFF_H
#define MOMENTUM_STUFF_H

void
convert_to_physmom( struct mom_info **mavg , 
		    double **x ,
		    const struct input_params *INPARAMS , 
		    const size_t NSLICES ) ;

#endif
