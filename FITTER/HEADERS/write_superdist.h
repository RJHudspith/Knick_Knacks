/**
   @file write_superdist.h
   @brief IO for creating a super-distribution
 */
#ifndef WRITE_SUPERDIST_H
#define WRITE_SUPERDIST_H

void
super_distribution( struct resampled **RAW ,
		    const struct mom_info **mominfo ,
		    const struct input_params INPARAMS ,
		    const int NSLICES ,
		    const int ND ,
		    const bool renormalise ) ;

#endif
