/**
   @file mass_splittings.c
   @brief computes ratios of masses
 */
#include "fitfunc.h"
#include "fit_and_plot.h"
#include "resampled_ops.h"

void
mass_extrap( double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT )
{
  size_t i ;
  for( i = 0 ; i < INPARAMS->NDATA[0] ; i++ ) {
    //mult( &bootavg[0][i] , bootavg[0][i] ) ;
    xavg[0][i] = INPARAMS -> quarks[i].ml ;
    printf( "[MASSextrap] %e :: %e +/- %e \n" , 
	    xavg[0][i] , bootavg[0][i].avg , bootavg[0][i].err ) ;
  }

  // fit the data
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bootavg ,
						    (const double**)xavg , 
						    (const struct mom_info **)mominfo ,
						    *INPARAMS , 1 , LT ) ;

  if( INPARAMS -> fittype == POLY1 ) {
    // evaluate the exact chiral point
    struct resampled zero ;
    zero.resampled = malloc( bootavg[0][0].NSAMPLES * sizeof( double ) ) ;

    equate( &zero , fitparams[0] ) ;
    mult_constant( &zero , -1 ) ;
    divide( &zero , fitparams[1] ) ;

    printf( "Zero point :: %1.12e +/- %1.12e \n" , zero.avg , zero.err ) ;

    equate( &zero , fitparams[1] ) ;
    mult_constant( &zero , INPARAMS -> fit_lo ) ;
    add( &zero , fitparams[0] ) ;

    printf( "Zero point :: %e %e %e\n" , zero.avg , 
	    zero.avg - zero.err_hi , zero.avg - zero.err_lo ) ;
    //raise( &zero , 0.5 ) ;

    mult_constant( &zero , INPARAMS -> dimensions[0][0] ) ;
    printf( "MRHOL :: %1.12e +/- %1.12e \n" , zero.avg , zero.err ) ;

    free( zero.resampled ) ;
  }

  free( (void*)fitparams ) ;

  return ;
}
