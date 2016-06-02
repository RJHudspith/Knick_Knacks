/**
   Analysis for the VPF
 */
#include "fitfunc.h"
#include "fitdata.h"
#include "fit_and_plot.h"
#include "run_boots.h"
#include "dispersions.h"
#include "Utils.h"
#include "write_distribution.h"
#include "read_distribution.h"
#include "effmass.h"
#include "GLU_bswap.h"

// tauvus computation
void
baryon_mass_split_eval( double **xavg ,
			struct resampled **bootavg ,
			struct mom_info **mominfo ,
			struct mom_info *moms ,
			struct input_params *INPARAMS ,
			const int NSLICES ,
			const int LT )
{
  if( NSLICES != 2 ) {
    printf( "Expected two correlators got %d \n" , NSLICES ) ;
    return ;
  }

  // fit and plot our correlators
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bootavg , 
						    (const double**)xavg , 
						    (const struct mom_info **)mominfo ,
						    *INPARAMS , NSLICES , LT ) ;


  // give delta_M for exp and cosh
  subtract( &fitparams[0] , fitparams[2] ) ;
  divide( &fitparams[1] , fitparams[3] ) ;

  printf( "Delta_M :: %e +/- %e \n" , fitparams[0].avg , fitparams[0].err ) ;
  printf( "Amplitde Ratio :: %e +/- %e \n" , fitparams[1].avg , fitparams[1].err ) ;

  return ;
}

