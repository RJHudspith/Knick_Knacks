#include "fitfunc.h"
#include "fitdata.h" // fit function

// for the set_mu and set_loops functions
#include "alphas.h"
#include "run_boots.h" // running of the coupling

#if 0
void
run_fit_and_plot( struct resampled **RES ,
		  double **flatx ,
		  double **sigma ,
		  struct input_params INPARAMS ,
		  const int NSLICES ,
		  const int LT )
{
  int i ;

  const int loops = 3 ;
  set_loops( loops ) ;

  // set up a minimum and a maximum ...
  const double mumin = 1.0 ;
  const double mumax = 3.0 ;
  const double gran  = 0.1 ;
  const int NMU = (int)( ( mumax - mumin ) / gran ) ;
  int counter = 0 ;

  // allocate X and Y
  double *X = malloc( NMU * sizeof( double ) ) ;
  struct resampled *MZ = malloc( NMU * sizeof( struct resampled ) ) ;

  double mu ;
  for( mu = mumin ; mu < mumax ; mu+=gran ) {

    X[ counter ] = mu ;

    // important parameters!!
    set_mu( mu ) ;
  
    // make sure the fit makes sense
    INPARAMS.fit_lo = ( INPARAMS.fit_lo < 0 ) ? 0 : INPARAMS.fit_lo ;
    INPARAMS.fit_hi = ( INPARAMS.fit_hi > INPARAMS.NDATA - 1 ) ?	\
      INPARAMS.NDATA-1 : INPARAMS.fit_hi ;

    // perform the fit
    int *NPARAMS = malloc( NSLICES * sizeof( int ) ) ;
    fitfunc fit ;
    struct resampled **fitparams = perform_bootfit( &fit ,
						    NPARAMS ,
						    (const struct resampled**)RES , // the bootstrapped data 
						    (const double**)flatx , // the (sorted) x-axis data
						    (const double**)sigma , // the s.d of the y-data
						    INPARAMS.quarks , // chiral data
						    INPARAMS.fittype , // the type of fit(s) we are using ...
						    1 , // number of BOOTs to fit
						    INPARAMS.NDATA , // the length of the X-array
						    INPARAMS.fit_hi ,
						    INPARAMS.fit_lo ,
						    NSLICES ,
						    LT ) ;

    // unpack the simultaneous fit into something more palatable
    struct resampled fit1 ;

    fit1.resampled = (double*)malloc( fitparams[0][0].NSAMPLES * sizeof( double ) ) ;

    equate( &fit1 , fitparams[0][0] ) ;

    mult_constant( &fit1 , M_PI ) ;

    printf( "Running from %f :: %f %f \n" , mu , fit1.avg , fit1.err ) ;

    // run fit_0 to MZ

    MZ[ counter ] = boot_run_MS_quick( mu , fit1 , loops , 3 , true ) ;

    compute_err( &MZ ) ;

    printf( "%d-LOOP :: %f %f \n" , loops ,  MZ[counter].avg , MZ[counter].err ) ;

    free( fit1.resampled ) ;

    counter++ ;
  }


  // free( X ) ;
  //free( MZ.resampled ) ;

  return ;
}
#endif

  // OK, use the fit to create a semi-continuum dataset
  /*
  double x ;
  const double XMIN = 1.0 ;
  const double XMAX = 3.0 ;
  const double GRAN = 0.01 ;
  const int NMEAS = (int)( ( XMAX - XMIN ) / GRAN ) ;
  struct resampled *Y = malloc( NMEAS * sizeof( struct resampled ) ) ;
  double *XXX = malloc( NMEAS * sizeof( double ) ) ;
  int idx = 0 ;
  for( x = XMIN ; x <= XMAX ; x+= GRAN ) {

    // extrapolate the fit parameters
    Y[idx] = extrapolate( fit1 , x ,
			  NPARAMS[0] , 
			  INPARAMS.quarks[0] , 
			  fit , LT ) ;

    
    XXX[idx] = x ;
    printf( "%g %g %g \n" , XXX[idx] , Y[idx].avg , Y[idx].err ) ;

    idx++ ;
  }
  */
  //write_distribution_arr2( XXX , Y , "b1.37_m0.001_fit.500.boot" , NMEAS ) ;
