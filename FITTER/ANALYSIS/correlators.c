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
#include "blackbox.h"

//#define SIGN_FLIP
//#define BLACKBOX_EFFMASS

void
write_corr_data( struct resampled data , 
		 struct mom_info *moms ,
		 const char *name )
{
  // write out the mass to a file
  FILE *outfile = fopen( name , "wb" ) ;

  // write out the momentum
  uint32_t num_mom[1] = { 1 } ;
  if( !BigEndian ) bswap_32( 1 , num_mom ) ;
  fwrite( num_mom , sizeof( uint32_t ) , 1 , outfile ) ;
  uint32_t mom[4] = { 3 , moms[0].n[0] , moms[0].n[1] , moms[0].n[2] } ;
  if( !BigEndian ) bswap_32( 4 , mom ) ;
  fwrite( mom , sizeof( uint32_t ) , 4 , outfile ) ;
 
  // write the distribution and its size
  write_singledist( data , outfile ) ;

  fclose( outfile ) ;
  return ;
}

// decay constants
static void
decay_constant( const struct resampled *fitparams ,
		struct input_params *INPARAMS ,
		const size_t selection ) 
{
  enum { MASS = 0 , 
	 PWPL = 1 , AWAL = 2 , PWAL = 3 , AWPL = 4 , 
	 PWPW = 5 , AWAW = 6 , PWAW = 7 , AWPW = 8 }  ;

  // amplitude helper
  struct resampled res ;
  res.resampled = malloc( fitparams[0].NSAMPLES * sizeof(double) ) ;

  const double Two_VOL = 2.0/ (INPARAMS->dimensions[0][0]*	\
			       INPARAMS->dimensions[0][1]*	\
			       INPARAMS->dimensions[0][2] ) ;
  switch( selection ) {
    // I don't think we can use this, relies on PCAC?
    // f = 1/m * sqrt( 2 * ( PP^{WL} PP^{WL} / ( V m PP^{WW} ) )
  case 0 :
    equate( &res , fitparams[ PWPL ] ) ;
    mult( &res , fitparams[ PWPL ] ) ;
    mult_constant( &res , Two_VOL ) ;
    divide( &res , fitparams[ MASS ] ) ;
    divide( &res , fitparams[ PWPW ] ) ;
    raise( &res , 0.5 ) ;
    divide( &res , fitparams[ MASS ] ) ;
    break ;
    // f = ZA * sqrt( 2 * AP^{WL} PA^{WL} / ( V m PP^{WW} ) )
  case 1 :
    equate( &res , fitparams[ PWAL ] ) ;
    mult( &res , fitparams[ PWAL ] ) ;
    mult_constant( &res , Two_VOL ) ;
    divide( &res , fitparams[ MASS ] ) ;
    divide( &res , fitparams[ PWPW ] ) ;
    raise( &res , 0.5 ) ;
    break ;
  }
  root( &res ) ;
  printf( "Decay^2 :: %e +/- %e \n" , res.avg , res.err ) ;

  free( res.resampled ) ;
  return ;
}

// decay constants
static void
decay_constant2( const struct resampled *fitparams ,
		 struct input_params *INPARAMS ,
		 const size_t selection ) 
{
  // matrix element counter
  enum { MASS = 0 } ;

  // amplitude helper
  struct resampled res ;
  res.resampled = malloc( fitparams[0].NSAMPLES * sizeof(double) ) ;

  const double Two_VOL = 2.0/ (INPARAMS->dimensions[0][0]*	\
			       INPARAMS->dimensions[0][1]*	\
			       INPARAMS->dimensions[0][2] ) ;

  equate( &res , fitparams[ selection ] ) ;
  mult( &res , fitparams[ selection ] ) ;
  mult_constant( &res , Two_VOL ) ;
  divide( &res , fitparams[ MASS ] ) ;
  root( &res ) ;

  printf( "Decay :: %e +/- %e \n" , res.avg , res.err ) ;

  free( res.resampled ) ;
  return ;
}

// compute a matrix element
static void
matrix_element( const struct resampled *fitparams ,
		struct input_params *INPARAMS ,
		const size_t index ,
		const char *message ,
		struct mom_info *moms ) 
{
  // amplitude helper
  struct resampled res ;
  res.resampled = malloc( fitparams[0].NSAMPLES * sizeof(double) ) ;

  const double Two_VOL = 2.0/ (INPARAMS->dimensions[0][0]*	\
			       INPARAMS->dimensions[0][1]*	\
			       INPARAMS->dimensions[0][2] ) ;

  // compute a matrix element as f_1 * sqrt( 2 * mass / L^3 )
  equate( &res , fitparams[ 0 ] ) ;
  mult_constant( &res , Two_VOL ) ;
  root( &res ) ;
  mult( &res , fitparams[ index ] ) ;

  printf( "%s Matrix Element :: %e +/- %e \n" , message , res.avg , res.err ) ;

  /*
  char str[ 256 ] ;
  sprintf( str , "%s.%s.matel" , INPARAMS->output_file , message ) ;
    
  printf( "\n--> Writing out the mass to %s <--\n" , str ) ;
    
  write_corr_data( fitparams[0] , moms , str ) ;
  */

  return ;
}

// tauvus computation
void
correlator_eval( double **xavg ,
		 struct resampled **bootavg ,
		 struct mom_info **mominfo ,
		 struct mom_info *moms ,
		 struct input_params *INPARAMS ,
		 const int NSLICES ,
		 const int LT )
{
#ifdef SIGN_FLIP
  if( INPARAMS -> type == PP_CORRELATORS || INPARAMS -> type == AMA_PP_CORRELATORS ) {
    size_t i ;
    for( i = 0 ; i < NSLICES ; i++ ) {
      size_t j ;
      for( j = 0 ; j < INPARAMS -> NDATA[ i ] ; j++ ) {
	mult_constant( &bootavg[i][j] , -1  ) ;
      }
    }
  }
#endif

  // print effective mass to file
  FILE *file = fopen( "effmass.dat" , "w" ) ;

  // effective masses
  {
#ifdef BLACKBOX_EFFMASS
    size_t i ;
    for( i = 0 ; i < NSLICES ; i++ ) {
      size_t k , j ;
      const size_t NSTATES = 2 ;
      for( k = 0 ; k < NSTATES ; k++ ) {
	struct resampled **effmass = prony_effmass( (const struct resampled**)bootavg ,
						    INPARAMS -> NDATA ,
						    NSLICES , NSTATES , k ) ;
	
	for( j = 0 ; j < INPARAMS -> NDATA[i] ; j++ ) {
	  fprintf( file , "%e %e %e \n" , xavg[i][j] , effmass[i][j].avg , effmass[i][j].err ) ;
	}
	fprintf( file , "\n" ) ;
	// free that mass
	free_resampled_dp( effmass , NSLICES , INPARAMS->NDATA ) ;
      }
    }
#else
    struct resampled **effmass = effective_mass( (const struct resampled**)bootavg ,
						 INPARAMS->NDATA ,
						 NSLICES ,
						 //LOG_EFFMASS ) ;
                                                 //LOG2_EFFMASS ) ;
						 //ACOSH_EFFMASS ) ;
						 //ASINH_EFFMASS ) ;
                                                 ATANH_EFFMASS ) ;
    size_t i , j ;
    for( i = 0 ; i < NSLICES ; i++ ) {
      for( j = 0 ; j < INPARAMS -> NDATA[i] ; j++ ) {
	fprintf( file , "%e %e %e \n" , xavg[i][j] , effmass[i][j].avg , effmass[i][j].err ) ;
      }
      fprintf( file , "\n" ) ;
    }
    free_resampled_dp( effmass , NSLICES , INPARAMS->NDATA ) ;
#endif
  }
  // plot the correlator(s)
  struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bootavg , 
						    (const double**)xavg , 
						    (const struct mom_info **)mominfo ,
						    *INPARAMS , NSLICES , LT ) ;

  if( fitparams != NULL ) {
    int j ;
    for( j = 0 ; j < INPARAMS -> NDATA[0] ; j++ ) {
      if( j >= INPARAMS->fit_lo && j <= INPARAMS->fit_hi ) {
	fprintf( file , "%e %e %e\n" , xavg[0][j] , fitparams[0].err_lo , 0.0 ) ;
      }
    }
    fprintf( file , "\n" ) ;
    for( j = 0 ; j < INPARAMS -> NDATA[0] ; j++ ) {
      if( j >= INPARAMS->fit_lo && j <= INPARAMS->fit_hi ) {
	fprintf( file , "%e %e %e\n" , xavg[0][j] , fitparams[0].avg , 0.0 ) ;
      }
    }
    fprintf( file , "\n" ) ;
    for( j = 0 ; j < INPARAMS -> NDATA[0] ; j++ ) {
      if( j >= INPARAMS->fit_lo && j <= INPARAMS->fit_hi ) {
	fprintf( file , "%e %e %e\n" , xavg[0][j] , fitparams[0].err_hi , 0.0 ) ;
      }
    }

    // compute F_PI using
    if( INPARAMS -> fittype == PPAA ) {
      // fitparam[0] is the mass
      decay_constant( fitparams , INPARAMS , 1 ) ;
    } else if( INPARAMS -> fittype == PPAA_WW ) {
      // is the Local Axial
      decay_constant2( fitparams , INPARAMS , 2 ) ;
      matrix_element( fitparams , INPARAMS , 1 , "Pseudo" , moms ) ;
      matrix_element( fitparams , INPARAMS , 2 , "Axial" , moms ) ;
    } else {
      // is the Local Vector
      decay_constant2( fitparams , INPARAMS , 1 ) ;
      matrix_element( fitparams , INPARAMS , 1 , "Vector" , moms ) ;
    }	       


    // write out the amplitudes as well?
    char str[ 256 ] ;
    sprintf( str , "%s.mass" , INPARAMS->output_file ) ;
    printf( "\n--> Writing out the mass to %s <--\n" , str ) ;
    write_corr_data( fitparams[0] , moms , str ) ;
  }

  fclose( file ) ;

  return ;
}

