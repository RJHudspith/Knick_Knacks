/**
   @file io_wrapper.c
   @brief wrapper for IO
 */
#include "fitfunc.h"

#include "equivalents.h"
#include "fake_data.h"
#include "GLU_bswap.h"
#include "GLU_ensemble_read.h"
#include "GLUdata_glueprop.h"
#include "GLUdata_config.h"
#include "read_hirep_corr.h"
#include "GLUdata_polycorr.h"
#include "read_simcorrs.h"
#include "read_tetras.h"
#include "read_gevpcorr.h"
#include "read_distribution.h"
#include "read_tmoments.h"
#include "readflow_couple.h"

// call this for fake data
static struct mom_info
zero_mominfo( void )
{
  struct mom_info mominfo ;
  mominfo.p6 = mominfo.p4 = mominfo.p2 = 0.0 ;
  return mominfo ;
}

static struct mom_info **
full_zeromom( double **X , 
	      const int *NDATA ,
	      const int NSLICES )
{
  struct mom_info **mominfo = malloc( NSLICES * sizeof( struct mom_info* ) ) ;
  int i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    mominfo[ i ] = (struct mom_info*)malloc( NDATA[i] * sizeof( struct mom_info ) ) ;
    int j ;
    for( j = 0 ; j < NDATA[i] ; j++ ) {
      mominfo[ i ][ j ] = zero_mominfo( ) ;
      mominfo[ i ][ j ].p2 = X[ i ][ j ] ;
    }
  }
  return mominfo ;
}

int
read_data( double ***X ,
	   struct resampled ***RAW ,
	   struct mom_info ***mominfo ,
	   struct mom_info **moms ,
	   struct input_params *INPARAMS ,
	   int *NSLICES ,
	   const int LT )
{	   
  switch( INPARAMS->type ) {
  case DEFAULT :
  case FAKE :
    printf( "\n--> Fake Data Generation <--\n" ) ;
    *RAW = fake_data( X , (double)LT , INPARAMS -> fittype , INPARAMS -> quarks ,
		     *NSLICES , INPARAMS->NDATA , INPARAMS -> NRAW[0] , LT ,
		     INPARAMS -> sim_params ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case KK_SS_EXTRAP :
  case KK_ML_EXTRAP :
    printf( "\n--> UKhadron Data Reading <--\n" ) ;
    *RAW = read_SUSY( X , INPARAMS , NSLICES ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case UKHADRON :
    printf( "\n--> UKhadron Data Reading <--\n" ) ;
    *RAW = read_SUSY( X , INPARAMS , NSLICES ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case SUPERDIST :
    printf( "\n--> SUPERDIST Reading <--\n" ) ;
    *RAW = read_superdist( mominfo , X , INPARAMS , NSLICES , INPARAMS->mom_type ) ;
    break ;
  case ALPHA_S :
    /*
    printf( "\n--> GLU Data Reading <--\n" ) ;
    *RAW = read_GLUprop( mominfo , X , INPARAMS , NSLICES , INPARAMS -> mom_type ) ;
    */
    printf( "\n--> SUPERDIST Reading <--\n" ) ;
    *RAW = read_superdist( mominfo , X , INPARAMS , NSLICES , INPARAMS->mom_type ) ;

    break ;
  case TAUVUS3 :
    printf( "\n--> Time Moments Reading <--\n" ) ;
    *RAW = read_tmoments( X , mominfo , INPARAMS , NSLICES ) ;
    break ;
  case TAUVUS :
  case TAUVUS2 :
  case PADE_AMU :
  case CONF_AMU :
  case J0 :
    printf( "\n--> SUPERDIST Reading <--\n" ) ;
    *RAW = read_superdist( mominfo , X , INPARAMS , NSLICES , INPARAMS->mom_type ) ;
    break ;
  case FLAVOUR_COMBINATION:
    printf( "\n--> GLU Data Reading <--\n" ) ;
    *RAW = read_GLUprop( mominfo , X , INPARAMS , NSLICES , INPARAMS -> mom_type ) ;
    break ;
  case MASS_SPLITTINGS :
  case DISPERSIONS :
  case SPEED_OF_LIGHT :
  case CORRELATIONS :
    printf( "\n--> Dispersion Reading <--\n" ) ;
    *RAW = read_dispersions( X , INPARAMS , mominfo , NSLICES ) ;
    break ;
  case VV_CORRELATORS :
  case AMA_VV_CORRELATORS :
  case ZV_EVALUATE :
    printf( "\n--> CORR Data Reading <--\n" ) ;
    *RAW = read_simveccorr( X , INPARAMS , moms , NSLICES , INPARAMS->tfold  ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case AA_CORRELATORS :
  case AMA_AA_CORRELATORS :
    printf( "\n--> CORR Data Reading <--\n" ) ;
    *RAW = read_simaxcorr( X , INPARAMS , moms , NSLICES , INPARAMS->tfold  ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case VV_BARYONS :
    printf( "\n--> CORR Data Reading <--\n" ) ;
    *RAW = read_simvecbarcorr( X , INPARAMS , moms , NSLICES , INPARAMS->tfold  ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case PP_BARYONS :
    printf( "\n--> CORR Data Reading <--\n" ) ;
    *RAW = read_simbarcorr( X , INPARAMS , moms , NSLICES , INPARAMS->tfold  ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case PP_CORRELATORS :
  case AMA_PP_CORRELATORS :
    printf( "\n--> CORR Data Reading <--\n" ) ;
    *RAW = read_simcorr( X , INPARAMS , moms , NSLICES , INPARAMS->tfold  ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case CORRELATORS :
    printf( "\n--> CORR Data Reading <--\n" ) ;
    *RAW = read_corr( X , INPARAMS , moms , NSLICES , 2 , 0 , INPARAMS->tfold ) ;
    //*RAW = read_hirep_corr( X , INPARAMS , moms , NSLICES , INPARAMS->tfold ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case TETRA_CORRELATORS :
    printf( "\n--> CORR Data Reading <--\n" ) ;
    *RAW = read_tetras( X , INPARAMS , moms , NSLICES , INPARAMS->tfold  ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case GEVP1_CORRELATORS :
    printf( "\n--> CORR Data Reading <--\n" ) ;
    *RAW = read_gevpcorr( X , INPARAMS , moms , NSLICES , INPARAMS->tfold  ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case GEVP2_CORRELATORS :
    printf( "\n--> CORR Data Reading <--\n" ) ;
    *RAW = read_gevpcorr2( X , INPARAMS , moms , NSLICES , INPARAMS->tfold  ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case TMOMENTS :
    printf( "\n--> Time Moments Reading <--\n" ) ;
    *RAW = read_tmoments( X , mominfo , INPARAMS , NSLICES ) ;
    break ;
  case FLOW_COUPLE :
    printf( "\n--> FLOW couple Reading <--\n" ) ;
    *RAW = read_flowcouple( X , INPARAMS , NSLICES ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  case WILSON_FLOW :
    printf( "Choice Not implemented \n" ) ;
    return -1 ;
  case STATIC_POTENTIAL :
    printf( "\n--> GLU Data Reading <--\n" ) ;
    *RAW = read_GLUprop( mominfo , X , INPARAMS , NSLICES , INPARAMS -> mom_type ) ;
    break ;
  case TOPOLOGICAL_SUSCEPTIBILITY :
    printf( "\n--> GLU Data Reading <--\n" ) ;
    *RAW = read_GLUprop_config( mominfo , X , INPARAMS , NSLICES , INPARAMS -> mom_type ) ;
    //*RAW = read_GLUprop( mominfo , X , INPARAMS , NSLICES , INPARAMS -> mom_type ) ;
    *mominfo = full_zeromom( *X , INPARAMS->NDATA , *NSLICES ) ;
    break ;
  }
  return SUCCESS ;
}
