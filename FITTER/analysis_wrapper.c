/**
   @file analysis_wrapper.c
   @brief wrapper for the analysis part of the code
 */
#include "fitfunc.h"

// analysis codes
#include "fit_and_plot.h"
#include "wflow_analysis.h"
#include "alphas_OPE.h"
#include "j0contrib.h"
#include "tauvus.h"
#include "tauvus_2.h"
#include "tauvus_3.h"
#include "dispersions.h"
#include "flavour_combination.h"
#include "correlators.h"
#include "speed_of_light.h"
#include "generalised_evp.h"
#include "generalised_evp2.h"
#include "time_moments.h"
#include "pade_amu.h"
#include "conformal_amu.h"
#include "bk_ss_extrap.h"
#include "bk_ll_extrap.h"
#include "mass_splittings.h"
#include "gradflow_inv.h"
#include "Adler.h"
#include "alphas_sub.h"
#include "arbitrary_mom.h"
#include "statpot.h"
#include "baryon_del_M.h"
#include "alphas_ope_v1.h"
#include "mass_extrap.h"
#include "mass_extrap2.h"
#include "ZV_evaluate.h"
#include "CORRelation.h"
#include "ama.h"
#include "tetra_eval.h"
#include "baryon_eval.h"
#include "b_correlators.h"
#include "tetra_eval_meson.h"

int
perform_analysis( double **xavg ,
		  struct resampled **BOOT , 
		  struct mom_info **mominfoavg ,
		  struct mom_info *moms ,
		  struct input_params *INPARAMS ,
		  const int NSLICES ,
		  const int LT )
{
  switch( INPARAMS -> type ) {
  case ALPHA_S :
    //alphas_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT , false ) ;
    alphas_v1_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT , false ) ;
    //alphas_sub( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT , false ) ;
    //adler_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT , false ) ;
    break ;
  case J0 :
    j0_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT , true ) ;
    break ;
  case FLAVOUR_COMBINATION :
    flavour_combination_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT , true ) ;
    break ;
  case TAUVUS :
    tauvus_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT , true ) ;
    break ;
  case TAUVUS2 :
    tauvus2_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT , true ) ;
    break ;
  case TAUVUS3 :
    tauvus3_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT , true ) ;
    break ;
  case MASS_SPLITTINGS :
    mass_splittings_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT ) ;
    //mass_extrap( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT ) ;
    //mass_extrap2( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT ) ;
    break ;
  case DISPERSIONS :
    dispersions_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT ) ;
    break ;
  case SPEED_OF_LIGHT :
    speed_of_light_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT ) ;
    break ;
  case GEVP1_CORRELATORS :
    gevp_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , 2 , LT ) ;
    break ;
  case GEVP2_CORRELATORS :
    gevp2_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , 2 , LT ) ;
    break ;
  case ZV_EVALUATE :
    ZV_eval( xavg , BOOT , mominfoavg , moms , INPARAMS , NSLICES , LT ) ;
    break ;
    // AMAS
  case AMA_AA_CORRELATORS :
  case AMA_VV_CORRELATORS :
  case AMA_PP_CORRELATORS :
    ama_eval( xavg , BOOT , mominfoavg , moms , INPARAMS , NSLICES , LT ) ;
    break ;
  case CORRELATORS :
  case AA_CORRELATORS :
  case PP_CORRELATORS :
  case VV_CORRELATORS :
    correlator_eval( xavg , BOOT , mominfoavg , moms , INPARAMS , NSLICES , LT ) ;
    //bmeson_eval( xavg , BOOT , mominfoavg , moms , INPARAMS , NSLICES , LT ) ;
    break ;
  case PP_BARYONS :
  case VV_BARYONS :
    baryon_eval( xavg , BOOT , mominfoavg , moms , INPARAMS , NSLICES , LT ) ;
    break ;
  case TETRA_FULL :
    tetra_eval( xavg , BOOT , mominfoavg , moms , INPARAMS , NSLICES , LT ) ;
    break ;
  case TETRA_MESONS :
    tetra_eval_meson( xavg , BOOT , mominfoavg , moms , INPARAMS , NSLICES , LT ) ;
    break ;
  case TETRA_CORRELATORS :
    tetra_eval_corr( xavg , BOOT , mominfoavg , moms , INPARAMS , NSLICES , LT ) ;
    break ;
  case FLOW_COUPLE :
    gradflow_inv( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT ) ;
    break ;
  case FAKE :
  case DEFAULT :
  case SUPERDIST :
  case UKHADRON :
    break ;
  case STATIC_POTENTIAL :
    statpot_eval( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT ) ;
    break ;
  case KK_ML_EXTRAP :
    ll_extrap_eval( xavg , BOOT , mominfoavg , moms , 
		    INPARAMS , NSLICES , LT ) ;
    break ;
  case KK_SS_EXTRAP :
    ss_extrap_eval( xavg , BOOT , mominfoavg , moms , 
		    INPARAMS , NSLICES , LT ) ;
    break ;
  case TMOMENTS :
    /*
    time_moment_eval( xavg , BOOT , mominfoavg ,
    		      INPARAMS , NSLICES , LT , false ) ;
    */
    arbitrary_mom_eval( xavg , BOOT , mominfoavg ,
			INPARAMS , NSLICES , LT , false ) ;
    break ;
  case PADE_AMU :
    pade_amu_eval( xavg , BOOT , mominfoavg ,
		   INPARAMS , NSLICES , LT , false ) ;
    break ;
  case CONF_AMU :
    conformal_amu_eval( xavg , BOOT , mominfoavg ,
			INPARAMS , NSLICES , LT , false ) ;
    break ;
  case WILSON_FLOW :
    break ;
  case TOPOLOGICAL_SUSCEPTIBILITY :
    correlator_eval( xavg , BOOT , mominfoavg , moms , INPARAMS , NSLICES , LT ) ;
    break ;
  case CORRELATIONS :
    correlation( xavg , BOOT , mominfoavg , INPARAMS , NSLICES , LT ) ;
    break ;
  }
  return SUCCESS ;
}
