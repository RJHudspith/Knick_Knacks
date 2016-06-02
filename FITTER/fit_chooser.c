#include "fitfunc.h"
#include "fit_chooser.h" // for the enum

// fit types
#include "expfit.h"
#include "const_exp.h"
#include "multiexp.h"
#include "padefit.h"
#include "polyfit.h"
#include "gluefit.h"
#include "prop.h"
#include "log_plus_c.h"
#include "coshfit.h"
#include "sinhfit.h"
#include "gaussian.h"

// found in ./ALPHAS/
#include "D0_diff.h"
#include "D0_diff_v1.h"
#include "D0_diff_multi.h"

// found in ./STATPOT/
#include "StaticPotfit.h"
#include "polystatfit.h"
#include "polystatalpha.h"

// found in ./CHIRAL_SUSY/
#include "chiral_continuum.h"
#include "chiral_continuum_R23.h"
#include "chiral_R23.h"
#include "chiral_R23_sim.h"

// found in ./CORRELATORS/
#include "pp_ap_pa_aa.h" 
#include "pp_ap_pa_aa_WW.h" 
#include "vv_WW.h" 
#include "vv_vtvt_WW.h" 

static gsl_multifit_function_fdf
init_pade_f( const fit_type FITTYPE , 
	     fitfunc *fit , 
	     int *NPARAMS , 
	     int *NCOMMON )
{
  gsl_multifit_function_fdf tmp ;
  // standard pade fit options
  tmp.f   = pade_f ;
  tmp.df  = pade_df ;
  tmp.fdf = pade_fdf ;
  fit->f     = pade_eval ;
  fit->guess = pade_guess ;
  fit->description = pade_description ;
  switch( FITTYPE ) {
  case PADE11 : pade_set_nm( 1 , 1 ) ; *NPARAMS = 3 ; break ;
  case PADE12 : pade_set_nm( 1 , 2 ) ; *NPARAMS = 4 ; break ;
  case PADE21 : pade_set_nm( 2 , 1 ) ; *NPARAMS = 4 ; break ;
  case PADE31 : pade_set_nm( 3 , 1 ) ; *NPARAMS = 5 ; break ;
  case PADE22 : pade_set_nm( 2 , 2 ) ; *NPARAMS = 5 ; break ;
  case PADE13 : pade_set_nm( 1 , 3 ) ; *NPARAMS = 5 ; break ;
  case PADE41 : pade_set_nm( 4 , 1 ) ; *NPARAMS = 6 ; break ;
  case PADE32 : pade_set_nm( 3 , 2 ) ; *NPARAMS = 6 ; break ; 
  case PADE23 : pade_set_nm( 2 , 3 ) ; *NPARAMS = 6 ; break ;
  case PADE14 : pade_set_nm( 1 , 4 ) ; *NPARAMS = 6 ; break ;
  case PADE33 : pade_set_nm( 3 , 3 ) ; *NPARAMS = 7 ; break ;
  default :
    pade_set_nm( 1 , 1 ) ; *NPARAMS = 3 ; break ;
    break ;
  }
  return tmp ;
}

static gsl_multifit_function_fdf
init_poly_f( const fit_type FITTYPE , 
	     fitfunc *fit , 
	     int *NPARAMS , 
	     int *NCOMMON )
{
  gsl_multifit_function_fdf tmp ;
  // standard pade fit options
  tmp.f   = poly_f ;
  tmp.df  = poly_df ;
  tmp.fdf = poly_fdf ;
  fit->f     = poly_eval ;
  fit->guess = poly_guess ;
  fit->description = poly_description ;
  switch( FITTYPE ) {
  case POLY9 : *NPARAMS = 9 ; break ; 
  case POLY8 : *NPARAMS = 8 ; break ;
  case POLY7 : *NPARAMS = 7 ; break ;
  case POLY6 : *NPARAMS = 6 ; break ;
  case POLY5 : *NPARAMS = 6 ; break ;
  case POLY4 : *NPARAMS = 5 ; break ;
  case POLY3 : *NPARAMS = 4 ; break ;
  case POLY2 : *NPARAMS = 3 ; break ;
  case POLY1 : *NPARAMS = 2 ; break ;
  case POLY0 : *NPARAMS = 1 ; break ;
  default :
    *NPARAMS = 2 ; break ;
  }
  return tmp ;
}

// susy bk chiral perturbation theory fits
static gsl_multifit_function_fdf
init_chipt_f( const fit_type FITTYPE , 
	      fitfunc *fit , 
	      int *NPARAMS , 
	      int *NCOMMON )
{ 
  gsl_multifit_function_fdf tmp ;
  tmp.f   = chiral_cont_R23_f ;
  tmp.df  = chiral_cont_R23_df ;
  tmp.fdf = chiral_cont_R23_fdf ;
  fit->f     = chiral_cont_R23_eval ;
  fit->guess = chiral_cont_R23_guess ;
  fit->description = chiral_cont_R23_description ;
  *NPARAMS = 3 ;
  switch( FITTYPE ) {
    // RENORM RATIOS
  case CHIPT_R145_RENORM : chiral_cont_R23_setfac(  3./2. ) ; break ;
  case CHIPT_R23_RENORM : chiral_cont_R23_setfac(  5./2. ) ; break ;
    // RENORM BS
  case CHIPT_B145_RENORM : chiral_cont_R23_setfac( -1./2. ) ; break ;
  case CHIPT_B23_RENORM : chiral_cont_R23_setfac(  1./2. ) ; break ;
    // SUSY RS
  case CHIPT_R123_SUSY :  chiral_cont_R23_setfac(  3./2. ) ; break ;
  case CHIPT_R45_SUSY :  chiral_cont_R23_setfac(  5./2. ) ; break ;
    // SUSY BS
  case CHIPT_B123_SUSY : chiral_cont_R23_setfac( -1./2. ) ; break ;
  case CHIPT_B45_SUSY :  chiral_cont_R23_setfac(  1./2. ) ; break ;
  default :
    chiral_cont_R23_setfac(  1./2. ) ;
    break ;
  }
  return tmp ;
}

// susy bk chiral perturbation theory fits
static gsl_multifit_function_fdf
init_chipt2_sim_f( const fit_type FITTYPE , 
		  fitfunc *fit , 
		  int *NPARAMS , 
		  int *NCOMMON )
{ 
  gsl_multifit_function_fdf tmp ;
  tmp.f   = chiral_R23_f ;
  tmp.df  = chiral_R23_df ;
  tmp.fdf = chiral_R23_fdf ;
  fit->f     = chiral_R23_eval ;
  fit->guess = chiral_R23_guess ;
  fit->description = chiral_R23_description ;
  *NPARAMS = 2 ;
  switch( FITTYPE ) {
    // RENORM RATIOS
  case CHIPT2_SIM_R23_RENORM : chiral_R23_setfac( 5./2. ) ; break ;
  case CHIPT2_SIM_R45_RENORM : chiral_R23_setfac( 3./2. ) ; break ;
    // RENORM BS
  case CHIPT2_SIM_B23_RENORM : chiral_R23_setfac( 1./2. ) ; break ;
  case CHIPT2_SIM_B45_RENORM : chiral_R23_setfac( -1./2. ) ; break ;
    // SUSY RS
  case CHIPT2_SIM_R23_SUSY : chiral_R23_setfac( 3./2. ) ; break ;
  case CHIPT2_SIM_R45_SUSY : chiral_R23_setfac( 5./2. ) ; break ;
    // SUSY BS
  case CHIPT2_SIM_B23_SUSY : chiral_R23_setfac( -1./2. ) ; break ;
  case CHIPT2_SIM_B45_SUSY : chiral_R23_setfac( 1./2. ) ; break ;
  default :
    printf( "WHAT? \n" ) ; exit(1) ;
    break ;
  }
  return tmp ;
}

// susy bk chiral perturbation theory fits
static gsl_multifit_function_fdf
init_chipt_sim_f( const fit_type FITTYPE , 
		  fitfunc *fit , 
		  int *NPARAMS , 
		  int *NCOMMON )
{ 
  gsl_multifit_function_fdf tmp ;
  tmp.f   = chiral_R23_sim_f ;
  tmp.df  = chiral_R23_sim_df ;
  tmp.fdf = chiral_R23_sim_fdf ;
  fit->f     = chiral_R23_sim_eval ;
  fit->guess = chiral_R23_sim_guess ;
  fit->description = chiral_R23_sim_description ;
  *NPARAMS = 9 ; 
  switch( FITTYPE ) {
    // RENORM RATIOS
  case CHIPT_SIM_R23_RENORM : chiral_R23_sim_setfac( 5./2. ) ; break ;
  case CHIPT_SIM_R45_RENORM : chiral_R23_sim_setfac( 3./2. ) ; break ;
    // RENORM BS
  case CHIPT_SIM_B23_RENORM : chiral_R23_sim_setfac( 1./2. ) ; break ;
  case CHIPT_SIM_B45_RENORM : chiral_R23_sim_setfac( -1./2. ) ; break ;
    // SUSY RS
  case CHIPT_SIM_R23_SUSY : chiral_R23_sim_setfac( 3./2. ) ; break ;
  case CHIPT_SIM_R45_SUSY : chiral_R23_sim_setfac( 5./2. ) ; break ;
    // SUSY BS
  case CHIPT_SIM_B23_SUSY : chiral_R23_sim_setfac( -1./2. ) ; break ;
  case CHIPT_SIM_B45_SUSY : chiral_R23_sim_setfac( 1./2. ) ; break ;
  default :
    printf( "WHAT? \n" ) ; exit(1) ;
    break ;
  }
  return tmp ;
}

gsl_multifit_function_fdf
initialise_f( const fit_type FITTYPE , 
	      fitfunc *fit , 
	      int *NPARAMS , 
	      int *NCOMMON )
{
  gsl_multifit_function_fdf tmp ;

  // standard polynomial fit
  tmp.f   = poly_f ;
  tmp.df  = poly_df ;
  tmp.fdf = poly_fdf ;
  fit->f     = poly_eval ;
  fit->guess = poly_guess ;
  fit->description = poly_description ;

  // switch for the fit
  switch( FITTYPE ) {
  case D0_FIT : 
    /*
    tmp.f   = D0_diff_f ;
    tmp.df  = D0_diff_df ;
    tmp.fdf = D0_diff_fdf ;
    fit->f     = D0_diff_eval ;
    fit->guess = D0_diff_guess ;
    fit->description = D0_diff_description ;
    *NPARAMS = D0_diff_n( ) ;
    */
    tmp.f   = D0_diff_v1_f ;
    tmp.df  = D0_diff_v1_df ;
    tmp.fdf = D0_diff_v1_fdf ;
    fit->f     = D0_diff_v1_eval ;
    fit->guess = D0_diff_v1_guess ;
    fit->description = D0_diff_v1_description ;
    *NPARAMS = D0_diff_v1_n( ) ;
    /*
    tmp.f   = D0_diff_multi_f ;
    tmp.df  = D0_diff_multi_df ;
    tmp.fdf = D0_diff_multi_fdf ;
    fit->f     = D0_diff_multi_eval ;
    fit->guess = D0_diff_multi_guess ;
    fit->description = D0_diff_multi_description ;
    *NPARAMS = D0_diff_multi_n( ) ;
    */
    break ;
  case COSH :
    tmp.f   = cosh_f ;
    tmp.df  = cosh_df ;
    tmp.fdf = cosh_fdf ;
    fit->f     = cosh_eval ;
    fit->guess = cosh_guess ;
    fit->description = cosh_description ;
    *NPARAMS = 2 ;
    break ;
  case PPAA :
    tmp.f   = PPAA_f ;
    tmp.df  = PPAA_df ;
    tmp.fdf = PPAA_fdf ;
    fit->f     = PPAAc_eval ;
    fit->guess = PPAA_guess ;
    fit->description = PPAA_description ;
    *NPARAMS = 2 ;
    break ;
  case PPAA_WW :
    tmp.f   = PPAA_WW_f ;
    tmp.df  = PPAA_WW_df ;
    tmp.fdf = PPAA_WW_fdf ;
    fit->f     = PPAA_WW_eval ;
    fit->guess = PPAA_WW_guess ;
    fit->description = PPAA_WW_description ;
    *NPARAMS = 5 ;
    break ;
  case VV_WW :
    tmp.f   = VV_WW_f ;
    tmp.df  = VV_WW_df ;
    tmp.fdf = VV_WW_fdf ;
    fit->f     = VV_WW_eval ;
    fit->guess = VV_WW_guess ;
    fit->description = VV_WW_description ;
    *NPARAMS = 3 ;
    break ;
  case VV_VTVT_WW :
    tmp.f   = VV_VTVT_WW_f ;
    tmp.df  = VV_VTVT_WW_df ;
    tmp.fdf = VV_VTVT_WW_fdf ;
    fit->f     = VV_VTVT_WW_eval ;
    fit->guess = VV_VTVT_WW_guess ;
    fit->description = VV_VTVT_WW_description ;
    *NPARAMS = 5 ;
    break ;
  case SINH :
    tmp.f   = sinh_f ;
    tmp.df  = sinh_df ;
    tmp.fdf = sinh_fdf ;
    fit->f     = sinh_eval ;
    fit->guess = sinh_guess ;
    fit->description = sinh_description ;
    *NPARAMS = 2 ;
    break ;
  case EXP :
    tmp.f   = exp_f ;
    tmp.df  = exp_df ;
    tmp.fdf = exp_fdf ;
    fit->f     = exp_eval ;
    fit->guess = exp_guess ;
    fit->description = exp_description ;
    *NPARAMS = 2 ;
    break ;
  case CONST_EXP :
    tmp.f   = const_exp_f ;
    tmp.df  = const_exp_df ;
    tmp.fdf = const_exp_fdf ;
    fit->f     = const_exp_eval ;
    fit->guess = const_exp_guess ;
    fit->description = const_exp_description ;
    *NPARAMS = 3 ;
    break ;
  case GAUSSIAN :
    tmp.f   = gaussian_f ;
    tmp.df  = gaussian_df ;
    tmp.fdf = gaussian_fdf ;
    fit->f     = gaussian_eval ;
    fit->guess = gaussian_guess ;
    fit->description = gaussian_description ;
    *NPARAMS = 2 ;
    break ;
  case PADE11 : case PADE12 : case PADE21 : case PADE31 :
  case PADE22 : case PADE13 : case PADE41 : case PADE32 :
  case PADE23 : case PADE14 : case PADE33 :
    return init_pade_f( FITTYPE , fit , NPARAMS , NCOMMON ) ;
  case POLY9 : case POLY8 : case POLY7 : case POLY6 :
  case POLY5 : case POLY4 : case POLY3 : case POLY2 :
  case POLY1 : case POLY0 :
    return init_poly_f( FITTYPE , fit , NPARAMS , NCOMMON ) ;
  case GLUEFIT :
    tmp.f   = polystatalpha_f ;
    tmp.df  = polystatalpha_df ;
    tmp.fdf = polystatalpha_fdf ;
    fit->f     = polystatalpha_eval ;
    fit->guess = polystatalpha_guess ;
    fit->description = polystatalpha_description ;
    *NCOMMON = 2 ;
    *NPARAMS = 3 ;
    break ;
  case LOG_PLUS_C :
    tmp.f   = log_plus_c_f ;
    tmp.df  = log_plus_c_df ;
    tmp.fdf = log_plus_c_fdf ;
    fit->f     = log_plus_c_eval ;
    fit->guess = log_plus_c_guess ;
    *NPARAMS = 2 ;
    break ;
  case BKCHIRAL :
    tmp.f   = chiral_cont_f ;
    tmp.df  = chiral_cont_df ;
    tmp.fdf = chiral_cont_fdf ;
    fit->f     = chiral_cont_eval ;
    fit->guess = chiral_cont_guess ;
    fit->description = chiral_cont_description ;
    *NPARAMS = 3 ;
    break ;
    // chiral fits
  case CHIPT_R145_RENORM : case CHIPT_R23_RENORM : 
  case CHIPT_B145_RENORM : case CHIPT_B23_RENORM : 
  case CHIPT_R123_SUSY : case CHIPT_R45_SUSY : 
  case CHIPT_B123_SUSY : case CHIPT_B45_SUSY : 
    return init_chipt_f( FITTYPE , fit , NPARAMS , NCOMMON ) ;
  case CHIPT2_SIM_R45_RENORM : case CHIPT2_SIM_R23_RENORM : 
  case CHIPT2_SIM_B45_RENORM : case CHIPT2_SIM_B23_RENORM : 
  case CHIPT2_SIM_R23_SUSY : case CHIPT2_SIM_R45_SUSY : 
  case CHIPT2_SIM_B23_SUSY : case CHIPT2_SIM_B45_SUSY :
    return init_chipt2_sim_f( FITTYPE , fit , NPARAMS , NCOMMON ) ;
  case CHIPT_SIM_R45_RENORM : case CHIPT_SIM_R23_RENORM : 
  case CHIPT_SIM_B45_RENORM : case CHIPT_SIM_B23_RENORM : 
  case CHIPT_SIM_R23_SUSY : case CHIPT_SIM_R45_SUSY : 
  case CHIPT_SIM_B23_SUSY : case CHIPT_SIM_B45_SUSY :
    return init_chipt_sim_f( FITTYPE , fit , NPARAMS , NCOMMON ) ;
  default :
    printf( "\n--> fit_chooser.c <--\n\n" ) ;
    printf( "Fit type %d not recognised \n" , FITTYPE ) ;
    return tmp ;
  }
  // sanity check
  if( *NPARAMS < *NCOMMON ) {
    printf( "\n--> fit_chooser.c <--\n\n" ) ;
    printf( "There are more (%d) Sim. params than fit params (%d) ... Leaving\n" ,
	    *NCOMMON , *NPARAMS ) ;
    exit(1) ;
  }
  return tmp ;
}
