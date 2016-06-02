#ifndef FITFUNC_H
#define FITFUNC_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h> // for memcpy
#include <stdint.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

//#define EPSILON_FIT

#define HAVE_SYS_TIME_H

// effective mass types
typedef enum {
  LOG_EFFMASS ,
  LOG2_EFFMASS ,
  ACOSH_EFFMASS , 
  ASINH_EFFMASS ,
  ATANH_EFFMASS } effmass_type ;

// fold selection
typedef enum {
  NOFOLD , PlPl , PlMi , MiPl , MiMi 
} foldselection ;

// add fit types to this as you see fit, hahaha
typedef enum{ NOFIT , EXP , CONST_EXP ,
	      POLY0 , POLY1 , POLY2 , POLY3 , POLY4 , POLY5 ,
	      POLY6 , POLY7 , POLY8 , POLY9 ,
	      PADE11 , 
	      PADE12 , PADE21 ,
	      PADE31 , PADE22 , PADE13 ,
	      PADE41 , PADE23 , PADE32 , PADE14 ,
	      PADE33 ,
	      GLUEFIT , ALPHAS_OPE , D0_FIT , LOG_PLUS_C ,
	      COSH , SINH , GAUSSIAN , BKCHIRAL , 
	      PPAA , PPAA_WW , VV_WW , VV_VTVT_WW ,
              CHIPT_R123_SUSY , CHIPT_R45_SUSY , 
	      CHIPT_B123_SUSY , CHIPT_B45_SUSY ,
	      CHIPT_R145_RENORM , CHIPT_R23_RENORM , 
	      CHIPT_B145_RENORM , CHIPT_B23_RENORM ,
	      CHIPT2_SIM_R23_SUSY , CHIPT2_SIM_R45_SUSY , 
	      CHIPT2_SIM_B23_SUSY , CHIPT2_SIM_B45_SUSY ,
              CHIPT2_SIM_R23_RENORM , CHIPT2_SIM_R45_RENORM ,
	      CHIPT2_SIM_B23_RENORM , CHIPT2_SIM_B45_RENORM ,
	      CHIPT_SIM_R23_SUSY , CHIPT_SIM_R45_SUSY , 
	      CHIPT_SIM_B23_SUSY , CHIPT_SIM_B45_SUSY ,
              CHIPT_SIM_R23_RENORM , CHIPT_SIM_R45_RENORM ,
	      CHIPT_SIM_B23_RENORM , CHIPT_SIM_B45_RENORM } fit_type ;

typedef enum{
  TWOSIN_MOM , PSQ_MOM , SIN_MOM , RSQ } momtype ;

// what analysis are we doing?
typedef enum {
  DEFAULT ,
  ALPHA_S ,
  CORRELATIONS ,
  CORRELATORS , AA_CORRELATORS , VV_CORRELATORS , PP_CORRELATORS ,
  ZV_EVALUATE ,
  TETRA_CORRELATORS , GEVP1_CORRELATORS , GEVP2_CORRELATORS ,
  AMA_PP_CORRELATORS , AMA_AA_CORRELATORS , AMA_VV_CORRELATORS ,
  PP_BARYONS , VV_BARYONS ,
  FAKE ,
  J0 ,
  STATIC_POTENTIAL ,
  TOPOLOGICAL_SUSCEPTIBILITY ,
  WILSON_FLOW , FLOW_COUPLE ,
  UKHADRON ,
  TAUVUS , TAUVUS2 , TAUVUS3 ,
  DISPERSIONS , SPEED_OF_LIGHT , MASS_SPLITTINGS ,
  SUPERDIST ,
  FLAVOUR_COMBINATION ,
  TMOMENTS ,
  PADE_AMU ,
  CONF_AMU ,
  KK_SS_EXTRAP , KK_ML_EXTRAP } analysis_type ;

// enum for momentum average
typedef enum {
  PSQ_AVERAGE ,
  Z_Nm1_AVERAGE ,
  NONE ,
} momavg_type ;

// error types
enum{ ERR , HI , LO , AVE } ;

// SUCCESS and FAILURE for the returns 
#define FAILURE -1
#define SUCCESS 1

// this is the same order as UKhadron
// where did the 1 go?
typedef enum { 
  RAWDATA = 0 ,
  JACKDATA = 2 ,
  BOOTDATA = 3 } resample_type ;

// chiral data stuff ...
struct chiral {
  double m_pi     ;
  double m_k      ;
  double m_eta    ;
  double f_pi     ;
  double f_k      ;
  double ml       ;
  double ms       ;
  double ZA       ;
  double ZV       ;
  double mu       ; // some renormalisation scale
  double ainverse ; // inverse lattice spacing
} ;

// momentum information
struct mom_info {
  int n[ 4 ] ; // number of Fourier modes
  double p2 ;
  double p4 ;
  double p6 ;
  double p8 ;
  double p12 ;
} ;

// wrapper for the fit evaluations
struct x_descriptor {
  double X ;
  struct chiral quark ;
  struct mom_info mom ;
  int LT ;
} ;

// data has space for quarks, need to worry about allocation?
struct data {
  size_t n ;
  const double *X ;
  const double *y ;
  const double *sigma ;
  bool sim_params[ 12 ] ; // hacked in the simparams
  int NPARAMS ;
  double FIT_LO ;
  double FIT_HI ;
  const int *NDATA ; 
  const struct chiral *quarks ;
  struct mom_info *mom ;
  int data_idx ;
  int SIMS ;    // number of simultaneous datasets 
  int NCOMMON ; // number of common parameters
  int LOGICAL_NPARS ; // the number of logical fit parameters
  int LT ;
};

// struct containing our statistics
struct resampled {
  double *resampled ;
  double avg ;
  double err_hi ;
  double err_lo ;
  double err ;
  int NSAMPLES ;
  resample_type restype ;
} ;

// fitfunction callback ...
typedef struct {
  double ( *f ) ( const double *x , const struct x_descriptor X , const int NPARAMS ) ;
  void ( *guess )( double *params , void *data , const int NPARAMS ) ;
  void ( *description )( const char *message , const double *params , const struct x_descriptor X , const int NPARAMS ) ; 
} fitfunc ;

#define MAX_NFILES 25

struct input_params {
  resample_type resample ;
  int NBOOTS ;
  int NDATA[ MAX_NFILES ] ;
  int NRAW[ MAX_NFILES] ;
  fit_type fittype ;
  double fit_hi ;
  double fit_lo ;
  bool sim_params[ 12 ] ; // simultaneous fit parameters
  int binning[ MAX_NFILES ] ;
  int traj_begin[ MAX_NFILES ] ;
  int traj_end[ MAX_NFILES ]   ;
  int traj_increment[ MAX_NFILES ] ;
  int dimensions[ MAX_NFILES ][ MAX_NFILES ] ; // to accommodate for our supergravity simulations
  analysis_type type ;
  char traj_file[ MAX_NFILES ][ 256 ] ; // ten files can be opened at oncea
  long unsigned int Seed ; // RNG seed
  char graph_name[ 128 ] ;
  char graph_xaxis[ 128 ] ;
  char graph_yaxis[ 128 ] ;
  momavg_type momavg ;
  momtype mom_type ;
  bool tfold ;
  struct chiral quarks[ MAX_NFILES ] ; // hard limit on the size of these
  int NFILES ; // number of actual files we have == number of actual chiral dat
  char output_file[ 128 ] ;
  int NCHANNELS ;
} ;

// for the endianess filter ...
int BigEndian ;

// for the operations on the resampled data
#include "resampled_ops.h"
#include "stats.h"

#endif
