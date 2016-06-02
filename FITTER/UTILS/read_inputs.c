// stupid input file parser

#include "fitfunc.h"
#include "fit_chooser.h"

// we might want to change this at some point
#define STR_LENGTH 256

#define FAILURE -1

// tokenize the input file
struct inputs {
  char TOKEN[ STR_LENGTH ] ;
  char VALUE[ STR_LENGTH ] ;
} ;

static struct inputs *INPUT ;

// counter for the number of tags
static int NTAGS = 0 ;

// allocate the input file struct
void 
pack_inputs( FILE *setup )
{
  INPUT = ( struct inputs* ) malloc ( STR_LENGTH * sizeof( struct inputs ) ) ;
  while( NTAGS++ , fscanf( setup , "%s = %s" , INPUT[ NTAGS ].TOKEN , INPUT[ NTAGS ].VALUE )  != EOF ) { }
  return ;
}

// frees up the struct
void
unpack_inputs( void )
{
  free( INPUT ) ;
  return ;
}

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) 
{
  return ( strcmp( str_1 , str_2 ) != 0 ) ? 0 : 1 ;
}

// look for a tag and return the index if found, else return -1
static int
tag_search( const char *tag ) 
{
  int i ;
  for( i = 0 ; i < NTAGS ; i++ ) {
    if( are_equal( INPUT[i].TOKEN , tag ) ) return i ;
  }
  return FAILURE ;
}

// prints out the problematic tag
static int
tag_failure( const char *tag )
{
  printf( "[IO] Failure looking for tag %s in input file ... Leaving\n" , tag ) ;
  return FAILURE ;
}

// I use this all over the place ...
static char*
set_token( char *search_term )
{
  // graph name
  const int idx = tag_search( search_term ) ;
  if( idx == FAILURE ) { tag_failure( search_term ) ; exit(-1) ; }
  return INPUT[idx].VALUE ;
}

static fit_type
get_fittype( void )
{
  fit_type fittype = NOFIT ;
  // choose the fittype
  const int fittype_idx = tag_search( "FITTYPE" ) ;
  if( fittype_idx == FAILURE ) { return tag_failure( "FITTYPE" ) ; }
  if( are_equal( INPUT[fittype_idx].VALUE , "EXP" ) ) {
    fittype = EXP ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CONST_EXP" ) ) {
    fittype = CONST_EXP ;
    // polynomial options start here
  } else if( are_equal( INPUT[fittype_idx].VALUE , "POLY0") ) {
    fittype = POLY0 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "POLY1") ) {
    fittype = POLY1 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "POLY2") ) {
    fittype = POLY2 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "POLY3") ) {
    fittype = POLY3 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "POLY4") ) {
    fittype = POLY4 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "POLY5") ) {
    fittype = POLY5 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "POLY6") ) {
    fittype = POLY6 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "POLY7") ) {
    fittype = POLY7 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "POLY8") ) {
    fittype = POLY8 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "POLY9") ) {
    fittype = POLY9 ;
    // pade options start here
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PADE11") ) {
    fittype = PADE11 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PADE21") ) {
    fittype = PADE21 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PADE12") ) {
    fittype = PADE12 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PADE31") ) {
    fittype = PADE31 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PADE22") ) {
    fittype = PADE22 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PADE13") ) {
    fittype = PADE13 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PADE41") ) {
    fittype = PADE41 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PADE32") ) {
    fittype = PADE32 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PADE23") ) {
    fittype = PADE23 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PADE14") ) {
    fittype = PADE14 ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PADE33") ) {
    fittype = PADE33 ;
    // pade options end here
  } else if( are_equal( INPUT[fittype_idx].VALUE , "GLUEFIT") ) {
    fittype = GLUEFIT ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "ALPHAS_OPE") ) {
    fittype = ALPHAS_OPE ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "D0_FIT") ) {
    fittype = D0_FIT ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "LOG_PLUS_C") ) {
    fittype = LOG_PLUS_C ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "COSH") ) {
    fittype = COSH ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PPAA") ) {
    fittype = PPAA ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "PPAA_WW") ) {
    fittype = PPAA_WW ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "VV_WW") ) {
    fittype = VV_WW ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "VV_VTVT_WW") ) {
    fittype = VV_VTVT_WW ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "SINH") ) {
    fittype = SINH ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "GAUSSIAN") ) {
    fittype = GAUSSIAN ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "BKCHIRAL") ) {
    fittype = BKCHIRAL ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_R123_SUSY") ) {
    fittype = CHIPT_R123_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_R45_SUSY") ) {
    fittype = CHIPT_R45_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_B123_SUSY") ) {
    fittype = CHIPT_B123_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_B45_SUSY") ) {
    fittype = CHIPT_B45_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_R145_RENORM") ) {
    fittype = CHIPT_R145_RENORM ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_R23_RENORM") ) {
    fittype = CHIPT_R23_RENORM ;
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_B145_RENORM") ) {
    fittype = CHIPT_B145_RENORM ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_B23_RENORM") ) {
    fittype = CHIPT_B23_RENORM ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT2_SIM_R23_SUSY") ) {
    fittype = CHIPT2_SIM_R23_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT2_SIM_R45_SUSY") ) {
    fittype = CHIPT2_SIM_R45_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT2_SIM_B23_SUSY") ) {
    fittype = CHIPT2_SIM_B23_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT2_SIM_B45_SUSY") ) {
    fittype = CHIPT2_SIM_B45_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT2_SIM_R23_RENORM") ) {
    fittype = CHIPT2_SIM_R23_RENORM ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT2_SIM_R45_RENORM") ) {
    fittype = CHIPT2_SIM_R45_RENORM ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT2_SIM_B23_RENORM") ) {
    fittype = CHIPT2_SIM_B23_RENORM ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT2_SIM_B45_RENORM") ) {
    fittype = CHIPT2_SIM_B45_RENORM ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_SIM_R23_SUSY") ) {
    fittype = CHIPT_SIM_R23_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_SIM_R45_SUSY") ) {
    fittype = CHIPT_SIM_R45_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_SIM_B23_SUSY") ) {
    fittype = CHIPT_SIM_B23_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_SIM_B45_SUSY") ) {
    fittype = CHIPT_SIM_B45_SUSY ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_SIM_R23_RENORM") ) {
    fittype = CHIPT_SIM_R23_RENORM ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_SIM_R45_RENORM") ) {
    fittype = CHIPT_SIM_R45_RENORM ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_SIM_B23_RENORM") ) {
    fittype = CHIPT_SIM_B23_RENORM ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "CHIPT_SIM_B45_RENORM") ) {
    fittype = CHIPT_SIM_B45_RENORM ; 
  } else if( are_equal( INPUT[fittype_idx].VALUE , "NOFIT") ) {
    fittype = NOFIT ; 
  } else {
    printf( "FITTYPE %s not recognised \n" , INPUT[fittype_idx].VALUE ) ;
    printf( "Not performing a fit\n" ) ;
    return fittype ;
  }
  printf( "Attempting to fit %s \n" , INPUT[fittype_idx].VALUE ) ;
  return fittype ;
}

static void
set_chiral( struct input_params *INPARAMS )
{
  // read the chiral data
  int counter = 0 ;
  while( counter < MAX_NFILES ) {
    char str[ 16 ] ;
    int chicount = 0 ;
    //////////////////////////////////////////////////////////////////
    sprintf( str , "MPI_%d" , counter ) ;
    const int mpi_idx = tag_search( str ) ;
    if( mpi_idx != FAILURE ) chicount++ ;
    INPARAMS -> quarks[ counter ].m_pi = atof( INPUT[mpi_idx].VALUE ) ;
    //
    sprintf( str , "MU_%d" , counter ) ;
    const int mu_idx = tag_search( str ) ;
    if( mu_idx != FAILURE ) chicount++ ;
    INPARAMS -> quarks[ counter ].mu = atof( INPUT[mu_idx].VALUE ) ;
    //////////////////////////////////////////////////////////////////
    sprintf( str , "MK_%d" , counter ) ;
    const int mk_idx = tag_search( str ) ;
    if( mk_idx != FAILURE ) chicount++ ;
    INPARAMS -> quarks[ counter ].m_k = atof( INPUT[mk_idx].VALUE ) ;
    //////////////////////////////////////////////////////////////////
    sprintf( str , "META_%d" , counter ) ;
    const int meta_idx = tag_search( str ) ;
    if( meta_idx != FAILURE ) chicount++ ;
    INPARAMS -> quarks[ counter ].m_eta = atof( INPUT[meta_idx].VALUE ) ;
    //////////////////////////////////////////////////////////////////
    sprintf( str , "FPI_%d" , counter ) ;
    const int fpi_idx = tag_search( str ) ;
    if( fpi_idx != FAILURE ) chicount++ ;
    INPARAMS -> quarks[ counter ].f_pi = atof( INPUT[fpi_idx].VALUE ) ;
    //////////////////////////////////////////////////////////////////
    sprintf( str , "FK_%d" , counter ) ;
    const int fk_idx = tag_search( str ) ;
    if( fk_idx != FAILURE ) chicount++ ;
    INPARAMS -> quarks[ counter ].f_k = atof( INPUT[fk_idx].VALUE ) ;
    //////////////////////////////////////////////////////////////////
    sprintf( str , "ML_%d" , counter ) ;
    const int ml_idx = tag_search( str ) ;
    if( ml_idx != FAILURE ) chicount++ ;
    INPARAMS -> quarks[ counter ].ml = atof( INPUT[ml_idx].VALUE ) ;
    //////////////////////////////////////////////////////////////////
    sprintf( str , "MS_%d" , counter ) ;
    const int ms_idx = tag_search( str ) ;
    if( ms_idx != FAILURE ) chicount++ ;
    INPARAMS -> quarks[ counter ].ms = atof( INPUT[ms_idx].VALUE ) ;
    //////////////////////////////////////////////////////////////////
    sprintf( str , "ZA_%d" , counter ) ;
    const int ZA_idx = tag_search( str ) ;
    if( ZA_idx != FAILURE ) chicount++ ;
    INPARAMS -> quarks[ counter ].ZA = atof( INPUT[ZA_idx].VALUE ) ;
    //////////////////////////////////////////////////////////////////
    sprintf( str , "ZV_%d" , counter ) ;
    const int ZV_idx = tag_search( str ) ;
    if( ZV_idx != FAILURE ) chicount++ ;
    INPARAMS -> quarks[ counter ].ZV = atof( INPUT[ZV_idx].VALUE ) ;
    //////////////////////////////////////////////////////////////////
    sprintf( str , "ainverse_%d" , counter ) ;
    const int ainv_idx = tag_search( str ) ;
    if( ainv_idx != FAILURE ) chicount++ ;
    INPARAMS -> quarks[ counter ].ainverse = atof( INPUT[ainv_idx].VALUE ) ;
    //////////////////////////////////////////////////////////////////
    if( chicount == 0 ) break ;
    counter++ ;
  }
  return ;
}

// set the simultaneous parameters
static void
set_simultaneous( struct input_params *INPARAMS )
{
  // read the chiral data
  int counter = 0 ;
  for( counter = 0 ; counter < 12 ; counter++ ) {
    INPARAMS -> sim_params[ counter ] = false ;
  }

  char str[ 256 ] ;
  counter = 0 ;
  while( counter < 12 ) {
    sprintf( str , "SIMPARAM_%d" , counter ) ;
    const int sim_idx = tag_search( str ) ;
    if( sim_idx == FAILURE ) break ;
    INPARAMS -> sim_params[ atoi( INPUT[sim_idx].VALUE ) ] = true ;
    printf( "Simultaneous fit index %d set\n" , atoi( INPUT[sim_idx].VALUE ) ) ;
    counter++ ;
  }
}

// set the simulation dimensions
static void
set_dimensions( struct input_params *INPARAMS )
{
  // read lattice dimensions
  int file_counter = 0 ;
  while( file_counter < MAX_NFILES ) {
    int mu = 0 ;
    for( mu = 0 ; mu < 4 ; mu++ ) { // hack for the moment
      char tmp[ STR_LENGTH ] ;
      sprintf( tmp , "DIM%d_%d" , file_counter , mu ) ;
      const int dim_idx = tag_search( tmp ) ;
      if( dim_idx == FAILURE ) { 
	return ;
      }
      INPARAMS -> dimensions[ file_counter ][ mu ] = atoi( INPUT[dim_idx].VALUE ) ;
    }
    printf( "DIMENSIONS_%d :: (" , file_counter ) ;
    for( mu = 0 ; mu < 4 ; mu++ ) { // hack for the moment
      printf( " %d " , INPARAMS -> dimensions[file_counter][ mu ] ) ;
    }
    printf( ")\n" ) ;
    file_counter++ ;
  }
  return ;
}

/// read the input file
int
read_inputs( struct input_params *INPARAMS ,
	     const char *filename )
{
  FILE *file = fopen( filename , "r" ) ;
  if( file == NULL ) {
    printf( "Cannot find input file %s \n" , filename ) ;
    return FAILURE ;
  }

  pack_inputs( file ) ;

  // look for the number of RAW data points
  // trajectory info
  int file_counter = 0 ;
  while( file_counter < MAX_NFILES ) {
    char str[ 128 ] ;
    sprintf( str , "TRAJ_BEG_%d" , file_counter ) ;
    const int trajbeg_idx = tag_search( str ) ;
    if( trajbeg_idx == FAILURE ) { break ; }
    INPARAMS -> traj_begin[file_counter] = atoi( INPUT[trajbeg_idx].VALUE ) ;

    sprintf( str , "TRAJ_END_%d" , file_counter ) ;
    const int trajend_idx = tag_search( str ) ;
    if( trajend_idx == FAILURE ) { break ; }
    INPARAMS -> traj_end[file_counter] = atoi( INPUT[trajend_idx].VALUE ) ;
    
    sprintf( str , "TRAJ_INC_%d" , file_counter ) ;
    const int trajinc_idx = tag_search( str ) ;
    if( trajinc_idx == FAILURE ) { break ; }
    INPARAMS -> traj_increment[file_counter] = atoi( INPUT[trajinc_idx].VALUE ) ;

    sprintf( str , "TRAJ_FILE_%d" , file_counter ) ;
    const int trajfile_idx = tag_search( str ) ;
    if( trajfile_idx == FAILURE ) { break ; }
    sprintf( INPARAMS -> traj_file[ file_counter ] , "%s" , INPUT[trajfile_idx].VALUE ) ;
    printf( "FILE :: %s \n" , INPARAMS -> traj_file[ file_counter ] ) ;

    sprintf( str , "BIN_%d" , file_counter ) ;
    const int bin_idx = tag_search( str ) ;
    INPARAMS -> binning[ file_counter ] = atoi( INPUT[bin_idx].VALUE ) ;

    // just to clear up some non-sensical possibilities
    if( ( INPARAMS -> traj_increment[file_counter] ) == 0 ) return FAILURE ;

    // compute NRAW from the trajectory info
    INPARAMS -> NRAW[ file_counter ] = 
      abs( ( INPARAMS -> traj_begin[ file_counter ] - 
	     INPARAMS -> traj_end[ file_counter ] ) 
	   / ( INPARAMS -> traj_increment[ file_counter ] ) ) ;
    printf( "Number of RAW points %d \n" , INPARAMS -> NRAW[ file_counter ] ) ;

    file_counter++ ;
  }

  // the number of files we have open
  INPARAMS -> NFILES = file_counter ;
  printf( "Analysing %d files \n" , INPARAMS -> NFILES ) ;

  // look for the data length, only used in the fake cases!
  INPARAMS -> NDATA[0] = atoi( set_token( "NDATA" ) ) ;

  // look at the reampling
  const int resample_idx = tag_search( "RESAMPLING" ) ;
  if( resample_idx == FAILURE ) { return tag_failure( "RESAMPLING" ) ; }
  if( are_equal( INPUT[resample_idx].VALUE , "BOOTSTRAP" ) ) {
    INPARAMS -> resample = BOOTDATA ;
    const int nboots_idx = tag_search( "NRESAMPLES" ) ;
    if( nboots_idx == FAILURE ) { return tag_failure( "NRESAMPLES" ) ; }
    INPARAMS -> NBOOTS = atoi( INPUT[nboots_idx].VALUE ) ;
  } else if( are_equal( INPUT[resample_idx].VALUE , "JACKKNIFE" ) ) {
    INPARAMS -> resample = JACKDATA ;
    INPARAMS -> NBOOTS = INPARAMS -> NRAW[0] ;
  } else {
    INPARAMS -> resample = RAWDATA ;
    INPARAMS -> NBOOTS = INPARAMS -> NRAW[0] ;
  }
  printf( "Attempting to %s resample, with %d resamples \n" ,
	  INPUT[resample_idx].VALUE , INPARAMS -> NBOOTS ) ;

  INPARAMS -> fit_hi = atof( set_token( "FIT_HI" ) ) ;
  INPARAMS -> fit_lo = atof( set_token( "FIT_LO" ) ) ;

  // sanity check fit ranges
#if 0
  if( INPARAMS -> fit_hi <= INPARAMS -> fit_lo ) {
    printf( "Cannot have upper fitrange %f \n" , INPARAMS -> fit_hi ) ;
    printf( "Below or equal to lower fitrange %f \n" , INPARAMS -> fit_lo ) ;
    printf( "That's just silly \n" ) ;
    return FAILURE ;
  } else if( INPARAMS -> fit_hi < 0 || INPARAMS -> fit_lo < 0 ) {
    printf( "Cannot have upper fitrange %f \n" , INPARAMS -> fit_hi ) ;
    printf( "Or lower fitrange %f less than 0\n" , INPARAMS -> fit_lo ) ;
    printf( "That's just silly \n" ) ;
    return FAILURE ;
  } else {
    printf( "Fitting over the range %f to %f \n" , 
	    INPARAMS -> fit_lo , INPARAMS -> fit_hi ) ;
  }
#endif
  printf( "Fitting over the range %f to %f \n" , 
	  INPARAMS -> fit_lo , INPARAMS -> fit_hi ) ;

  // set the chiral parameters
  set_chiral( INPARAMS ) ;

  // have a look at the fittype
  INPARAMS -> fittype = get_fittype( ) ;
  if( INPARAMS -> fittype == FAILURE ) return FAILURE ;

  // get the filetype default is fake data
  {
    const int analysis_idx = tag_search( "ANALYSIS" ) ;
    if( analysis_idx == FAILURE ) { return tag_failure( "ANALYSIS" ) ; }
    INPARAMS -> type = FAKE ;
    if( are_equal( INPUT[analysis_idx].VALUE , "ALPHA_S" ) ) {
      INPARAMS -> type = ALPHA_S ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "UKHADRON" ) ) {
      INPARAMS -> type = UKHADRON ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "KK_SS_EXTRAP" ) ) {
      INPARAMS -> type = KK_SS_EXTRAP ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "KK_ML_EXTRAP" ) ) {
      INPARAMS -> type = KK_ML_EXTRAP ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "TOPOLOGICAL_SUSCEPTIBILITY" ) ) {
      INPARAMS -> type = TOPOLOGICAL_SUSCEPTIBILITY ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "CORRELATORS" ) ) {
      INPARAMS -> type = CORRELATORS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "AA_CORRELATORS" ) ) {
      INPARAMS -> type = AA_CORRELATORS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "GEVP1_CORRELATORS" ) ) {
      INPARAMS -> type = GEVP1_CORRELATORS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "GEVP2_CORRELATORS" ) ) {
      INPARAMS -> type = GEVP2_CORRELATORS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "AMA_PP_CORRELATORS" ) ) {
      INPARAMS -> type = AMA_PP_CORRELATORS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "AMA_VV_CORRELATORS" ) ) {
      INPARAMS -> type = AMA_VV_CORRELATORS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "AMA_AA_CORRELATORS" ) ) {
      INPARAMS -> type = AMA_AA_CORRELATORS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "PP_CORRELATORS" ) ) {
      INPARAMS -> type = PP_CORRELATORS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "PP_BARYONS" ) ) {
      INPARAMS -> type = PP_BARYONS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "VV_CORRELATORS" ) ) {
      INPARAMS -> type = VV_CORRELATORS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "VV_BARYONS" ) ) {
      INPARAMS -> type = VV_BARYONS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "ZV_EVALUATE" ) ) {
      INPARAMS -> type = ZV_EVALUATE ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "TETRA_CORRELATORS" ) ) {
      INPARAMS -> type = TETRA_CORRELATORS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "FAKE" ) ) {
      INPARAMS -> type = FAKE ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "J0" ) ) {
      INPARAMS -> type = J0 ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "STATIC_POTENTIAL" ) ) {
      INPARAMS -> type = STATIC_POTENTIAL ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "TAUVUS" ) ) {
      INPARAMS -> type = TAUVUS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "TAUVUS2" ) ) {
      INPARAMS -> type = TAUVUS2 ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "TAUVUS3" ) ) {
      INPARAMS -> type = TAUVUS3 ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "DISPERSIONS" ) ) {
      INPARAMS -> type = DISPERSIONS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "CORRELATIONS" ) ) {
      INPARAMS -> type = CORRELATIONS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "MASS_SPLITTINGS" ) ) {
      INPARAMS -> type = MASS_SPLITTINGS ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "SPEED_OF_LIGHT" ) ) {
      INPARAMS -> type = SPEED_OF_LIGHT ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "SUPERDIST" ) ) {
      INPARAMS -> type = SUPERDIST ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "TMOMENTS" ) ) {
      INPARAMS -> type = TMOMENTS ;
    }  else if( are_equal( INPUT[analysis_idx].VALUE , "PADE_AMU" ) ) {
      INPARAMS -> type = PADE_AMU ;
    }  else if( are_equal( INPUT[analysis_idx].VALUE , "CONF_AMU" ) ) {
      INPARAMS -> type = CONF_AMU ;
    } else if(  are_equal( INPUT[analysis_idx].VALUE , "FLAVOUR_COMBINATION" ) ) {
      INPARAMS -> type = FLAVOUR_COMBINATION ;
    } else if( are_equal( INPUT[analysis_idx].VALUE , "FLOW_COUPLE" ) ) {
      INPARAMS -> type = FLOW_COUPLE ;
    }
  }

  // momentum average
  {
    const int momavg_idx = tag_search( "MOMAVG" ) ;
    if( momavg_idx == FAILURE ) { return tag_failure( "MOMAVG" ) ; }
    INPARAMS -> momavg = NONE ;
    if( are_equal( INPUT[momavg_idx].VALUE , "ZNM1_AVERAGE" ) ) {
      INPARAMS -> momavg = Z_Nm1_AVERAGE ;
    } else if( are_equal( INPUT[momavg_idx].VALUE , "PSQ_AVERAGE" ) ) {
      INPARAMS -> momavg = PSQ_AVERAGE ;
    }
  }

  // momentum definition
  {
    const int momtype_idx = tag_search( "MOMTYPE" ) ;
    if( momtype_idx == FAILURE ) { return tag_failure( "MOMTYPE" ) ; }
    INPARAMS -> mom_type = TWOSIN_MOM ;
    if( are_equal( INPUT[momtype_idx].VALUE , "TWOSIN_MOM" ) ) {
    printf( "[MOM] 2 Sin( p/2 ) lattice momentum definition\n" ) ;
      INPARAMS -> mom_type = TWOSIN_MOM ;
    } else if( are_equal( INPUT[momtype_idx].VALUE , "PSQ_MOM" ) ) {
     printf( "[MOM] Fourier mode lattice momentum definition\n" ) ;
      INPARAMS -> mom_type = PSQ_MOM ;
    } else if( are_equal( INPUT[momtype_idx].VALUE , "SIN_MOM" ) ) {
      printf( "[MOM] sin(p) lattice momentum definition\n" ) ;
      INPARAMS -> mom_type = SIN_MOM ;
    } else if( are_equal( INPUT[momtype_idx].VALUE , "RSQ" ) ) {
      printf( "[MOM] r^2 lattice momentum definition\n" ) ;
      INPARAMS -> mom_type = RSQ ;
    } else {
      printf( "I don't understand your momentum type %s\n" ,
	      INPUT[momtype_idx].VALUE ) ;
      return FAILURE ;
    }
  }

  // time fold
  {
    const int tfold_idx = tag_search( "TFOLD" ) ;
    if( tfold_idx == FAILURE ) { return tag_failure( "TFOLD" ) ; }
    INPARAMS -> tfold = false ;
    if( are_equal( INPUT[tfold_idx].VALUE , "TRUE" ) ) {
      INPARAMS -> tfold = true ;
    }
  }

  // check what simultaneous fit params we have
  set_simultaneous( INPARAMS ) ;

  // set the dimensions
  set_dimensions( INPARAMS ) ;

  // seed number :: 0 generates a Seed from /dev/urandom !
  sscanf( set_token( "RNG_SEED" ) , "%lu" , &(INPARAMS -> Seed ) ) ;

  // graph name
  sprintf( INPARAMS -> graph_name , "%s" ,  set_token( "GRAPHNAME" ) ) ;

  // graph name
  sprintf( INPARAMS -> graph_xaxis, "%s" ,  set_token( "XAXIS" ) ) ;

  // graph name
  sprintf( INPARAMS -> graph_yaxis , "%s" ,  set_token( "YAXIS" ) ) ;

  // graph name
  sprintf( INPARAMS -> output_file , "%s" ,  set_token( "OUTFILE" ) ) ;

  // close and unpack
  unpack_inputs( ) ;
  fclose( file ) ;

  return SUCCESS ;
}
	     
