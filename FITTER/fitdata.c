#include "fitfunc.h"
#include "fit_chooser.h"
#include "Utils.h"

#include "chisq_check.h"

#define USE_AVERAGE

//#define VERBOSE

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

#ifdef VERBOSE
// status functions
static void
print_state( const size_t iter , 
	     gsl_multifit_fdfsolver *s ,
	     const int NPARAMS )
{
  int i ;
  printf ("\niter: %lu \n" , iter ) ;
  for( i = 0 ; i < NPARAMS ; i++ ) {
    printf( "param[%d] %15.8f \n" , i , FIT(i) ) ;
  } 
  printf( "|f(x)| = %e\n\n" , gsl_blas_dnrm2( s->f ) ) ;
  return ;
}
#endif

// bit that performs the fit
static int
perform_fit( gsl_vector_view *__restrict x , 
	     gsl_multifit_function_fdf f ,
	     double *__restrict params ,
	     double *__restrict chisq ,
	     struct data d )
{
  const int p = d.SIMS * ( d.NPARAMS - d.NCOMMON ) + d.NCOMMON ;

  // compute the range
  int j , xposit = 0 , range = 0 ;
  for( j = 0 ; j < d.SIMS ; j++ ) {
    const int fitlo_idx = find_idx( d.FIT_LO , d.X , xposit + d.NDATA[j] , xposit ) ;
    const int fithi_idx = find_idx( d.FIT_HI , d.X , xposit + d.NDATA[j] , xposit ) ;
    range += ( 1 + fithi_idx - fitlo_idx ) ;
    xposit += d.NDATA[j] ;
  }

  // covariance matrix
  gsl_matrix *covar = gsl_matrix_alloc ( p , p ) ;

  ///////////////////////////////////////////////////////

  const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
  //const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmder;
  gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc ( T , range , p ) ;

  //gsl_multifit_function_fdf f ;
  int status = 0 ;

  // and set up the data ...
  f.n = range ;
  f.p = p ;
  f.params = &d ;

  gsl_multifit_fdfsolver_set( s , &f , &(x->vector) ) ;

  int iter = 0 ;
  do
    {
      iter++;

      // sometimes the solver complains too much
      status = gsl_multifit_fdfsolver_iterate( s ) ;

      #ifdef VERBOSE
        print_state( iter , s , p ) ;
	printf( "status = %s\n" , gsl_strerror( status ) ) ;
      #endif

      status = gsl_multifit_test_delta ( s->dx , s->x ,
					 1E-7 , 1E-7 ) ;
      
      #ifdef VERBOSE
        printf( "status = %s\n" , gsl_strerror( status ) ) ;
      #endif
    }
  while( status != GSL_SUCCESS && iter < 10000 ) ;

  if( iter == 10000 ) {
    printf( "fit is failing to converge ... Leaving \n" ) ;
    printf( "status = %s\n" , gsl_strerror( status ) ) ;
    return GSL_FAILURE ;
  }

  gsl_multifit_covar( s->J , 0.0 , covar ) ;

  *chisq = pow( gsl_blas_dnrm2( s->f ) , 2 ) / ( range - p ) ;
#ifndef QUIET
  printf( "CHISQ/dof %e :: iter %d \n\n" , *chisq , iter ) ;
#endif

  for( j = 0 ; j < p ; j++ ) { params[ j ] = FIT( j ) ; }

  gsl_multifit_fdfsolver_free( s ) ;
  gsl_matrix_free( covar ) ;

  return GSL_SUCCESS ;
}

static struct data
setdata( const double *X ,
	 const double *sigma ,
	 const struct mom_info *mominfo ,
	 const struct chiral *quarks ,
	 const double FIT_HI ,
	 const double FIT_LO ,
	 const int *NPARAMS ,
	 const int SIMS , 
	 const int NCOMMON ,
	 const int NFLAT ,
	 const int NPARS ,
	 const int LT ,
	 const int *NDATA ,
	 const bool *sim_params )
{
  // pack the local data struct
  struct data d ;
  // allocations
  d.X = (const double*)X ;
  d.sigma = (const double*)sigma ;
  d.mom = (struct mom_info*)mominfo ;
  d.SIMS = SIMS ;
  d.NCOMMON = NCOMMON ;
  d.n = NFLAT ;
  d.NDATA = NDATA ;
  d.NPARAMS = NPARAMS[0] ;
  d.quarks = (const struct chiral*)quarks ;
  d.data_idx = 0 ;
  // and the fit range ...
  d.FIT_HI = FIT_HI ; 
  d.FIT_LO = FIT_LO ;
  // fit params
  d.LOGICAL_NPARS = NPARS ;
  d.LT = LT ;
  // set up the simultaneous fit shared param(s)?
  int j ;
  for( j = 0 ; j < d.NPARAMS ; j++ ) {
    d.sim_params[ j ] = sim_params[ j ] ;
  }
  return d ;
}

// fits the data, returns an allocated fitparams array
struct resampled *
perform_bootfit( fitfunc *fit ,
		 int *NPARAMS ,
		 const struct resampled *boots , // the bootstrapped data 
		 const double *X , // the (sorted) x-axis data
		 const double *sigma , // the s.d of the y-data
		 const struct mom_info *mominfo ,
		 const struct chiral *quarks , // chiral data
		 const fit_type fittype ,
		 const int *NDATA , // the length of the X-array
		 const double FIT_HI ,
		 const double FIT_LO ,
		 const int SIMS , // number of simultaneous datasets
		 const int LT ,
		 const bool sim_params[ 12 ] ,
		 const int NCOMMON )
{
  // compute number of flat parameters ...
  int NFLAT = 0 , i ;
  for( i = 0 ; i < SIMS ; i++ ) {
    NFLAT += NDATA[i] ;
  }

  // initialise the fit
  int replicas , NCO = NCOMMON ;

  // set up the fitfunc
  const gsl_multifit_function_fdf f = initialise_f( fittype , fit , &NPARAMS[0] , &NCO ) ;

  const int NPARS = SIMS * ( NPARAMS[0] - NCOMMON ) + NCOMMON ;
    
  // allocate the array
  struct resampled *fitparams = malloc( NPARS * sizeof( struct resampled ) ) ;
  int j ;
  for( j = 0 ; j < NPARS ; j++ ) {
    fitparams[j].resampled = malloc( ( boots[0].NSAMPLES ) * sizeof( double ) ) ;
    fitparams[j].NSAMPLES  = boots[0].NSAMPLES ;
    fitparams[j].restype   = boots[0].restype ;
  }

  const int NBOOTS = boots[0].NSAMPLES ;

  // allocate chisq ...
  struct resampled chisq ;
  chisq.resampled = malloc( NBOOTS * sizeof( double ) ) ;
  chisq.NSAMPLES  = NBOOTS ;
  chisq.restype   = boots[0].restype ;

  double ave_params[ NPARS ] ;
  // compute the value for the average
  {
    // pack the local data struct
    struct data d = setdata( X , sigma , mominfo , quarks ,
			     FIT_HI , FIT_LO , NPARAMS , SIMS , 
			     NCOMMON , NFLAT , NPARS , LT ,
			     NDATA , sim_params ) ;
      
    // set the local y to the bootstrapped one ...
    double *ytemp = malloc( NFLAT * sizeof( double ) ) ;
    for( j = 0 ; j < NFLAT ; j++ ) {
      ytemp[ j ] = boots[ j ].avg ;
    }
    d.y = (const double*)ytemp ;
      
    // fit params can go here ... Initialise to 0
    fit -> guess( ave_params , &d , NPARAMS[0] ) ;

    // provide some guesses, use the results of previous fits
    gsl_vector_view x = gsl_vector_view_array( ave_params , NPARS ) ;

    // y = x ;
    if( perform_fit( &x , f , ave_params , &chisq.avg , d ) 
	== GSL_FAILURE ) {
      exit( 1 ) ;
    }

    // set params
    for( j = 0 ; j < NPARS ; j++ ) {
      fitparams[ j ].avg = ave_params[ j ] ;
    }
      
    // end of the boots loop
    free( ytemp ) ;
  }

  // loop over the number of boots performing the fit
#pragma omp parallel for private(replicas)
  for( replicas = 0 ; replicas < NBOOTS ; replicas++ ) {

    // set the data struct ...
    struct data d  = setdata( X , sigma , mominfo , quarks ,
			      FIT_HI , FIT_LO , NPARAMS , SIMS , 
			      NCOMMON , NFLAT , NPARS , LT ,
			      NDATA , sim_params ) ;
      
    // set the local y to the bootstrapped one ...
    double *ytemp = malloc( NFLAT * sizeof( double ) ) ;
    for( j = 0 ; j < NFLAT ; j++ ) {
      ytemp[ j ] = boots[ j ].resampled[ replicas ] ;
    }
    d.y = (const double*)ytemp ;
      
    // fit params can go here ... Initialise to 0
    double params[ NPARS ] ;
#if defined USE_AVERAGE
    for( j = 0 ; j < NPARS ; j++ ) {
      params[ j ] = ave_params[ j ] ;
    }
#else
    fit -> guess( params , &d , NPARS ) ;
#endif

    // provide some guesses, use the results of previous fits
    gsl_vector_view x = gsl_vector_view_array( params , NPARS ) ;
      
    // y = x ;
    if( perform_fit( &x , f , params , &chisq.resampled[ replicas ] , d ) 
	== GSL_FAILURE ) {
      exit( 1 ) ;
    }

    // set params
    for( j = 0 ; j < NPARS ; j++ ) {
      fitparams[ j ].resampled[ replicas ] = params[ j ] ;
    }
      
    // end of the boots loop
    free( ytemp ) ;
  }

  // tell us what the distribution of error is 
  compute_err( &chisq ) ;
  printf( "FULL_CHISQ :: %f +/- %f \n" , chisq.avg , chisq.err ) ;

  // tell us what the average of the fit params is
  for( replicas = 0 ; replicas < NPARS ; replicas++ ) {
    compute_err( &fitparams[ replicas ] ) ;
    if( sim_params[replicas] == true ) {
      printf( "[SIMULTANEOUS] " ) ;
    }
    printf( "FIT PARAM %d :: %1.10f +/- %1.10f \n" , 
	    replicas ,
	    fitparams[replicas].avg , 
	    0.5 * ( fitparams[replicas].err_hi - 
		    fitparams[replicas].err_lo ) ) ;
  }

  free( chisq.resampled ) ;

  return fitparams ;
}
