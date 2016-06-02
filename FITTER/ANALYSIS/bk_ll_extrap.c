/**
   @file bk_ss_extrap.c
   @brief linear extrapolation in the strange, write out the result
 */
#include "fitfunc.h"

#include "fit_and_plot.h"
#include "GLU_bswap.h"
#include "svd.h"
#include "write_distribution.h"
#include "graph_data.h"
#include "fit_chooser.h"

#ifdef EPSILON_FIT
static const double extrap = 0.01821868254756 / ( 0.13041 * 0.13041 ) ;
//static const double extrap = 0.01946025 / ( 0.13041 * 0.13041 ) ;
#else
static const double extrap = 0.01821868254756 / ( 0.12229 * 0.12229 ) ;
//static const double extrap = 0.01946025 / ( 0.12229 * 0.12229 ) ;
#endif

static void
write_data( struct resampled data , 
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

static void
plot_extrap( struct resampled *fit1 , 
	     double **xavg ,
	     struct resampled **bootavg ,
	     struct mom_info **mominfo ,
	     struct mom_info *moms ,
	     struct input_params *INPARAMS ,
	     const int NSLICES ,
	     const int LT ,
	     const char *type , 
	     const int N )
{
  struct input_params temp = *INPARAMS ;
  fitfunc fit ;
  int NPARAMS , NCOMMON ;
  struct resampled y_cont , y_fine ;
  initialise_f( INPARAMS->fittype , &fit , &NPARAMS , &NCOMMON ) ;
  plot_fitfunc2( fit1 , bootavg[0] , xavg[0] , mominfo[0] , fit , 3 ,
		 INPARAMS->quarks[0] , INPARAMS -> fit_lo , INPARAMS -> fit_hi ,
		 LT , INPARAMS -> NDATA[1] , 0 ) ;
  plot_fitfunc2( fit1 , bootavg[1] , xavg[1] , mominfo[1] , fit , 3 ,
		 INPARAMS->quarks[3] , INPARAMS -> fit_lo , INPARAMS -> fit_hi ,
		 LT , INPARAMS -> NDATA[1] , 1 ) ;
  y_fine = extrapolate( fit1 , extrap , 3 , temp.quarks[0] , mominfo[0][0] , fit , LT , 3 ) ;

  // plot the continuum extrapolation
  temp.quarks[0].ainverse = 1E9 ;
  plot_fitfunc2( fit1 , bootavg[0] , xavg[0] , mominfo[0] , fit , 3 ,
		 temp.quarks[0] , INPARAMS -> fit_lo , INPARAMS -> fit_hi ,
		 LT , INPARAMS -> NDATA[0] , 1 ) ;
  y_cont = extrapolate( fit1 , extrap , 3 , temp.quarks[0] , mominfo[0][0] , fit , LT , 3 ) ;
  printf( "[%s] CHIRAL-CONTINUUM %d :: %1.12f +/- %f \n" , type , N , y_cont.avg , y_cont.err ) ;

  subtract( &y_cont , y_fine ) ;
  mult_constant( &y_cont , 0.5 ) ;

  printf( "[%s] DISCRETISATION %d :: %1.12f +/- %f \n" , type , N , y_cont.avg , y_cont.err ) ;

  free( y_fine.resampled ) ;
  free( y_cont.resampled ) ;
  return ;
}

void
ll_extrap_eval( double **xavg ,
		struct resampled **bootavg ,
		struct mom_info **mominfo ,
		struct mom_info *moms ,
		struct input_params *INPARAMS ,
		const int NSLICES ,
		const int LT )
{
  // do the fit
  const struct resampled *fitparams = fit_data_plot_data( (const struct resampled**)bootavg , 
							  (const double**)xavg , 
							  (const struct mom_info **)mominfo ,
							  *INPARAMS , NSLICES , LT ) ;
  graph_reset_color( ) ;

  // initialise the xmgrace graph
  make_xmgrace_graph( INPARAMS -> graph_name ,
		      INPARAMS -> graph_xaxis , 
		      INPARAMS -> graph_yaxis ) ;
    
  // plot the data first but keep the graph file open
  size_t i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    plot_data( bootavg[i] , xavg[i] , INPARAMS -> NDATA[i] ) ;
  }
  printf( "Plotted\n" ) ;
  graph_reset_color( ) ;

  if( INPARAMS -> fittype == CHIPT_SIM_R45_RENORM ||
      INPARAMS -> fittype == CHIPT_SIM_R23_RENORM ||
      INPARAMS -> fittype == CHIPT_SIM_B45_RENORM ||
      INPARAMS -> fittype == CHIPT_SIM_B23_RENORM ||
      INPARAMS -> fittype == CHIPT_SIM_R23_SUSY ||
      INPARAMS -> fittype == CHIPT_SIM_R45_SUSY ||
      INPARAMS -> fittype == CHIPT_SIM_B23_SUSY ||
      INPARAMS -> fittype == CHIPT_SIM_B45_SUSY ) {
    
    // plot the fitfuncs
    struct resampled *fit1 = malloc( 3 * sizeof( struct resampled ) ) ;
    for( i = 0 ; i < 3 ; i++ ) {
      fit1[i].resampled = (double*)malloc( fitparams[0].NSAMPLES * sizeof( double ) ) ;
    }

    equate( &fit1[0] , fitparams[0] ) ;
    equate( &fit1[1] , fitparams[1] ) ;
    equate( &fit1[2] , fitparams[2] ) ;
    plot_extrap( fit1 , xavg , bootavg , mominfo , moms , INPARAMS , 
		 NSLICES , LT , "CHIPT3" , 1 ) ;


    equate( &fit1[0] , fitparams[3] ) ;
    equate( &fit1[1] , fitparams[1] ) ;
    equate( &fit1[2] , fitparams[4] ) ;
    plot_extrap( fit1 , xavg , bootavg , mominfo , moms , INPARAMS , 
		 NSLICES , LT , "CHIPT3" , 2 ) ;

    equate( &fit1[0] , fitparams[5] ) ;
    equate( &fit1[1] , fitparams[1] ) ;
    equate( &fit1[2] , fitparams[6] ) ;
    plot_extrap( fit1 , xavg , bootavg , mominfo , moms , INPARAMS , 
		 NSLICES , LT , "CHIPT3" , 2 ) ;

    equate( &fit1[0] , fitparams[7] ) ;
    equate( &fit1[1] , fitparams[1] ) ;
    equate( &fit1[2] , fitparams[8] ) ;
    plot_extrap( fit1 , xavg , bootavg , mominfo , moms , INPARAMS , 
		 NSLICES , LT , "CHIPT3" , 2 ) ;
    
    // close up the graph
    close_xmgrace_graph(  ) ;
    
    printf( "\n---> Graph plotted to %s <---\n" , INPARAMS -> graph_name) ;
  }

  // simfit over two parameters
  if( INPARAMS -> fittype == CHIPT2_SIM_R45_RENORM ||
      INPARAMS -> fittype == CHIPT2_SIM_R23_RENORM ||
      INPARAMS -> fittype == CHIPT2_SIM_B45_RENORM ||
      INPARAMS -> fittype == CHIPT2_SIM_B23_RENORM ||
      INPARAMS -> fittype == CHIPT2_SIM_R23_SUSY ||
      INPARAMS -> fittype == CHIPT2_SIM_R45_SUSY ||
      INPARAMS -> fittype == CHIPT2_SIM_B23_SUSY ||
      INPARAMS -> fittype == CHIPT2_SIM_B45_SUSY ) {
    
    // plot the fitfuncs
    struct resampled *fit1 = malloc( 3 * sizeof( struct resampled ) ) ;
    for( i = 0 ; i < 3 ; i++ ) {
      fit1[i].resampled = (double*)malloc( fitparams[0].NSAMPLES * sizeof( double ) ) ;
    }

    equate( &fit1[0] , fitparams[0] ) ;
    equate( &fit1[1] , fitparams[1] ) ;
    equate( &fit1[2] , fitparams[2] ) ;
    plot_extrap( fit1 , xavg , bootavg , mominfo , moms , INPARAMS , 
		 NSLICES , LT , "CHIPT2" , 1 ) ;


    equate( &fit1[0] , fitparams[3] ) ;
    equate( &fit1[1] , fitparams[1] ) ;
    equate( &fit1[2] , fitparams[4] ) ;
    plot_extrap( fit1 , xavg , bootavg , mominfo , moms , INPARAMS , 
		 NSLICES , LT , "CHIPT2" , 2 ) ;
    
    // close up the graph
    close_xmgrace_graph(  ) ;
    
    printf( "\n---> Graph plotted to %s <---\n" , INPARAMS -> graph_name) ;
  }

  // non simfit
  if( INPARAMS -> fittype == CHIPT_R145_RENORM ||
      INPARAMS -> fittype == CHIPT_R23_RENORM ||
      INPARAMS -> fittype == CHIPT_B145_RENORM ||
      INPARAMS -> fittype == CHIPT_B23_RENORM ||
      INPARAMS -> fittype == CHIPT_R123_SUSY ||
      INPARAMS -> fittype == CHIPT_R45_SUSY ||
      INPARAMS -> fittype == CHIPT_B123_SUSY ||
      INPARAMS -> fittype == CHIPT_B45_SUSY || 
      INPARAMS -> fittype == BKCHIRAL ) {

    plot_extrap( (struct resampled*)fitparams , xavg , bootavg , mominfo , moms , INPARAMS , 
		 NSLICES , LT , "CHIPT1" , 0 ) ;
    
    //
    struct mom_info moms[1] ;
    moms[0].n[0] = moms[0].n[1] = moms[0].n[2] = 0 ;
    write_data( fitparams[0] , moms , INPARAMS -> output_file ) ;
  }

  free( (void*)fitparams ) ;

  return ;
}
