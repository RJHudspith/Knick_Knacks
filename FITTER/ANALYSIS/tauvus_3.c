#include "fitfunc.h"

#include "graph_data.h"
#include "Utils.h"

// little enum for the type of knobs we like
typedef enum { CENTRE_POSITION , DELTA_SIZE , NWEIGHTS } knobtype ;

// tau mass squared in GeV^2, gets converted to lattice units at 
// the start of the code
static double mtsquared = 1.77682 * 1.77682 ;

// struct for storing our weight parameters
struct weight_params {
  size_t Nweights ;
  double centre ;
  double delta ;
  double MAX ;
  double MIN ;
  double STEP ;
} ;

// gets the appropriate weights
static void
get_weights( const size_t Nweights ,
	     double weights[ Nweights ] ,
	     const double centre , 
	     const double delta )
{
  int j , idx = 0 ;
  // N is even the mid two weights are either side of it
  if( (Nweights&1) == 0 ) {
    for( j = -(int)Nweights/2 ; j < (int)Nweights/2 ; j++ ) {
      weights[idx] = centre + j*delta + delta/2.0 ;
      idx++ ;
    } 
  // N is odd the centre is one of the weights
  } else {
    for( j = -(int)Nweights/2 ; j < (int)Nweights/2+1 ; j++ ) {
      weights[idx] = centre + j*delta ;
      idx++ ;
    }
  }
  return ;
}

// DFT for positive values of Q
static struct resampled 
arbq2( const struct resampled *corr , 
       const double Q2 ,
       const momtype mom_type ,
       const int L )
{
  // allocate result
  struct resampled res ;
  res.resampled = malloc( corr[0].NSAMPLES * sizeof( double ) ) ;
  equate_constant( &res , 0.0 , corr[0].NSAMPLES , corr[0].restype ) ;
  // compute twiddle factors
  double mom = 0.0 ;
  switch( mom_type ) {
  case TWOSIN_MOM : mom = 2.0 * asin( 0.5 * sqrt( Q2 ) ) ; break ;
  case PSQ_MOM    : mom = sqrt( Q2 )                     ; break ;
  default : printf( "UNRECOGNISED MOMTYPE ... Leaving \n" ) ; exit(1) ; break ;
  }
  // cache some important values for later
  double coscache[ L ] , tsqcache[ L ] ; 
  const int L_2 = L>>1 ;
  int t , tp ;
  for( t = 0 ; t < L ; t++ ) {
    tp = t < L_2 ? t : t - L ;
    coscache[ t ] = cos( mom * tp ) ;
    tsqcache[ t ] = (double)tp*tp   ;
  }
  // perform the DFT subtracting the two possible zeros
  size_t j ;
  for( j = 0 ; j < res.NSAMPLES ; j++ ) {
    register double sum1 = 0.0 , sum2 = 0.0 , sum3 = 0.0 , ct ;
    for( t = 0 ; t < L ; t++ ) {
      ct = corr[t].resampled[j] ;
      sum1 += ct * coscache[ t ] ;
      sum2 += ct ;
      sum3 += ct * tsqcache[t] ;
    }
    res.resampled[j] = ( Q2 != 0.0 ) ? ( sum1 - sum2 ) / ( Q2 ) + 0.5 * sum3 : 0.0 ;
  }
  compute_err( &res ) ;
  return res ;
} 

// computes
// \sum_{k} \Pi( Q(k) ) / ( \prod_{j\neq k} (Q(x^2) - Q(j)^2 ) )
// weights are (aq)^2 - valued
struct resampled
weighted_sum( const size_t Nweights , 
	      const double weights[ Nweights ] ,
	      const struct resampled *VpA0 ,
	      const struct resampled *VpA1 ,
	      const momtype mom_type ,
	      const size_t LT )
{
  struct resampled sum , temp0 , temp1 ;
  sum.resampled   = malloc( VpA0[0].NSAMPLES * sizeof( double ) ) ;
  temp0.resampled = malloc( VpA0[0].NSAMPLES * sizeof( double ) ) ;
  temp1.resampled = malloc( VpA0[0].NSAMPLES * sizeof( double ) ) ;
  equate_constant( &sum , 0.0 , VpA0[0].NSAMPLES , VpA0[0].restype ) ;
  size_t k , j ;
  for( k = 0 ; k < Nweights ; k++ ) {
    // compute the product
    register double prod = 1.0 ;
    for( j = 0 ; j < Nweights ; j++ ) {
      if( j != k ) {
	prod *= ( weights[j] - weights[k] ) ;
      }
    }
    // is the longitudinal component
    equate( &temp0 , arbq2( VpA0 , weights[k] , mom_type , LT ) ) ;
    // is the transverse gets multiplied by (1-2*(q/m_tau)^2)
    equate( &temp1 , arbq2( VpA1 , weights[k] , mom_type , LT ) ) ;
    mult_constant( &temp1 , ( 1.0 - 2.0 * weights[k] / mtsquared ) ) ;
    // add the two contributions together and divide by the weights product
    add( &temp0 , temp1 ) ;
    mult_constant( &temp0 , 1.0 / prod ) ;
    add( &sum , temp0 ) ;
  }
  free( temp0.resampled ) ;
  free( temp1.resampled ) ;
  return sum ;
}

// plot our momentum-space correlator
static void
plot_momentum_corr( const struct resampled *corr ,
		    const double *x ,
		    const size_t NDATA ,
		    const momtype mom_type ,
		    const int LT )
{
  struct resampled *ydata = malloc( NDATA * sizeof( struct resampled ) ) ;
  size_t j ;
#pragma omp parallel for private(j)
  for( j = 0 ; j < NDATA ; j++ ) {
    ydata[j] = arbq2( corr , x[j] , mom_type , LT ) ;
  }
  plot_data( ydata , x , NDATA ) ;
  free_resampled( ydata , NDATA ) ;
  return ;
}

// looks at the possible knobs we can twiddle
static void
loop_possibles( const struct resampled *VpA0 ,
		const struct resampled *VpA1 ,
		const struct weight_params w ,
		const double a2 ,
		const momtype mom_type , 
		const int LT ,
		const knobtype KNOB )
{
  const size_t nevals = (size_t)(( w.MAX - w.MIN )/ w.STEP ) ;
  switch( KNOB ) {
  case CENTRE_POSITION :
    make_xmgrace_graph( "Central.agr" , "Centre [GeV\\S2\\N]" , "Percentage error" ) ;
    break ;
  case DELTA_SIZE :
    make_xmgrace_graph( "Delta.agr" , "Delta [GeV\\S2\\N]" , "Percentage error" ) ;
    break ;
  case NWEIGHTS : 
    make_xmgrace_graph( "Weights.agr" , "Weight #" , "Percentage error" ) ;
    break ;
  }
  double *xdata = malloc( nevals * sizeof( double ) ) , *ydata = malloc( nevals * sizeof( double ) ) ;
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < nevals ; i++ ) {
    double weights[ w.Nweights ] , delta = w.delta , centre = w.centre ;
    int Nweights = w.Nweights ;
    // consider our knobs
    switch( KNOB ) {
    case CENTRE_POSITION : centre = (w.MIN + i*w.STEP) ; break ;
    case DELTA_SIZE :      delta = (w.MIN + i*w.STEP)  ; break ;
    case NWEIGHTS : Nweights = (int)(w.MIN/a2)+i       ; break ;
    }
    // get the various weights
    get_weights( Nweights , weights , centre , delta ) ;
    struct resampled res = weighted_sum( Nweights , weights , VpA0 , VpA1 , mom_type , LT ) ;
    // fractional error 
    ydata[i] = fabs( res.err / res.avg ) * 100 ;
    xdata[i] = (w.MIN + i*w.STEP)/a2 ; // convert to physical?
    free( res.resampled ) ;
  }
  plot_simple_data( ydata , xdata , nevals , 
		    "(1-2q\\S2\\N/m\\s\\xt\\f{}\\N\\h{-0.3}\\S2\\N)\\xP\\f{}\\S1:V+A\\N"
		    "+\\xP\\f{}\\S0:V+A\\N" ) ;
  free( xdata ) ; free( ydata ) ;
  close_xmgrace_graph( ) ;
  return ;
}

// set weight parameters to lattice values
// needs to happen because our arbitrary momentum evaluation
// works with lattice-valued quantities but we want the weights 
// in physical pole positions
static void
lattice_weight_parameters( struct weight_params *w , 
			   const double a2 )
{
  w -> centre *= a2 ;
  w -> delta  *= a2 ;
  w -> MAX    *= a2 ;
  w -> MIN    *= a2 ;
  w -> STEP   *= a2 ;
}

// 
void
tauvus3_eval( double **xavg ,
	      struct resampled **bootavg ,
	      struct mom_info **mominfo ,
	      struct input_params *INPARAMS ,
	      const int NSLICES ,
	      const int LT ,
	      const bool renormalise )
{
  printf( "\n--> Tau vus (3) evaluation <--\n" ) ;

  // lattice spacing squared
  const double a2 = 1.0 / pow( INPARAMS -> quarks[0].ainverse , 2 ) ;

  // set mtsquared to lattice units (am)^2
  mtsquared *= a2 ;

  // renormalise tcorrs, i.e. multiply everything by ZV
  printf( "\n--> Renormalising <--\n" ) ;
  size_t i , j ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    for( j = 0 ; j < INPARAMS->NDATA[i] ; j++ ) {
      #ifdef verbose
      printf( "[TMOMENTS] corr :: %zu %zu %f %e %e \n" , i , j , xavg[ i ][j] , 
	      bootavg[i][j].avg , bootavg[i][j].err ) ;
      #endif
      mult_constant( &bootavg[i][j] , INPARAMS->quarks[0].ZV ) ;
    }
  }

  printf( "\n--> Computing V+A Pis <--\n" ) ;
  struct resampled *VpA1 = malloc( INPARAMS -> NDATA[0] * sizeof( struct resampled ) ) ;
  struct resampled *VpA0 = malloc( INPARAMS -> NDATA[0] * sizeof( struct resampled ) ) ;
  for( j = 0 ; j < INPARAMS -> NDATA[0] ; j++ ) {
    VpA0[j].resampled = malloc( bootavg[ 0 ][ j ].NSAMPLES * sizeof( double ) ) ;
    VpA1[j].resampled = malloc( bootavg[ 0 ][ j ].NSAMPLES * sizeof( double ) ) ;
    // longitudinal
    equate( &VpA0[j] , bootavg[ 2 ][ j ] ) ;
    add( &VpA0[j] , bootavg[ 3 ][ j ] ) ;
    // transverse
    equate( &VpA1[j] , bootavg[ 0 ][ j ] ) ;
    add( &VpA1[j] , bootavg[ 1 ][ j ] ) ;
  }

  printf( "\n--> Plotting the VPFs <--\n" ) ;
  {
    make_xmgrace_graph( "HVP.agr" , "(aq)\\S2\\N" , "\\xP\\f{}" ) ;
    plot_momentum_corr( bootavg[0] , xavg[4] , INPARAMS->NDATA[4] , INPARAMS -> mom_type , LT ) ;
    plot_momentum_corr( bootavg[1] , xavg[5] , INPARAMS->NDATA[5] , INPARAMS -> mom_type , LT ) ;
    plot_momentum_corr( bootavg[2] , xavg[6] , INPARAMS->NDATA[6] , INPARAMS -> mom_type , LT ) ;
    plot_momentum_corr( bootavg[3] , xavg[7] , INPARAMS->NDATA[7] , INPARAMS -> mom_type , LT ) ;
    //plot_momentum_corr( VpA0 , xavg[6] , INPARAMS->NDATA[6] , INPARAMS -> mom_type , LT ) ;
    //plot_momentum_corr( VpA1 , xavg[7] , INPARAMS->NDATA[7] , INPARAMS -> mom_type , LT ) ;
    close_xmgrace_graph() ;
  }

  printf( "\n--> Measuring the sum <--\n" ) ;
  // central dependence
  {
    struct weight_params w = { .Nweights = 3 , .centre = 1.0 , .delta = 0.1 ,
			       .MAX = 4.0 , .MIN = 0.2 , .STEP = 0.05 } ;
    lattice_weight_parameters( &w , a2 ) ;
    loop_possibles( VpA0 , VpA1 , w , a2 , INPARAMS -> mom_type , LT , CENTRE_POSITION ) ;
  }
  // delta dependence
  {
    struct weight_params w = { .Nweights = 3 , .centre = 1.0 , .delta = 0.1 ,
			       .MAX = 1.0 , .MIN = 0.05 , .STEP = 0.01 } ;
    lattice_weight_parameters( &w , a2 ) ;
    loop_possibles( VpA0 , VpA1 , w , a2 , INPARAMS -> mom_type , LT , DELTA_SIZE ) ;
  }
  // Nweights dependence
  {
    struct weight_params w = { .Nweights = 8 , .centre = 1.0 , .delta = 0.05 ,
			       .MAX = 10 , .MIN = 2 , .STEP = 1 } ;
    lattice_weight_parameters( &w , a2 ) ;
    loop_possibles( VpA0 , VpA1 , w , a2 , INPARAMS -> mom_type , LT , NWEIGHTS ) ;
  }
  
  // free the V+A longitudinal and transverse components
  free_resampled( VpA0 , INPARAMS -> NDATA[0] ) ;
  free_resampled( VpA1 , INPARAMS -> NDATA[0] ) ;

  return ;
}
