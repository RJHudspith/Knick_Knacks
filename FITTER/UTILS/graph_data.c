#include "fitfunc.h"
#include "stats.h"
#include "Utils.h"
#include "D0_diff.h"
#include "D0_diff_multi.h"

static FILE *file ;
static int dataset = 0 , colorset = 1 ;

static void
PrintHeader( void )
{
  const int symbol = ( dataset + 1 ) % 11 ;
  fprintf( file , "@\ts%d hidden false\n" , dataset ) ;  
  fprintf( file , "@\ts%d symbol %d\n" , dataset , symbol ) ;  
  fprintf( file , "@\ts%d symbol size 0.33\n" , dataset ) ;  
  fprintf( file , "@\ts%d symbol color %d\n" , dataset , colorset ) ;  
  fprintf( file , "@\ts%d symbol fill color %d\n" , dataset , colorset ) ;  
  fprintf( file , "@\ts%d symbol fill pattern 1\n" , dataset ) ;  
  fprintf( file , "@\ts%d line type 0\n" , dataset ) ;  
  fprintf( file , "@\ts%d line color %d\n" , dataset , colorset ) ;  
  fprintf( file , "@\ts%d errorbar on\n" , dataset ) ;  
  fprintf( file , "@\ts%d errorbar place both\n" , dataset ) ;  
  fprintf( file , "@\ts%d errorbar color %d\n" , dataset , colorset ) ;  
  fprintf( file , "@\ts%d errorbar size 1\n" , dataset ) ; 
  return ;
}

static void
initialise_graph( const char *x_axis , 
		  const char *y_axis )
{
  fprintf( file , "@map font 0 to \"Times-Roman\", \"Times-Roman\"\n" ) ;
  fprintf( file , "@map font 1 to \"Symbol\", \"Symbol\"\n" ) ;
  fprintf( file , "@map font 2 to \"Helvetica\", \"Helvetica\"\n" ) ;
  fprintf( file , "@map color 0 to (255,255,255), \"white\"\n" ) ;
  fprintf( file , "@map color 1 to (0,0,0), \"black\"\n" ) ;
  fprintf( file , "@map color 2 to (255,0,0), \"red\"\n" ) ;
  fprintf( file , "@map color 3 to (0,255,0), \"green\"\n" ) ;
  fprintf( file , "@map color 4 to (0,0,255), \"blue\"\n" ) ;
  fprintf( file , "@map color 5 to (104,7,72), \"maroon\"\n" ) ;
  fprintf( file , "@map color 6 to (188,143,143), \"brown\"\n" ) ;
  fprintf( file , "@map color 7 to (255,165,0), \"orange\"\n" ) ;
  fprintf( file , "@map color 8 to (0,255,255), \"cyan\"\n" ) ;
  fprintf( file , "@default linestyle 1\n" ) ;
  fprintf( file , "@default linewidth 1.0\n" ) ;
  fprintf( file , "@default color 1\n" ) ;
  fprintf( file , "@default pattern 1\n" ) ;
  fprintf( file , "@default font 0\n" ) ;
  fprintf( file , "@background color 0\n" ) ;
  fprintf( file , "@page background fill on\n" ) ;
  fprintf( file , "@g0 on\n" ) ;
  fprintf( file , "@g0 type XY\n" ) ;
  fprintf( file , "@with g0\n" ) ;
  fprintf( file , "@\tview 0.180000, 0.180000, 1.250000, 0.950000\n" ) ;
  fprintf( file , "@\txaxes invert off\n" ) ;  
  fprintf( file , "@\tyaxes invert off\n" ) ;  
  fprintf( file , "@\txaxes scale Normal\n" ) ;  
  fprintf( file , "@\tyaxes scale Normal\n" ) ;  
  fprintf( file , "@\txaxis on\n" ) ;  
  fprintf( file , "@\txaxis\tlabel char size 1.490000\n" ) ;  
  fprintf( file , "@\txaxis\ttick minor ticks 0\n" ) ; 
  fprintf( file , "@\txaxis\ttick out\n" ) ;  
  fprintf( file , "@\tyaxis on\n" ) ;  
  fprintf( file , "@\tyaxis\tlabel char size 1.490000\n" ) ; 
  fprintf( file , "@\tyaxis\ttick minor ticks 0\n" ) ; 
  fprintf( file , "@\tyaxis\ttick out\n" ) ; 
  // X axis label
  fprintf( file , "@\txaxis label \"%s\" \n" , x_axis ) ;  
  // Y axis label
  fprintf( file , "@\tyaxis label \"%s\" \n" , y_axis ) ;  
  fprintf( file , "@\txaxis label place auto\n" ) ;  
  fprintf( file , "@\tyaxis label place auto\n" ) ;  
  fprintf( file , "@\taltxaxis\toff\n" ) ;  
  fprintf( file , "@\taltyaxis\toff\n" ) ;  
  fprintf( file , "@\tlegend on\n" ) ;  
  fprintf( file , "@\tlegend 0.75, 0.94\n" ) ;  
  fprintf( file , "@\tlegend char size 1.50000\n") ;  
  return ;
}

void
plot_data( const struct resampled *boots ,
	   const double *X ,
	   const int NDATA )
{
  // plot the data
  PrintHeader( ) ;
  fprintf( file , "@\ts%d legend \"\" \n" , dataset ) ;  
  fprintf( file , "@target G0.S%d\n" , dataset ) ;
  fprintf( file , "@type xydxdxdydy\n" ) ;

  // PACK the data for an average ...
  int i ;
  for( i = 0 ; i < NDATA ; i++ ) {
    fprintf( file , "%e %e 0.0 0.0 %e %e\n" ,
	     X[i] , boots[i].avg ,  
	     boots[i].err_hi - boots[i].avg ,
	     boots[i].avg - boots[i].err_lo ) ;
  }
  fprintf( file , "&\n" ) ;
  dataset ++ ;
  colorset = ( colorset + 1 ) % 15 ;
  return ;
}

void
plot_simple_data( const double *ydata ,
		  const double *xdata ,
		  const int NDATA ,
		  const char *legend )
{
  // plot the data
  PrintHeader( ) ;
  fprintf( file , "@\ts%d legend \"%s\" \n" , dataset , legend ) ;  
  fprintf( file , "@target G0.S%d\n" , dataset ) ;
  fprintf( file , "@type xy\n" ) ;
  // PACK the data for an average ...
  int i ;
  for( i = 0 ; i < NDATA ; i++ ) {
    fprintf( file , "%e %e \n" ,
	     xdata[i] , ydata[i] ) ;
  }
  fprintf( file , "&\n" ) ;
  dataset ++ ;
  colorset = ( colorset + 1 ) % 15 ;
  return ;
}

void
graph_reset_color( void )
{
  if( colorset != 1 ) colorset = 1 ;
}

void
graph_reset_symbols( void )
{
  if( dataset != 0 ) dataset = 1 ;
}

static void
draw_line( const double *X , 
	   const double *Y , 
	   const int length )
{
  fprintf( file , "@\ts%d hidden false\n" ,dataset ) ;  
  fprintf( file , "@\ts%d symbol 0\n" ,dataset ) ;  
  fprintf( file , "@\ts%d line color %d\n" , dataset , colorset ) ;  
  fprintf( file , "@\ts%d errorbar color %d\n" , dataset , colorset ) ;  
  fprintf( file , "@\ts%d legend \"\"\n" , dataset ) ;  
  fprintf( file , "@target G0.S%d\n" , dataset ) ;
  fprintf( file , "@type xy\n" ) ;

  int i ;
  for( i = 0 ; i <= length ; i++ ) {
    fprintf( file , "%1.12e %1.12e\n" , X[i] , Y[i] ) ;
  }
  fprintf( file , "&\n" ) ;
  dataset++ ;
  return ;
}

void
plot_fitfunc_french( const struct resampled *resample , 
		     const struct resampled *y ,
		     const double *x ,
		     const struct mom_info *mominfo ,
		     const fitfunc fit ,
		     const int NPARAMS ,
		     const struct chiral quarks ,
		     const double xmin ,
		     const double xmax ,
		     const int LT ,
		     const int NDATA )
{
  const int low = find_idx( xmin , x , NDATA , 0 ) ;
  const int upp = find_idx( xmax , x , NDATA , 0 ) ;
  const int NDIFF = upp-low ;

  double *X = malloc( NDIFF * sizeof( double ) ) ;
  double *Y = malloc( NDIFF * sizeof( double ) ) ;
  double *YMIN = malloc( NDIFF * sizeof( double ) ) ;
  double *YMAX = malloc( NDIFF * sizeof( double ) ) ;

  const int NBOOTS = resample[0].NSAMPLES ;

  int idx = 0 , i ;
  for( i = low ; i < upp ; i++ ) {

    struct resampled yy ;
    yy.resampled = (double*)malloc( NBOOTS * sizeof( double ) ) ;
    equate( &yy , y[i] ) ;

    // perform subtraction
    struct resampled temp ;
    temp.resampled = (double*)malloc( NBOOTS * sizeof( double ) ) ;
    equate( &temp , resample[3] ) ;
    mult_constant( &temp , mominfo[i].p4 / mominfo[i].p2 ) ;
    subtract( &yy , temp ) ;

    equate( &temp , resample[4] ) ;
    mult_constant( &temp , pow( mominfo[i].p4 / mominfo[i].p2 , 2 ) ) ;
    subtract( &yy , temp ) ;

    equate( &temp , resample[5] ) ;
    mult_constant( &temp , mominfo[i].p6 / mominfo[i].p2 ) ;
    subtract( &yy , temp ) ;

    equate( &temp , resample[6] ) ;
    mult_constant( &temp , mominfo[i].p4 ) ;
    subtract( &yy , temp ) ;

    X[ idx ] = x[ i ] ;

    Y[ idx ] = yy.avg ;
    YMAX[ idx ] = yy.err_hi ;
    YMIN[ idx ] = yy.err_lo ;

    free( temp.resampled ) ;
    free( yy.resampled ) ;

    idx ++ ;
  }

  // and draw the maximum and minimum of the fit
  draw_line( X , Y , NDIFF ) ;
  draw_line( X , YMIN , NDIFF ) ;
  draw_line( X , YMAX , NDIFF ) ;
  colorset = ( colorset + 1 ) % 15 ;

  free( X ) ;
  free( Y ) ;
  free( YMIN ) ;
  free( YMAX ) ;

  return ;
}

void
plot_fitfunc( const struct resampled *resample , 
	      const struct resampled *y ,
	      const double *x ,
	      const struct mom_info *mominfo ,
	      const fitfunc fit ,
	      const int NPARAMS ,
	      const struct chiral quarks ,
	      const double xmin ,
	      const double xmax ,
	      const int LT ,
	      const int NDATA ,
	      const int SLICE )
{
  printf( "minmax :: %f %f \n" , xmin , xmax ) ;

  const int low = find_idx( xmin , x , NDATA , 0 ) ;
  const int upp = find_idx( xmax , x , NDATA , 0 ) ;

  const int Q2_idx = find_idx( xmin , x , NDATA-1 , 0 ) ;
  set_q2q3_D0_diff_multi( mominfo[Q2_idx].p2 , 0 ) ;

  const int NDIFF = upp-low+1 ;

  double *X = malloc( NDIFF * sizeof( double ) ) ;
  double *Y = malloc( NDIFF * sizeof( double ) ) ;
  double *YMIN = malloc( NDIFF * sizeof( double ) ) ;
  double *YMAX = malloc( NDIFF * sizeof( double ) ) ;

  int idx = 0 , i ;
  for( i = low ; i <= upp ; i++ ) {

    // extrapolate the fit parameters
    struct resampled y = extrapolate( resample , x[i] ,
				      NPARAMS , quarks , 
				      mominfo[i] , fit , LT , SLICE ) ;

    X[ idx ] = x[ i ] ;
    Y[ idx ] = y.avg ;
    YMAX[ idx ] = y.err_hi ;
    YMIN[ idx ] = y.err_lo ;

    free( y.resampled ) ;

    idx ++ ;
  }

  // and draw the maximum and minimum of the fit
  draw_line( X , Y , NDIFF-1 ) ;
  draw_line( X , YMIN , NDIFF-1 ) ;
  draw_line( X , YMAX , NDIFF-1 ) ;
  colorset = ( colorset + 1 ) % 15 ;

  free( X ) ;
  free( Y ) ;
  free( YMIN ) ;
  free( YMAX ) ;

  return ;
}

void
plot_fitfunc2( const struct resampled *resample , 
	       const struct resampled *y ,
	       const double *x ,
	       const struct mom_info *mominfo ,
	       const fitfunc fit ,
	       const int NPARAMS ,
	       const struct chiral quarks ,
	       const double xmin ,
	       const double xmax ,
	       const int LT ,
	       const int NDATA ,
	       const int SLICE )
{
  printf( "minmax :: %f %f \n" , xmin , xmax ) ;


  const int granularity = 100 ;
  const double step = ( xmax - xmin ) / granularity ;
  double *X = malloc( granularity * sizeof( double ) ) ;
  double *Y = malloc( granularity * sizeof( double ) ) ;
  double *YMIN = malloc( granularity * sizeof( double ) ) ;
  double *YMAX = malloc( granularity * sizeof( double ) ) ;

  int idx = 0 ;

  double xx ;
  for( xx = xmin ; xx <= xmax ; xx+= step ) {

    // is a poor approximation
    struct mom_info tempmom ;

    // extrapolate the fit parameters
    struct resampled y = extrapolate( resample , xx ,
				      NPARAMS , quarks , 
				      tempmom , fit , LT , SLICE ) ;

    X[ idx ] = xx ;
    YMAX[ idx ] = y.err_hi ;
    YMIN[ idx ] = y.err_lo ;
    Y[ idx ] = y.avg ;

    idx ++ ;

    free( y.resampled ) ;
  }

  // and draw the maximum and minimum of the fit
  draw_line( X , Y , granularity-1 ) ;
  draw_line( X , YMIN , granularity-1 ) ;
  draw_line( X , YMAX , granularity-1 ) ;

  free( X ) ;
  free( Y ) ;
  free( YMIN ) ;
  free( YMAX ) ;

  colorset = ( colorset + 1 ) % 15 ;

  return ;
}


void
plot_chiral_fitfunc( const struct resampled *resample , 
		     const struct resampled *y ,
		     const double *x ,
		     const struct mom_info *mominfo ,
		     const fitfunc fit ,
		     const int NPARAMS ,
		     const struct chiral quarks ,
		     const double xmin ,
		     const double xmax ,
		     const int LT ,
		     const int NDATA ,
		     const double ainverse ,
		     const int SLICE )
{
  const int granularity = 100 ;
  const double step = ( xmax - xmin ) / granularity ;
  double *X = malloc( granularity * sizeof( double ) ) ;
  double *Y = malloc( granularity * sizeof( double ) ) ;
  double *YMIN = malloc( granularity * sizeof( double ) ) ;
  double *YMAX = malloc( granularity * sizeof( double ) ) ;

  int idx = 0 ;

  double xx ;
  for( xx = xmin ; xx <= xmax ; xx+= step ) {

    struct mom_info tempmom ;

    struct chiral tquark ;
    tquark.ml = xx ;
    tquark.ms = quarks.ZV ;
    tquark.ZA = quarks.ZA ;
    tquark.ZV = quarks.ZV ;
    tquark.ainverse = ainverse ;

    // extrapolate the fit parameters
    struct resampled y = extrapolate( resample , xx ,
				      NPARAMS , tquark , 
				      tempmom , fit , LT , SLICE ) ;

    compute_err( &y ) ;

    X[ idx ] = xx + quarks.ZV ;
    YMAX[ idx ] = y.err_hi ;
    YMIN[ idx ] = y.err_lo ;

    // and compute the average
    int k ;
    double data[ NPARAMS ] ;
    for( k = 0 ; k < NPARAMS ; k++ ) {
      data[ k ] = resample[ k ].avg ;
    }

    // CCCCCC-Callback for the fitfunction
    struct x_descriptor XX ; 
    XX.X = xx ; XX.quark = tquark ; XX.mom = tempmom ; XX.LT = LT ;
    Y[ idx ] = fit.f( data , XX , NPARAMS ) ;

    idx ++ ;

    free( y.resampled ) ;
  }

  // and draw the maximum and minimum of the fit
  draw_line( X , Y , granularity-1 ) ;
  draw_line( X , YMIN , granularity-1 ) ;
  draw_line( X , YMAX , granularity-1 ) ;

  free( X ) ;
  free( Y ) ;
  free( YMIN ) ;
  free( YMAX ) ;

  colorset = ( colorset + 1 ) % 15 ;

  return ;
}
// open
void
make_xmgrace_graph( const char *filename ,
		    const char *x_axis , 
		    const char *y_axis )
{
  printf( "\n[GRAPH] name :: %s \n" , filename ) ;
  printf( "[GRAPH] axes (x) %s (y) %s \n" , x_axis , y_axis ) ;
  file = fopen( filename , "w" ) ;
  initialise_graph( x_axis , y_axis ) ;
  return ;
}

// and close
void
close_xmgrace_graph( void )
{
  dataset = 0 ; colorset = 1 ;
  fclose( file ) ;
  return ;
}
