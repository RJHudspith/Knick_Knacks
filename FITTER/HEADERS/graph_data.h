#ifndef GRAPH_DATA_H
#define GRAPH_DATA_H

void
graph_reset_color( void ) ;

void
make_xmgrace_graph( const char *filename ,
		    const char *x_axis , 
		    const char *y_axis ) ;

void
close_xmgrace_graph( void ) ;

void
plot_simple_data( const double *ydata ,
		  const double *xdata ,
		  const int NDATA ,
		  const char *legend ) ;

void
plot_data( const struct resampled *boots ,
	   const double *X ,
	   const int NDATA ) ;


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
		     const int NDATA ,
		     const int SLICE ) ;
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
	      const int SLICE ) ;

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
	       const int SLICE ) ;

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
		     const int SLICE ) ;

#endif
