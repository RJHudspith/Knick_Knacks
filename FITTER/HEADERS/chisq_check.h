#ifndef CHISQ_CHECK_H
#define CHISQ_CHECK_H

struct resampled
chisq_check( const struct resampled *fitparams ,
	     const struct resampled *ydata , 
	     const double *xdata , 
	     const int NPARAMS ,
	     const struct chiral *quarks ,
	     const struct mom_info *mominfo ,
	     const fitfunc fit ,
	     const int Ndata ,
	     const int LT ) ;

#endif
