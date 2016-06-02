#ifndef DERIVATIVES_H
#define DERIVATIVES_H

struct resampled
fit_der( const struct resampled *fitparams ,
	 const double extrap_point ,
	 const int NPARAMS ,
	 const struct chiral quarks ,
	 const struct mom_info mominfo ,
	 const fitfunc fit ,
	 const int LT ,
	 const int SLICE ) ;

#endif
