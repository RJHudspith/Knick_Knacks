#ifndef FITDATA_H
#define FITDATA_H


struct resampled *
perform_bootfit( fitfunc *fit ,
		 int *NPARAMS ,
		 const struct resampled *boots , // the bootstrapped data 
		 const double *X , // the (sorted) x-axis data
		 const double *sigma , // the s.d of the y-data
		 const struct mom_info *mom_info ,
		 const struct chiral *quarks , // chiral data
		 const fit_type fittype , 
		 const int *NDATA , // the length of the X-array
		 const double FIT_HI ,
		 const double FIT_LO ,
		 const int SIMS ,
		 const int LT ,
		 const bool sim_params[ 12 ] ,
		 const int NCOMMON ) ;

#endif
