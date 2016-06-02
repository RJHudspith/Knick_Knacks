#ifndef FAKE_DATA_H
#define FAKE_DATA_H

double *
fake_X( const int NDATA , const double RANGE ) ;

struct resampled **
fake_data( double ***X ,
	   const double XRANGE ,
	   const fit_type fittype ,
	   const struct chiral *quarks ,
	   const int NSLICES ,
	   int *NDATA ,
	   const int NRAW ,
	   const int LT ,
	   const bool sim_params[12] ) ;

#endif
