/**
   Coefficients for the description of the OPE

   nf = 3 for now
 */

// compute the OPE for just the vector for now
static inline double
c0VV_0loop( const double t ) { return ( 0.0253303 * t ) ; }
static inline double
c0VV_1loop( const double t ) { return ( 0.0253303 * t ) ; }
static inline double
c0VV_2loop( const double t ) { return t * ( 0.0415372 - 0.0284966 * t ) ; }
static inline double
c0VV_3loop( const double t ) { return t * ( 0.16138 + t * ( -0.144119 + 0.0427449 * t ) ) ; }

/* D0 OPE 
   t   is log( Q^2 / mu^2 ) : 
   a   is the coupling divided by Pi -> fit param p[1]
   der is the order of the derivative w.r.t a
 */
double
D0_OPE( const double t , const double a , const int der , const int NLOOPS ) {
  switch( NLOOPS ) {
  case 0 :
    switch( der ) {
    case 0 :
      return ( c0VV_0loop( t ) ) ;
    default :
      return 0.0 ;
    }
  case 1 :
    switch( der ) {
    case 0 : 
      return ( c0VV_0loop( t ) + a * ( c0VV_1loop( t ) ) ) ;
    case 1 : 
      return ( c0VV_1loop( t ) ) ;
    }
  case 2 :
    switch( der ) {
    case 0 : 
      return ( c0VV_0loop( t ) + a * ( c0VV_1loop( t ) + a * ( c0VV_2loop( t ) ) ) ) ;
    case 1 : 
      return ( c0VV_1loop( t ) + a * ( 2.0 * c0VV_2loop( t ) ) ) ;
    }
  default :
    switch( der ) {
    case 0 : 
      return ( c0VV_0loop( t ) + a * ( c0VV_1loop( t ) + a * ( c0VV_2loop( t ) + a * c0VV_3loop( t ) ) ) ) ;
    case 1 : 
      return ( c0VV_1loop( t ) + a * ( 2.0 * c0VV_2loop( t ) + 3.0 * a * c0VV_3loop( t ) ) ) ;
    }
  }
  return 0.0 ;
}

/**
   a^2 corrections
 */
double
D0_OPE_a2( const double t , const double a , const int der , const int NLOOPS ) {
  switch( NLOOPS ) {
  default :
    switch( der ) {
    case 0 : 
      return ( c0VV_0loop( t ) + a * ( c0VV_1loop( t ) + a * ( c0VV_2loop( t ) + a * c0VV_3loop( t ) ) ) ) / ( 12.0 * t ) ;
    case 1 : 
      return ( c0VV_1loop( t ) + a * ( 2.0 * c0VV_2loop( t ) + 3.0 * a * c0VV_3loop( t ) ) ) / ( 12.0 * t ) ;
    }
  }
  return 0.0 ;
}

// the vector terms ~ m^2 / q^2
static inline double
c2VV_0loop( const double t ) { return -0.151982 ; }
static inline double
c2VV_1loop( const double t ) { return -0.405285 + 0.30396 * t ; }
static inline double
c2VV_2loop( const double t ) { return -3.6424592 + t * ( +2.8749886 + t*-0.6459225 ) ; }

//  this one is proportional to mq/q^2 and uses the vector
//  I have put in the strange quark terms by hand
double
D2_V_OPE( const double t , const double a , const int der , const int NLOOPS , const double ml , const double ms ) {
  switch( NLOOPS ) {
  case 0 :
    switch( der ) {
    case 0 : 
      return ml * ml * ( c2VV_0loop( t ) ) ;
    case 1 : 
      return 0.0 ;
    }
  case 1 :
    switch( der ) {
    case 0 : 
      return ml * ml * ( c2VV_0loop( t ) + a * ( c2VV_1loop( t ) ) ) ;
    case 1 : 
      return ml * ml * ( c2VV_1loop( t ) ) ;
    }
  default :
    switch( der ) {
    case 0 : 
      return ml * ml * ( c2VV_0loop( t ) + a * ( c2VV_1loop( t ) + a * ( c2VV_2loop( t ) ) ) ) + a * a * ms * ms * 0.026602167 ;
    case 1 : 
      return ml * ml * ( c2VV_1loop( t ) + a * ( 2.0 * c2VV_2loop( t ) ) ) + 2.0 * a * ms * ms * 0.026602167 ;
    }
  }
  return 0.0 ;
}

// the GG condensate ...
static inline double
c4GG_1loop( const double t ) { return 1.0/12. ; }
static inline double
c4GG_2loop( const double t ) { return ( 11./18. ) /12. ; }


// D = 4 I split up into the gluon condensate and the quark condensate bit
double
D4_GG_OPE( const double t , const double a , const int der , const int NLOOPS ) {
  switch( NLOOPS ) {
  case 0 :
    return 0.0 ;
  case 1 :
    switch( der ) {
    case 0 : 
      return ( a * ( c4GG_1loop( t ) ) ) ;
    case 1 : 
      return ( c4GG_1loop( t ) ) ;
    case 2 : 
      return 0.0 ;
    }
  default :
    switch( der ) {
    case 0 : 
      return ( a * ( c4GG_1loop( t ) + a * c4GG_2loop( t ) ) ) ;
    case 1 : 
      return ( c4GG_1loop( t ) + 2.0 * a * c4GG_2loop( t ) ) ;
    case 2 : 
      return 2.0 * c4GG_2loop( t ) ;
    }
  }
  return 0.0 ;
}

// the q-qbar vector terms for the D=4 OPE
static inline double
c4qqV_0loop( const double t , const double ml , const double ms ) { return 2. * ml ; }
static inline double
c4qqV_1loop( const double t , const double ml , const double ms ) { return 0.96296296 * ml + 0.14814815 * ms ; }
static inline double
c4qqV_2loop( const double t , const double ml , const double ms ) { return ml * ( +13.1479 - 2.16667 * t ) + ms * ( 1.07394 - 0.333333 * t ) ; }

// and the vector case
double
D4_Vqq_OPE( const double t , const double a , const int der , const int NLOOPS , const double ml_ll , const double ms_ss )
{
  switch( NLOOPS ) {
  case 0 :
    switch( der ) {
    case 0 : return c4qqV_0loop( t , ml_ll , ms_ss ) ;
    case 1 :
      return 0.0 ;
    }
  case 1 :
    switch( der ) {
    case 0 : 
      return ( c4qqV_0loop( t , ml_ll , ms_ss ) + a * ( c4qqV_1loop( t , ml_ll , ms_ss ) ) ) ;
    case 1 : 
      return ( c4qqV_1loop( t , ml_ll , ms_ss ) ) ;
    }
  default :
    switch( der ) {
    case 0 : 
      return ( c4qqV_0loop( t , ml_ll , ms_ss ) + a * ( c4qqV_1loop( t , ml_ll , ms_ss ) + a * c4qqV_2loop( t , ml_ll , ms_ss ) ) ) ;
    case 1 : 
      return ( c4qqV_1loop( t , ml_ll , ms_ss ) + 2.0 * a * c4qqV_2loop( t , ml_ll , ms_ss ) ) ;
    }
  }
  return 0.0 ;
}


