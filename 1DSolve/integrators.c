/**
   @file integrators.c
   @brief file where I keep all the integrators

   TODO :: compute the implicit numerical factors by 
   hand into a LUT - J
 */
#include "common.h" // common definitions

// some numerical constants for gaussian quadrature
#define r3 (1.7320508075688772)
#define r6 (2.449489742783178)
#define r15 (3.872983346207417)

// maximum fixed-point iterations
#define FXMAX (50)

// fixed point iteration tolerance
#define FXTOL (1E-12)

// trivial euler solution
static inline double
euler_step( double (*f)( const point P ) ,
	    const double h ,
	    const point P ) 
{
  return h * f( P ) ;
}

// backward implicit euler iteration
static inline double
backward_euler_step( double (*f)( const point P ) ,
		     const double h ,
		     const point P )
{
  // fixed-point iterative solution for K1
  size_t n = 0 ;
  double K1 = f( P ) , res = K1 , err = 1.0 ;
  while( err > FXTOL && n < FXMAX ) {
    res = f( (point){ P.x , P.y + h*K1 } ) ;
    err = fabs( K1 - res ) ;
    K1 = res ;
  }
  // tell us if we err
  if( n == FXMAX ) {
    printf( "Backward Euler fixed point iteration failure :: %e\n" , err ) ;
  }
  return h * K1 ;
}

// explicit midpoint method (RK2)
static inline double
midpoint_step( double (*f)( const point P ) ,
	       const double h ,
	       const point P ) 
{
  return h * f( (point){ P.x + h/2. , P.y + h/2.*f( P ) } ) ;
}

// implicit midpoint integration
static inline double
implicit_midpoint_step( double (*f)( const point P ) ,
			const double h ,
			const point P ) 
{
  size_t n = 0 ;
  double K1 = f( P ) , res , err = 1.0 ;
  while( err > FXTOL && n < FXMAX ) {
    res = f( (point){ P.x + h/2. , P.y + (h/2.)*K1 } ) ;
    err = fabs( K1 - res ) ;
    K1  = res ;
    n++ ;
  }
  // tell us if we err
  if( n == FXMAX ) {
    printf( "Implicit midpoint fixed point iteration failure :: %e\n" , err ) ;
  }
  return h*K1 ;
}

// heun's method
static inline double
heun_step( double (*f)( const point P ) ,
	   const double h ,
	   const point P )
{
  return 0.5 * h * ( f( P ) + 
		     f( (point){ P.x+h , P.y + h*f( P ) } ) ) ;
}

// at this level it is the same as a radau 1A third order method
static inline double
ralston_step( double (*f)( const point P ) ,
	      const double h ,
	      const point P )
{
  return (h/4.) * ( f( P ) + 
		   3*f( (point){ P.x + h*2./3. , 
			 P.y + (h*2./3.)*f( P ) } ) ) ;
}

// RK3
static inline double
RK3_step( double (*f)( const point P ) ,
	  const double h ,
	  const point P )
{
  register const double K1 = f( P ) ;
  register const double K2 = f( (point){ P.x + h/2. , P.y + (h/2.)*K1 } ) ;
  register const double K3 = f( (point){ P.x + h , P.y + h*( -K1 + 2*K2 ) } ) ; 
  return (h/6.0) * ( K1 + 4*K2 + K3 ) ;
}

// RK4
static inline double
RK4_step( double (*f)( const point P ) ,
	  const double h ,
	  const point P )
{
  register const double K1 = f( P ) ;
  register const double K2 = f( (point){ P.x + h/2. , P.y + (h/2.)*K1 } ) ;
  register const double K3 = f( (point){ P.x + h/2. , P.y + (h/2.)*K2 } ) ;
  register const double K4 = f( (point){ P.x + h , P.y + h*K3 } ) ;
  return (h/6.0) * ( K1 + 2*(K2+K3) + K4 ) ;
} 

// RK4
static inline double
RK4_38_step( double (*f)( const point P ) ,
	     const double h ,
	     const point P )
{
  register const double K1 = f( P ) ;
  register const double K2 = f( (point){ P.x + h/3. , P.y + (h/3.)*K1 } ) ;
  register const double K3 = f( (point){ P.x + (2./3.)*h , 
	                                 P.y + h/3.*( -K1 + 3*K2 ) } ) ;
  register const double K4 = f( (point){ P.x + h , 
	                                 P.y + h*( K1 - K2 + K3 ) } ) ;
  return (h/8.0) * ( K1 + 3*( K2 + K3 ) + K4 ) ;
}

// gaussian quadratures are special and need a fixed point solve
static inline double
gauss4_step( double (*f)( const point P ) ,
	     const double h ,
	     const point P )
{
  size_t n = 0 ;
  double K1 = f( P ) , res1 , res2 , err = 1.0 ;
  double K2 = f( (point){ P.x + h/2. , P.y + (h/2.)*K1 } ) ; 
  while( err > FXTOL && n < FXMAX ) {
    // first term
    res1 = f( (point){ P.x+(h/6.)*( 3 - r3 ) , 
	  P.y + (h/12.)*( 3 * K1 + K2 * ( 3 - 2 * r3 ) ) } ) ;
    // second term
    res2 = f( (point){ P.x+(h/6.)*( 3 + r3 ) , 
	  P.y + (h/12.)*( K1 * ( 3 + 2 * r3 ) + 3 * K2 ) } ) ;
    // compute the error
    err = fabs( K1 - res1 + K2 - res2 ) ;
    K1 = res1 ; K2 = res2 ;
    n++ ;
  }
  // tell us if we err
  if( n == FXMAX ) {
    printf( "Gauss4 fixed point iteration failure :: %e\n" , err ) ;
  }
  return (h/2.)*( K1 + K2 ) ;
}

// fifth order implicit gaussian step
static inline double
gauss5_step( double (*f)( const point P ) ,
	     const double h ,
	     const point P )
{
  size_t n = 0 ;
  double res1 , res2 , res3 , err = 1.0 ;
  double K1 = f( P ) ;
  double K2 = f( (point){ P.x + h/2. , P.y + (h/2.)*K1 } ) ; 
  double K3 = f( (point){ P.x + h , P.y + (h)*K2 } ) ; 
  while( err > FXTOL && n < FXMAX ) {
    // first row
    res1 = f( (point){ P.x+(h/10.)*( 5 - r15 ) , 
       P.y + h*( (5/36.)*K1 + (2/9. - r15/15.)*K2 + (5/36.-r15/30.)*K3 ) } ) ;
    // second row
    res2 = f( (point){ P.x + h/2. , 
       P.y + h*( (5/36.+r15/24.)*K1 + (2/9.)*K2 + (5/36.-r15/24.)*K3 ) } ) ;
    // third row
    res3 = f( (point){ P.x+(h/10.)*( 5 + r15 ) , 
       P.y + h*( (5/36.+r15/30.)*K1 + (2/9.+r15/15)*K2 + (5/36.)*K3  ) } ) ;
    // compute error
    err = fabs( K1 - res1 + K2 - res2 + K3 - res3 ) ;
    K1 = res1 ; K2 = res2 ; K3 = res3 ;
    n++ ;
  }
  // to err is human, to forgive divine
  if( n == FXMAX ) {
    printf( "Gauss5 fixed point iteration failure :: %e\n" , err ) ;
  }
  return (h/18.)*( 5 * K1 + 8 * K2 + 5 * K3 ) ;
}

// randau method of type 1a
static inline double
radau3_step( double (*f)( const point P ) ,
	     const double h ,
	     const point P )
{
  size_t n = 0 ;
  double err = 1.0 , res1 , res2 ;
  double K1 = f( P ) ;
  double K2 = f( (point){ P.x + (h/2.) , P.y + (h/2.)*K1 } ) ;
  while( err > FXTOL && n < FXMAX ) {
    // first row of the tableau
    res1 = f( (point){ P.x + h/3. , P.y + (h/12.)*( 5*K1 - K2 ) } ) ;
    // second row
    res2 = f( (point){ P.x + h , P.y + (h/4.)*( 3*K1 + K2 ) } ) ;
    // compute error
    err = fabs( K1 - res1 + K2 - res2 ) ;
    K1 = res1 ; K2 = res2 ;
    n++ ;
  }
  // to err is human
  if( n == FXMAX ) {
    printf( "Radau3 fixed point iteration failure :: %e\n" , err ) ;
  }
  return (h/4.)*( 3*K1 + K2 ) ;
}

// fifth order randau method of type1a
static inline double
radau5_step( double (*f)( const point P ) ,
	     const double h ,
	     const point P )
{
  size_t n = 0 ;
  double res1 , res2 , res3 , err = 1.0 ;
  double K1 = f( P ) ;
  double K2 = f( (point){ P.x + h/2. , P.y + (h/2.)*K1 } ) ; 
  double K3 = f( (point){ P.x + h , P.y + (h)*K2 } ) ; 
  while( err > FXTOL && n < FXMAX ) {
    // first row
    res1 = f( (point){ P.x + (h/10.)*( 4 - r6 ) , 
	  P.y + h*( (11/45.-7*r6/360.)*K1 + (37/225.-169*r6/1800.)*K2 + 
		    (-2/225.+r6/75.)*K3 ) } ) ;
    // second row
    res2 = f( (point){ P.x + (h/10.)*( 4 + r6 ) , 
	  P.y + h*( (37/225.+169*r6/1800.)*K1 + (11/45.+7*r6/360.)*K2 +
		    (-2/225.-r6/75.)*K3 ) } ) ;
    // third row
    res3 = f( (point){ P.x + h , 
	  P.y + h*( (4/9.-r6/36.)*K1 + (4/9.+r6/36.)*K2 + (1/9.)*K3 ) } ) ;
    // compute the error
    err = fabs( K1 - res1 + K2 - res2 + K3 - res3 ) ;
    K1 = res1 ; K2 = res2 ; K3 = res3 ;
    n++ ;
  }
  // to err is human, to forgive divine
  if( n == FXMAX ) {
    printf( "radau5 fixed point iteration failure :: %e\n" , err ) ;
  }
  return (h/36.)*( ( 16-r6 ) * K1 + ( 16+r6 ) * K2 + 4 * K3 ) ;
}

// initialise our integration scheme
Integrator
initialise_integrator( const integration_scheme schema , 
		       const double tolerance )
{
  Integrator integrator ;
  integrator.nsteps = 0 ;
  integrator.notoksteps = 0 ;
  switch( schema ) {
  case EULER : 
    integrator.step_fwd = euler_step ; 
    integrator.adaptive_growth = -0.04 ; 
    integrator.adaptive_shrink = -0.04 ; 
    integrator.maxsteps = 100 ;
    break ;
  case BACKWARD_EULER : 
    integrator.step_fwd = backward_euler_step ; 
    integrator.adaptive_growth = -0.5 ; 
    integrator.adaptive_shrink = -0.4 ; 
    integrator.maxsteps = 100 ;
    break ; 
  case MIDPOINT : 
    integrator.step_fwd = midpoint_step ; 
    integrator.adaptive_growth = -0.15 ; 
    integrator.adaptive_shrink = -0.12 ; 
    integrator.maxsteps = 50 ;
    break ; 
  case IMPLICIT_MIDPOINT : 
    integrator.step_fwd = implicit_midpoint_step ; 
    integrator.adaptive_growth = -0.4 ; 
    integrator.adaptive_shrink = -0.35 ; 
    integrator.maxsteps = 50 ;
    break ; 
  case HEUN : 
    integrator.step_fwd = heun_step ; 
    integrator.adaptive_growth = -0.09 ; 
    integrator.adaptive_shrink = -0.06 ; 
    integrator.maxsteps = 50 ;
    break ; 
  case RALSTON : 
    integrator.step_fwd = ralston_step ; 
    integrator.adaptive_growth = -0.12 ; 
    integrator.adaptive_shrink = -0.08 ; 
    integrator.maxsteps = 30 ;
    break ; 
  case RK3 : 
    integrator.step_fwd = RK3_step ; 
    integrator.adaptive_growth = -0.08 ; 
    integrator.adaptive_shrink = -0.05 ; 
    integrator.maxsteps = 20 ;
    break ;
  case RK4 : 
    integrator.step_fwd = RK4_step ; 
    integrator.adaptive_growth = -0.08 ; 
    integrator.adaptive_shrink = -0.05 ; 
    integrator.maxsteps = 20 ;
    break ;
  case RK4_38 : 
    integrator.step_fwd = RK4_38_step ; 
    integrator.adaptive_growth = -0.08 ; 
    integrator.adaptive_shrink = -0.05 ; 
    integrator.maxsteps = 20 ;
    break ;
  case GAUSS4 : // implicit gaussian quadrature
    integrator.step_fwd = gauss4_step ; 
    integrator.adaptive_growth = -0.04 ; 
    integrator.adaptive_shrink = -0.03 ; 
    integrator.maxsteps = 20 ;
    break ;
  case GAUSS5 : // implicit gaussian quadrature
    integrator.step_fwd = gauss5_step ; 
    integrator.adaptive_growth = -0.01 ; 
    integrator.adaptive_shrink = -0.0075 ; 
    integrator.maxsteps = 20 ;
    break ;
  case RADAU3 : // implicit radau type1a method
    integrator.step_fwd = radau3_step ; 
    integrator.adaptive_growth = -0.2 ; 
    integrator.adaptive_shrink = -0.1 ; 
    integrator.maxsteps = 20 ;
    break ;
  case RADAU5 : // implicit radau type1a method
    integrator.step_fwd = radau5_step ; 
    integrator.adaptive_growth = -0.04 ; 
    integrator.adaptive_shrink = -0.03 ; 
    integrator.maxsteps = 20 ;
    break ;
  }
  integrator.adaptive_safety = 0.9 ;
  integrator.adaptive_errcon = powl( 5.0/integrator.adaptive_safety , 
				     1.0/integrator.adaptive_growth ) ;
  integrator.tolerance = tolerance ;
  return integrator ;
}
