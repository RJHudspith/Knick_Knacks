
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static gsl_rng *r ;
static long unsigned int SEED = -1 ;

double
rng_double( void ) { return gsl_rng_uniform( r ) ; }

// between 0 and some number "max_idx"
int
rng_int( const int max_idx ) { 
  return gsl_rng_uniform_int( r , max_idx ) ;
}

double
rng_gaussian( const double sigma ) { 
  return gsl_ran_gaussian( r , sigma ) ; 
}

void
init_rng( long unsigned int Seed )
{
  gsl_rng_env_setup( ) ;
  
  const gsl_rng_type *type = gsl_rng_default ;
  r = gsl_rng_alloc (type) ;

  if( Seed == 0 ) {
    FILE *urandom = fopen( "/dev/urandom" , "r" ) ;
    if( urandom == NULL ) exit( 1 ) ;
    if( fread( &Seed , sizeof( Seed ) , 1 , urandom ) != 1 ) exit(1) ;
    fclose( urandom ) ;
  }

  SEED = Seed ;
  printf( "USING GSL seed %lu \n" , SEED ) ;

  gsl_rng_set( r , Seed ) ;

  return ;
}

void
rng_reseed( void ) 
{
  if( SEED != -1 ) {
    gsl_rng_set( r , SEED ) ;
  } else {
    printf( "GSL rng not seeded properly!\n" );
    exit(1) ;
  }
  return ;
}

void
free_rng( void )
{
  gsl_rng_free( r );
  return ;
}
