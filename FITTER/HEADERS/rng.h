#ifndef RNG_H
#define RNG_H

double
rng_double( void ) ;
int
rng_int( const int max_idx ) ;
double
rng_gaussian( const double sigma ) ;
void
init_rng( long unsigned int Seed ) ;
void
rng_reseed( void ) ;
void
free_rng( void ) ;

#endif
