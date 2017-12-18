//  This header contains declarations of some message output functions,
//  not all mpi-ranks need to write meassages (most of them are the same).

#include <iostream>
#include <string>

#ifndef __VERBOSE__
#define __VERBOSE__

void MPI_cout( std::string mesg, int rank=0 );


#endif // __VRBOSE__
