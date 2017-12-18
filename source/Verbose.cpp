#include "Verbose.hpp"

void MPI_cout( std::string mesg, int rank ){
    if( rank == 0 ){
        std::cout << "==> " << mesg << std::endl;
    }
}
