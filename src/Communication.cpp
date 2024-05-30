#include <mpi.h>
#include <iostream>
#include "Communication.hpp"
#include "Fields.hpp"

MPI_Comm MPI_COMMUNICATOR;

Communication::Communication(int argn, char **args){
    init_parallel(argn, args);
}

void Communication::init_parallel(int argn, char **args){
    MPI_Init(&argn, &args);

    // initialized to sequential execution
    int iproc{1};
    int jproc{1};

    if (argn > 2){
        iproc = *args[2] - 48;  //convert to ASCII
        jproc = *args[3] - 48;
    }

    int num_proc; // total number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    if (num_proc != iproc * jproc) {
        std::cerr << "Incompatible number of processors and domain decomposition!";
        MPI_Finalize();
        return;
    }

    // Ask MPI to decompose our processes in a 2D cartesian grid for us
    int dims[2] = {iproc, jproc};
    MPI_Dims_create(num_proc, 2, dims);     // #processes - #dimensions of cartesian grid - #dimensions vector

    // both dimensions are not periodic
    int periods[2] = {false, false};

    // Let MPI assign arbitrary ranks if it deems it necessary
    int reorder = true;

    // Create a communicator given the 2D torus topology.
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &MPI_COMMUNICATOR);

    // My rank in the new communicator
    int my_rank;
    MPI_Comm_rank(MPI_COMMUNICATOR, &my_rank);
    my_rank_global = my_rank;

    // Get my coordinates in the new communicator
    int my_coords[2];
    MPI_Cart_coords(MPI_COMMUNICATOR, my_rank, 2, my_coords);
    my_coords_global[0] = my_coords[0];
    my_coords_global[1] = my_coords[1];

    // Print my location in the 2D torus.
    printf("[MPI process %d] I am located at (%d, %d).\n", my_rank, my_coords[0],my_coords[1]);
}

void Communication::finalize(){
    MPI_Finalize();
}

void Communication::communicate(Fields &field){
    // Get my coordinates in the new communicator
    int my_coords[2];
    MPI_Cart_coords(MPI_COMMUNICATOR, my_rank_global, 2, my_coords);

    int left, right, up, down;
    MPI_Cart_shift(MPI_COMMUNICATOR, 0, 1, &left, &right);
    MPI_Cart_shift(MPI_COMMUNICATOR, 1, 1, &down, &up);

}