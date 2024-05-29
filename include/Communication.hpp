#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <mpi.h>
#include <iostream>

inline int my_rank_global;
inline int my_coords_global[2];

static void init_parallel(int argn, char **args){
    MPI_Init(&argn, &args);

    int iproc{*args[2]};
    int jproc{*args[3]};
    iproc = iproc - 48; // convert ascii to integer
    jproc = jproc - 48;

    int num_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    if (num_proc != iproc * jproc) {
        std::cerr << "Incompatible number of processors and domain decomposition!";
        return;
    }


    // Ask MPI to decompose our processes in a 2D cartesian grid for us
    int dims[2] = {iproc, jproc};
    MPI_Dims_create(num_proc, 2, dims);

    // Make both dimensions periodic
    int periods[2] = {false, false};

    // Let MPI assign arbitrary ranks if it deems it necessary
    int reorder = true;

    // Create a communicator given the 2D torus topology.
    MPI_Comm MPI_COMMUNICATOR;
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

static void finalize(){
    MPI_Finalize();
}

#endif