#include <mpi.h>
#include <iostream>
#include "Communication.hpp"
#include <vector>
#include "Fields.hpp"
#include "Datastructures.hpp"

MPI_Comm MPI_COMMUNICATOR;

/*
TODO -->  each process should store here information about its rank and communicator,
as well as its neighbors. If a process does not have any neighbor in some direction, simply
use MPI_PROC_NULL as rank to skip the respective communication steps
*/

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

void Communication::communicate(Matrix<double> &field){
    // Get my coordinates in the new communicator
    int my_coords[2];
    MPI_Cart_coords(MPI_COMMUNICATOR, my_rank_global, 2, my_coords);

    enum DIRECTIONS {DOWN, UP, LEFT, RIGHT};
    char* neighbours_names[4] = {"down", "up", "left", "right"};
    int neighbours_ranks[4];
 
    // check if neighbours to communicate with
    MPI_Cart_shift(MPI_COMMUNICATOR, 0, 1, &neighbours_ranks[LEFT], &neighbours_ranks[RIGHT]);
    MPI_Cart_shift(MPI_COMMUNICATOR, 1, 1, &neighbours_ranks[DOWN], &neighbours_ranks[UP]);
 
    // for(int i = 0; i < 4; i++)
    // {
    //     if(neighbours_ranks[i] == MPI_PROC_NULL)
    //         printf("[MPI process %d] I have no %s neighbour.\n", my_rank_global, neighbours_names[i]);
    //     else
    //         printf("[MPI process %d] I have a %s neighbour: process %d.\n", my_rank_global, neighbours_names[i], neighbours_ranks[i]);
    // }
 
    // for(auto &neighbour : neighbours_ranks){
    //     if(neighbour != MPI_PROC_NULL){
    //         MPI_Sendrecv(&field, );
    //     }
    // }

    MPI_Status status;

    // TODO --> execution gets stuck, potenetial deadlocks

    if(neighbours_ranks[LEFT]!= MPI_PROC_NULL){
        // std::cout << "COMM LEFT" << std::endl;
        std::vector<double> rcv_buffer, send_buffer;

        for(int j=0; j<field.num_rows(); j++){
            send_buffer.push_back(field(1,j));
            rcv_buffer.push_back(0);
        }

        int check = MPI_Sendrecv(send_buffer.data(), send_buffer.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0, 
                     rcv_buffer.data(), rcv_buffer.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,  MPI_COMMUNICATOR, &status);

        if (not(check)){
            std::cout << "MPI_Sendrecv failed " << std::endl;
        }

        for(int j=0; j<field.num_rows(); j++){
            std::cout << rcv_buffer.at(j) << " ";
            field(0,j) = rcv_buffer.at(j);
        }
        std::cout << std::endl;

        
    }

    // if(neighbours_ranks[RIGHT] != MPI_PROC_NULL){

    //     std::vector<double> rcv_buffer, send_buffer;
    //     int inner_index = field.num_cols() - 2;

    //     for(int j=0; j<field.num_rows(); j++){
    //         send_buffer.push_back(field(inner_index,j));
    //         // std::cout << field(inner_index,j) << " ";
    //         rcv_buffer.push_back(0);
    //     }
    //     // std::cout << "\n" << std::endl;

    //     // std::cout << "COMM RIGHT" << std::endl;             
    //     MPI_Sendrecv(send_buffer.data(), send_buffer.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0, 
    //                  rcv_buffer.data(), rcv_buffer.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,  MPI_COMMUNICATOR, &status);

    //     for(int j=0; j<field.num_rows(); j++){
    //         // std::cout << rcv_buffer.at(j) << " ";
    //         field(inner_index + 1,j) = rcv_buffer.at(j);
    //     }
    //     // std::cout << std::endl;
    // }

    // if(neighbours_ranks[UP]!= MPI_PROC_NULL){
    //     // std::cout << "COMM UP" << std::endl;
    //     MPI_Sendrecv(&field(1,field.num_rows() - 2), field.num_cols() - 2, MPI_DOUBLE, neighbours_ranks[UP], 0, 
    //                  &field(1,field.num_rows() - 1), field.num_cols() - 2, MPI_DOUBLE, neighbours_ranks[DOWN], 0,  MPI_COMMUNICATOR, &status);
    // }

    // if(neighbours_ranks[DOWN] != MPI_PROC_NULL){
    //     // std::cout << "COMM DOWN" << std::endl;
    //     MPI_Sendrecv(&field(1,1), field.num_cols() - 2, MPI_DOUBLE, neighbours_ranks[DOWN], 0, 
    //                  &field(1,0), field.num_cols() - 2, MPI_DOUBLE, neighbours_ranks[UP], 0,  MPI_COMMUNICATOR, &status);
    // }

}

double reduce_min(){
// // Finding the global Minimum
// // ...
// double localMin = *std::min_element(myVector.begin(),myVector.end());
// double globalMin;
// MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}

double reduce_sum(){

}