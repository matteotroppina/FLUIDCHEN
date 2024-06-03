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
 

    MPI_Status status;

    // TODO --> execution gets stuck, potenetial deadlocks
    std::vector<double> buffer_left;
    std::vector<double> buffer_right;
    std::vector<double> buffer_up;
    std::vector<double> buffer_down;
    int inner_index_cols = field.num_cols() - 2;
    int inner_index_rows = field.num_rows() - 2;
    
    for(int j=0; j<field.num_rows(); j++){
            buffer_left.push_back(field(1,j));
            buffer_right.push_back(field(inner_index_cols,j));
    }
    for(int i=0; i<field.num_cols(); i++){
            buffer_up.push_back(field(i,1));
            buffer_down.push_back(field(i,inner_index_rows));
    }

    if(neighbours_ranks[RIGHT]!= MPI_PROC_NULL){
        if(neighbours_ranks[LEFT]= MPI_PROC_NULL){
        MPI_Sendrecv(buffer_right.data(), buffer_right.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0, 
                     buffer_left.data(), buffer_left.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,  MPI_COMMUNICATOR, &status);
            for(int j=0; j<field.num_rows(); j++){
                field(inner_index_cols+1,j) = buffer_left.at(j);
            }        
        }
    }

    if(neighbours_ranks[RIGHT]= MPI_PROC_NULL){
        if(neighbours_ranks[LEFT]!= MPI_PROC_NULL){
        MPI_Sendrecv(buffer_left.data(), buffer_left.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0, 
                     buffer_right.data(), buffer_right.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,  MPI_COMMUNICATOR, &status);
            for(int j=0; j<field.num_rows(); j++){  
                field(0,j) = buffer_right.at(j);
            }
        }
    } 

    if(neighbours_ranks[UP]!= MPI_PROC_NULL){
        if(neighbours_ranks[DOWN]= MPI_PROC_NULL){
        MPI_Sendrecv(buffer_up.data(), buffer_up.size(), MPI_DOUBLE, neighbours_ranks[UP], 0, 
                     buffer_down.data(), buffer_down.size(), MPI_DOUBLE, neighbours_ranks[UP], 0,  MPI_COMMUNICATOR, &status);
            for(int i=0; i<field.num_cols(); i++){
                field(i,0) = buffer_down.at(i);
            }        
        }
    }

    if(neighbours_ranks[UP]= MPI_PROC_NULL){
        if(neighbours_ranks[DOWN]!= MPI_PROC_NULL){
        MPI_Sendrecv(buffer_down.data(), buffer_down.size(), MPI_DOUBLE, neighbours_ranks[DOWN], 0, 
                     buffer_up.data(), buffer_up.size(), MPI_DOUBLE, neighbours_ranks[DOWN], 0,  MPI_COMMUNICATOR, &status);
            for(int i=0; i<field.num_cols(); i++){  
                field(i,inner_index_rows+1) = buffer_up.at(i);
            }
        }
    }    

}


//double reduce_min(){
    // // Finding the global Minimum
// // ...
// double localMin = *std::min_element(myVector.begin(),myVector.end());
// double globalMin;
// MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//}



//double reduce_sum(){

//}