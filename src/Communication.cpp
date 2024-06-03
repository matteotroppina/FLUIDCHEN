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
    int inner_index_cols = field.num_cols() - 2;
    int inner_index_rows = field.num_rows() - 2;

    int num_proc; // total number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    std::vector<double> send1 (field.num_rows());
    std::vector<double> rcv1 (field.num_rows());


    if(neighbours_ranks[RIGHT]!= MPI_PROC_NULL){
        if(neighbours_ranks[LEFT]= MPI_PROC_NULL){
            for(int j=0; j<field.num_rows(); j++){
                send1[j] = field(inner_index_cols,j);
            }
            MPI_Sendrecv(&send1[0], send1.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,
                         &rcv1[0], rcv1.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<field.num_rows(); j++){
                field(inner_index_cols+1,j) = rcv1[j];
            }
        }
    }

    if(neighbours_ranks[LEFT]!= MPI_PROC_NULL){
        if(neighbours_ranks[RIGHT]= MPI_PROC_NULL){
            for(int j=0; j<field.num_rows(); j++){
                send1[j] = field(1,j);
            }
            MPI_Sendrecv(&send1[0], send1.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,
                         &rcv1[0], rcv1.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<field.num_rows(); j++){
                field(0,j) = rcv1[j];
            }
        }
    }
    
}

//     if(neighbours_ranks[RIGHT]!= MPI_PROC_NULL){
//         if(neighbours_ranks[LEFT]= MPI_PROC_NULL){
//         send[1] = 5;
//         MPI_Sendrecv(&send[0], send.size(), MPI_INT, neighbours_ranks[RIGHT], 0,
//                      &receive[0], receive.size(), MPI_INT, neighbours_ranks[RIGHT], 0,  MPI_COMMUNICATOR, &status);
//         }
//         std::cout << my_rank_global << ".rank" << " receive = " << receive[0] << " " << receive[1] << std::endl;
//     }
//     if(neighbours_ranks[LEFT]!= MPI_PROC_NULL){
//         if(neighbours_ranks[RIGHT]= MPI_PROC_NULL){
//         send[1] = 6;
//         MPI_Sendrecv(&send[0], send.size(), MPI_INT, neighbours_ranks[LEFT], 0,
//                      &receive[0], receive.size(), MPI_INT, neighbours_ranks[LEFT], 0,  MPI_COMMUNICATOR, &status);
//         }
//         std::cout << my_rank_global << ".rank" << " receive = " << receive[0] << " " << receive[1] << std::endl;
//     }

    //TODO -- make also num proc global
//     std::vector<std::vector<double>> all_left (num_proc);
//     std::vector<std::vector<double>> all_right (num_proc);
//     std::vector<std::vector<double>> all_up (num_proc);
//     std::vector<std::vector<double>> all_down (num_proc);

//     for(int j=0; j<field.num_rows(); j++){
//             all_left[my_rank_global].push_back(field(1,j));
//             all_right[my_rank_global].push_back(field(inner_index_cols,j));
//     }
//     for(int i=0; i<field.num_cols(); i++){
//             all_up[my_rank_global].push_back(field(i,1));
//             all_down[my_rank_global].push_back(field(i,inner_index_rows));
//     }
//     if(neighbours_ranks[RIGHT]!= MPI_PROC_NULL){
//         if(neighbours_ranks[LEFT]= MPI_PROC_NULL){
//         MPI_Sendrecv(&all_right[my_rank_global], all_right[my_rank_global].size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0, 
//                      &all_left[neighbours_ranks[RIGHT]], all_left[neighbours_ranks[RIGHT]].size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,  MPI_COMMUNICATOR, &status);
//             for(int j=0; j<field.num_rows(); j++){
//                 //field(inner_index_cols+1,j) = all_left[neighbours_ranks[RIGHT]].at(j);
//                 field(inner_index_cols+1,j) = all_left[neighbours_ranks[RIGHT]][j];
//             }        
//         }
//     }
//     if(neighbours_ranks[RIGHT]= MPI_PROC_NULL){
//         if(neighbours_ranks[LEFT]!= MPI_PROC_NULL){
//         MPI_Sendrecv(&all_left[my_rank_global], all_left[my_rank_global].size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0, 
//                      &all_right[neighbours_ranks[LEFT]], all_right[neighbours_ranks[LEFT]].size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,  MPI_COMMUNICATOR, &status);
//             for(int j=0; j<field.num_rows(); j++){  
//                 field(0,j) = all_right[neighbours_ranks[LEFT]][j];
//             }
//         }
//     } 

//     if(neighbours_ranks[UP]!= MPI_PROC_NULL){
//         if(neighbours_ranks[DOWN]= MPI_PROC_NULL){
//         MPI_Sendrecv(&all_up[my_rank_global], all_up[my_rank_global].size(), MPI_DOUBLE, neighbours_ranks[UP], 0, 
//                      &all_down[neighbours_ranks[UP]], all_down[neighbours_ranks[UP]].size(), MPI_DOUBLE, neighbours_ranks[UP], 0,  MPI_COMMUNICATOR, &status);
//             for(int i=0; i<field.num_cols(); i++){
//                 field(i,0) = all_down[neighbours_ranks[UP]][i];
//             }        
//         }
//     }

//     if(neighbours_ranks[UP]= MPI_PROC_NULL){
//         if(neighbours_ranks[DOWN]!= MPI_PROC_NULL){
//         MPI_Sendrecv(&all_down[my_rank_global], all_down[my_rank_global].size(), MPI_DOUBLE, neighbours_ranks[DOWN], 0, 
//                      &all_up[neighbours_ranks[DOWN]],  all_up[neighbours_ranks[DOWN]].size(), MPI_DOUBLE, neighbours_ranks[DOWN], 0,  MPI_COMMUNICATOR, &status);
//             for(int i=0; i<field.num_cols(); i++){  
//                 field(i,inner_index_rows+1) =  all_up[neighbours_ranks[DOWN]][i];
//             }
//         }
//     }    

// }


//double reduce_min(){
    // // Finding the global Minimum
// // ...
// double localMin = *std::min_element(myVector.begin(),myVector.end());
// double globalMin;
// MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//}



//double reduce_sum(){

//}