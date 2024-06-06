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

void Communication::communicateU(Fields &field, Matrix<double> &matrix){
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
    int inner_index_cols = matrix.num_cols() - 2;
    int inner_index_rows = matrix.num_rows() - 2;
    std::vector<double> send1 (matrix.num_rows());
    std::vector<double> rcv1 (matrix.num_rows());

    if(neighbours_ranks[RIGHT]!= MPI_PROC_NULL){
        if(neighbours_ranks[LEFT] == MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send1[j] = field.u(inner_index_cols,j);
            }
            MPI_Sendrecv(&send1[0], send1.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,
                         &rcv1[0], rcv1.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.u(inner_index_cols+1,j) = rcv1[j];
            }
        }
    }
    if(neighbours_ranks[LEFT]!= MPI_PROC_NULL){
        if(neighbours_ranks[RIGHT] == MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send1[j] = field.u(1,j);
            }
            MPI_Sendrecv(&send1[0], send1.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,
                         &rcv1[0], rcv1.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.u(0,j) = rcv1[j];
            }
        }
    }
}

void Communication::communicateV(Fields &field, Matrix<double> &matrix){
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
    int inner_index_cols = matrix.num_cols() - 2;
    int inner_index_rows = matrix.num_rows() - 2;
    std::vector<double> send2 (matrix.num_rows());
    std::vector<double> rcv2 (matrix.num_rows());

    if(neighbours_ranks[RIGHT]!= MPI_PROC_NULL){
        if(neighbours_ranks[LEFT]= MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send2[j] = field.v(inner_index_cols,j);
            }
            MPI_Sendrecv(&send2[0], send2.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,
                         &rcv2[0], rcv2.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.v(inner_index_cols+1,j) = rcv2[j];
            }
        }
    }
    if(neighbours_ranks[LEFT]!= MPI_PROC_NULL){
        if(neighbours_ranks[RIGHT]= MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send2[j] = field.v(1,j);
            }
            MPI_Sendrecv(&send2[0], send2.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,
                         &rcv2[0], rcv2.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.v(0,j) = rcv2[j];
            }
        }
    }
}

void Communication::communicateT(Fields &field, Matrix<double> &matrix){
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
    int inner_index_cols = matrix.num_cols() - 2;
    int inner_index_rows = matrix.num_rows() - 2;
    std::vector<double> send3 (matrix.num_rows());
    std::vector<double> rcv3 (matrix.num_rows());

    if(neighbours_ranks[RIGHT]!= MPI_PROC_NULL){
        if(neighbours_ranks[LEFT]= MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send3[j] = field.T(inner_index_cols,j);
            }
            MPI_Sendrecv(&send3[0], send3.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,
                         &rcv3[0], rcv3.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.T(inner_index_cols+1,j) = rcv3[j];
            }
        }
    }
    if(neighbours_ranks[LEFT]!= MPI_PROC_NULL){
        if(neighbours_ranks[RIGHT]= MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send3[j] = field.T(1,j);
            }
            MPI_Sendrecv(&send3[0], send3.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,
                         &rcv3[0], rcv3.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.T(0,j) = rcv3[j];
            }
        }
    }
}

void Communication::communicateP(Fields &field, Matrix<double> &matrix){
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
    int inner_index_cols = matrix.num_cols() - 2;
    int inner_index_rows = matrix.num_rows() - 2;
    std::vector<double> send4 (matrix.num_rows());
    std::vector<double> rcv4 (matrix.num_rows());

    if(neighbours_ranks[RIGHT]!= MPI_PROC_NULL){
        if(neighbours_ranks[LEFT]= MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send4[j] = field.p(inner_index_cols,j);
            }
            MPI_Sendrecv(&send4[0], send4.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,
                         &rcv4[0], rcv4.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.p(inner_index_cols+1,j) = rcv4[j];
            }
        }
    }
    if(neighbours_ranks[LEFT]!= MPI_PROC_NULL){
        if(neighbours_ranks[RIGHT]= MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send4[j] = field.p(1,j);
            }
            MPI_Sendrecv(&send4[0], send4.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,
                         &rcv4[0], rcv4.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.p(0,j) = rcv4[j];
            }
        }
    }
}

void Communication::communicateF(Fields &field, Matrix<double> &matrix){
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
    int inner_index_cols = matrix.num_cols() - 2;
    int inner_index_rows = matrix.num_rows() - 2;
    std::vector<double> send1 (matrix.num_rows());
    std::vector<double> rcv1 (matrix.num_rows());

    if(neighbours_ranks[RIGHT]!= MPI_PROC_NULL){
        if(neighbours_ranks[LEFT]= MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send1[j] = field.f(inner_index_cols,j);
            }
            MPI_Sendrecv(&send1[0], send1.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,
                         &rcv1[0], rcv1.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.f(inner_index_cols+1,j) = rcv1[j];
            }
        }
    }
    if(neighbours_ranks[LEFT]!= MPI_PROC_NULL){
        if(neighbours_ranks[RIGHT]= MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send1[j] = field.f(1,j);
            }
            MPI_Sendrecv(&send1[0], send1.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,
                         &rcv1[0], rcv1.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.f(0,j) = rcv1[j];
            }
        }
    }
}

void Communication::communicateG(Fields &field, Matrix<double> &matrix){
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
    int inner_index_cols = matrix.num_cols() - 2;
    int inner_index_rows = matrix.num_rows() - 2;
    std::vector<double> send1 (matrix.num_rows());
    std::vector<double> rcv1 (matrix.num_rows());

    if(neighbours_ranks[RIGHT]!= MPI_PROC_NULL){
        if(neighbours_ranks[LEFT]= MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send1[j] = field.g(inner_index_cols,j);
            }
            MPI_Sendrecv(&send1[0], send1.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,
                         &rcv1[0], rcv1.size(), MPI_DOUBLE, neighbours_ranks[RIGHT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.g(inner_index_cols+1,j) = rcv1[j];
            }
        }
    }
    if(neighbours_ranks[LEFT]!= MPI_PROC_NULL){
        if(neighbours_ranks[RIGHT]= MPI_PROC_NULL){
            for(int j=0; j<matrix.num_rows(); j++){
                send1[j] = field.g(1,j);
            }
            MPI_Sendrecv(&send1[0], send1.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,
                         &rcv1[0], rcv1.size(), MPI_DOUBLE, neighbours_ranks[LEFT], 0,  MPI_COMMUNICATOR, &status);

            for(int j=0; j<matrix.num_rows(); j++){
                field.g(0,j) = rcv1[j];
            }
        }
    }
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


//double reduce_min(){
    // // Finding the global Minimum
// // ...
// double localMin = *std::min_element(myVector.begin(),myVector.end());
// double globalMin;
// MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//}

double Communication::reduce_min(double value){
    double global_min ; 
    MPI_Allreduce(&value, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMMUNICATOR);
    return global_min;
}


double Communication::reduce_sum(double residual){
    double globalsum ; 
    MPI_Allreduce(&residual, &globalsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMMUNICATOR);
    return globalsum;
    }

