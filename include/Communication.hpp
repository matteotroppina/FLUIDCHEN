#pragma once

#include <mpi.h>
#include <iostream>
#include "Fields.hpp"
#include "Datastructures.hpp"
#include <array>

// stores the rank of the current process in the custom communicator
inline int my_rank_global;
// stores the 2D coordinates of the current process in the cartesian grid
inline int my_coords_global[2];

extern MPI_Comm MPI_COMMUNICATOR;
class Communication{
    public:

        Communication() = default;
        
        /**
        * @brief Communication constructor
        *
        * @param[in] argn number of arguments from command line
        * @param[in] args arguments from command line
        *
        */ 
        Communication(int argn, char **args);

        /**
        * @brief initialize communication
        *
        * @param[in] argn number of arguments from command line
        * @param[in] args arguments from command line
        *
        */ 
        static void init_parallel(int argn, char **args);

        /**
        * @brief finalize communication
        *
        */ 
        static void finalize();

        /**
                 * @brief get the rank of the current process
        *
        */
        static int get_rank();

        /**
        * @brief get the number of processes
        *
        */
        static int get_size();

        /**
        * @brief get the cartesian communicator
        *
        */
        static MPI_Comm get_communicator();

        /**
        * @brief get the 2D coordinates of the current process
        *
        */
        static std::array<int, 2> get_coords();

        /**
        * @brief get the neighbours of the current process
        *
        */
        static std::array<int, 4> get_neighbours();

        /**
        * @brief communicate a matrix
        *
        * @param[in] matrix
        *
        */ 
        static void communicate(Matrix<double> &matrix);

        /**
        * @brief find minimum value across all processes
        *
        * @param[in] value
        *
        */ 
        static double reduce_min(double value);

        /**
        * @brief find total sum across all processes
        *
        * @param[in] value
        *
        */ 
        static double reduce_sum(double residual);


        ~Communication() = default;
        
};
