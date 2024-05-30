#pragma once

#include <mpi.h>
#include <iostream>

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
        * @brief communicate a field
        *
        * @param[in] field 
        *
        */ 
        static void communicate();

        /**
        * @brief find minimum value across all processes
        *
        * @param[in] value
        *
        */ 
        static double reduce_min();

        /**
        * @brief find total sum across all processes
        *
        * @param[in] value
        *
        */ 
        static double reduce_sum();

        ~Communication() = default;
        
};
