#include <mpi.h>
#include <iostream>

static void init_parallel(int iproc, int jproc){
    MPI_Init(nullptr, nullptr);

    int num_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    if (num_proc != iproc * jproc) {
        std::cerr << "Incompatible number of processors and domain decomposition!";
        return;
    }

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    std::cout << "Hello world from rank " << my_rank
              << "of" << num_proc << std::endl;
    MPI_Finalize();
}