//
// Created by dan on 23.06.24.
//

#ifndef CFDLAB_UTILSGPU_H
#define CFDLAB_UTILSGPU_H

#ifdef __CUDACC__
#warning "CUDA enabled"
#include <cuda_runtime.h>


typedef struct gridParams {
    int imax;
    int jmax;
    double dx;
    double dy;
    unsigned long size_fluid_cells;
} gridParams;

#endif
#endif // CFDLAB_UTILSGPU_H


