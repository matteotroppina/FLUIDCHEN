#ifndef CFDLAB_UTILSGPU_H
#define CFDLAB_UTILSGPU_H

#include <cstdint>

#ifdef __CUDACC__
#warning "CUDA enabled"
#include <cuda_runtime.h>

typedef struct gridParams {
    int imax;
    int jmax;
    double dx;
    double dy;
    size_t size_fluid_cells;
} gridParams;

class Grid;

void generate_gpu_masks(Grid &_grid, bool* &fluid_mask, uint8_t* &boundary_type, uint8_t* &border_position, bool * &d_fluid_mask, uint8_t * &d_boundary_type, uint8_t * &d_border_position);

void allocate_gpu_memory(Grid &_grid, double * &d_p_matrix_new, double * &d_p_matrix, double * &d_rs_matrix);

void free_gpu_memory(double * &d_p_matrix_new, double * &d_p_matrix, double * &d_rs_matrix, bool * &d_fluid_mask, uint8_t * &d_boundary_type, uint8_t * &d_border_position);

#endif
#endif // CFDLAB_UTILSGPU_H