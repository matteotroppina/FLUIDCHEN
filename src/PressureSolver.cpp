#include <cmath>
#include "Communication.hpp"
#include "PressureSolver.hpp"
#include "UtilsGPU.h"

SOR::SOR(double omega) : _omega(omega) {}

double SOR::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {

    double dx = grid.dx();
    double dy = grid.dy();

    double coeff = _omega / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) = (1.0 - _omega) * field.p(i, j) +
                        coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
    }

    double res = 0.0;
    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += (val * val);
    }
    {
        res = rloc / (grid.fluid_cells().size());
        res = std::sqrt(res);
    }

    for (auto &b : boundaries) {
        b->applyPressure(field);
    }

    return res;
}

__global__ void jacobiKernel(double *d_p_matrix_new, double *d_p_matrix, const double *d_rs_matrix, const bool *d_fluid_mask,
                             const double coeff, const int imax, const int jmax, const double dx, const double dy) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int idx = i + j * (imax + 2);
    int idx_left = idx - 1;
    int idx_right = idx + 1;
    int idx_top = idx + (imax + 2);
    int idx_bottom = idx - (imax + 2);

    if (i <= imax && j <= jmax && d_fluid_mask[idx]) {

        d_p_matrix_new[idx] =
            coeff * ((d_p_matrix[idx_left] + d_p_matrix[idx_right]) / (dx * dx) +
                     (d_p_matrix[idx_top] + d_p_matrix[idx_bottom]) / (dy * dy) - d_rs_matrix[idx]);
    }
}


////call with <<<1, 1>>> to swap pointers
//__global__ void swapPointers(double *d_p_matrix, double *d_p_matrix_new) {
//    __syncthreads();
//    double *temp = d_p_matrix;
//    d_p_matrix = d_p_matrix_new;
//    d_p_matrix_new = temp;
//}

#pragma acc routine seq
void swapPointers(double *&d_p_matrix, double *&d_p_matrix_new) {
    double *temp = d_p_matrix;
    d_p_matrix = d_p_matrix_new;
    d_p_matrix_new = temp;
}

double gpu_psolve(double *d_p_matrix, double *d_p_matrix_new, const double *d_rs_matrix, const bool *d_fluid_mask,
                  const uint8_t *d_boundary_type, const gridParams grid, const int num_iterations) {

    double dx = grid.dx;
    double dy = grid.dy;
    int imax = grid.imax;
    int jmax = grid.jmax;
    int size_linear = (imax + 2) * (jmax + 2);
    double coeff = 1.0 / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
    double res = 0;
    int size_fluid_cells = grid.size_fluid_cells;

    //present means data present on GPU
    // JACOBI ITERATION
    #pragma acc loop seq
    for (int iter = 0; iter < num_iterations; iter++) {
        #pragma acc parallel loop collapse(2) vector_length(128) \
        present(d_rs_matrix[0 : size_linear], d_fluid_mask[0 : size_linear], d_boundary_type[0 : size_linear]) \
        deviceptr(d_p_matrix_new, d_p_matrix)
        for (int i = 1; i <= imax; i++) {
            for (int j = 1; j <= jmax; j++) {
                int idx = i + j * (imax + 2);
                int idx_left = idx - 1;
                int idx_right = idx + 1;
                int idx_top = idx + (imax + 2);
                int idx_bottom = idx - (imax + 2);

                d_p_matrix_new[idx] =
                    coeff * ((d_p_matrix[idx_left] + d_p_matrix[idx_right]) / (dx * dx) +
                             (d_p_matrix[idx_top] + d_p_matrix[idx_bottom]) / (dy * dy) - d_rs_matrix[idx]);
                d_p_matrix_new[idx] *= d_fluid_mask[idx]; // set to zero if not fluid cell
            }
        }

        swapPointers(d_p_matrix, d_p_matrix_new);

    }

    //TODO apply boundary conditions

    // CALCULATE RESIDUAL
    #pragma acc parallel loop collapse(2) reduction(+:res) present(d_p_matrix[0:size_linear], d_rs_matrix[0:size_linear], d_fluid_mask[0:size_linear])
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            int idx = i + j * (imax + 2);
            int idx_left = idx - 1;
            int idx_right = idx + 1;
            int idx_top = idx + (imax + 2);
            int idx_bottom = idx - (imax + 2);

            double val = (d_p_matrix[idx_left] - 2.0 * d_p_matrix[idx] + d_p_matrix[idx_right]) / (dx * dx) +
                         (d_p_matrix[idx_bottom] - 2.0 * d_p_matrix[idx] + d_p_matrix[idx_top]) / (dy * dy) -
                         d_rs_matrix[idx];
            res += val * val * d_fluid_mask[idx];
        }
    }

    return res / size_fluid_cells;

}