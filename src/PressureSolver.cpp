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

__global__ void jacobiKernel(double **d_p_matrix_new_dev, double **d_p_matrix_dev, const double *d_rs_matrix, const bool *d_fluid_mask,
                             const double coeff, const int imax, const int jmax, const double dx, const double dy) {

    double* d_p_matrix_new = *d_p_matrix_new_dev;
    double* d_p_matrix = *d_p_matrix_dev;
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


//call with <<<1, 1>>> to swap pointers
__global__ void swapPointers(double **d_p_matrix, double **d_p_matrix_new) {
    double *temp = *d_p_matrix;
    *d_p_matrix = *d_p_matrix_new;
    *d_p_matrix_new = temp;
}


double gpu_psolve(double *p_matrix, double *p_matrix_new, const double *rs_matrix, const bool *fluid_mask,
                  const uint8_t *boundary_type, const gridParams grid, const int num_iterations) {

    double dx = grid.dx;
    double dy = grid.dy;
    int imax = grid.imax;
    int jmax = grid.jmax;
    int size_linear = (imax + 2) * (jmax + 2);
    double coeff = 1.0 / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
    double res = 0;
    int size_fluid_cells = grid.size_fluid_cells;

    double **d_p_matrix, **d_p_matrix_new;
    cudaMalloc(&d_p_matrix, sizeof(double *));
    cudaMalloc(&d_p_matrix_new, sizeof(double *));

    cudaMemcpy(d_p_matrix, &p_matrix, sizeof(double *), cudaMemcpyHostToDevice);
    cudaMemcpy(d_p_matrix_new, &p_matrix_new, sizeof(double *), cudaMemcpyHostToDevice);

    dim3 threadsPerBlock(8, 8);
    dim3 numBlocks((imax + 2 + threadsPerBlock.x - 1) / threadsPerBlock.x, (jmax + 2 + threadsPerBlock.y - 1) / threadsPerBlock.y);

    for (int iter = 0; iter < num_iterations; iter++) {
        jacobiKernel<<<numBlocks, threadsPerBlock>>>(d_p_matrix_new, d_p_matrix, rs_matrix, fluid_mask, coeff, imax, jmax, dx, dy);
        swapPointers<<<1, 1>>>(d_p_matrix, d_p_matrix_new);
    }

    //TODO apply boundary conditions

    // CALCULATE RESIDUAL
    #pragma acc parallel loop collapse(2) reduction(+:res) present(p_matrix[0:size_linear], rs_matrix[0:size_linear], fluid_mask[0:size_linear])
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            int idx = i + j * (imax + 2);
            int idx_left = idx - 1;
            int idx_right = idx + 1;
            int idx_top = idx + (imax + 2);
            int idx_bottom = idx - (imax + 2);

            double val = (p_matrix[idx_left] - 2.0 * p_matrix[idx] + p_matrix[idx_right]) / (dx * dx) +
                         (p_matrix[idx_bottom] - 2.0 * p_matrix[idx] + p_matrix[idx_top]) / (dy * dy) - rs_matrix[idx];
            res += val * val * fluid_mask[idx];
        }
    }

    cudaFree(d_p_matrix);
    cudaFree(d_p_matrix_new);

    return res / size_fluid_cells;

}