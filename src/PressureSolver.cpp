#include "PressureSolver.hpp"
#include "Communication.hpp"
#include <cmath>

SOR::SOR(double omega) : _omega(omega) {}

double SOR::solve(Fields &field, Grid &grid) {

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

    res = rloc / (grid.fluid_cells().size());

    return res;
}

#ifdef __CUDACC__
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
        d_p_matrix[idx] = d_p_matrix_new[idx];
    }
}

void apply_pressure_bcs(double *d_p_matrix, const bool *d_fluid_mask, const uint8_t *d_boundary_type, const uint8_t *d_border_position,
                const gridParams grid) {
    int imax = grid.imax;
    int jmax = grid.jmax;


    uint8_t ZERO_GRADIENT = static_cast<uint8_t>(cell_type::ZERO_GRADIENT); // zero gradient in velocity
    //uint8_t FIXED_VELOCITY = static_cast<uint8_t>(cell_type::FIXED_VELOCITY);
    //uint8_t FIXED_WALL = static_cast<uint8_t>(cell_type::FIXED_WALL);

    uint8_t TOP = static_cast<uint8_t>(border_position::TOP); // 0
    uint8_t BOTTOM = static_cast<uint8_t>(border_position::BOTTOM); // 1
    uint8_t LEFT = static_cast<uint8_t>(border_position::LEFT); // 2
    uint8_t RIGHT = static_cast<uint8_t>(border_position::RIGHT); // 3
    uint8_t TOP_LEFT = static_cast<uint8_t>(border_position::TOP) + 4; // 4
    uint8_t BOTTOM_RIGHT = static_cast<uint8_t>(border_position::BOTTOM) + 4; // 5
    uint8_t TOP_RIGHT = static_cast<uint8_t>(border_position::TOP) + 4; // 6
    uint8_t BOTTOM_LEFT = static_cast<uint8_t>(border_position::BOTTOM) + 4; // 7

    #pragma acc parallel loop collapse(2) deviceptr(d_p_matrix, d_fluid_mask, d_boundary_type, d_border_position) async
    for (int i = 0; i <= imax+1; i++) {
        for (int j = 0; j <= jmax+1; j++) {
            int idx = i + j * (imax + 2);
            int idx_left = idx - 1;
            int idx_right = idx + 1;
            int idx_top = idx + (imax + 2);
            int idx_bottom = idx - (imax + 2);

            if (!d_fluid_mask[idx]) {

                    double value_right = d_p_matrix[idx_right];
                    double value_left = d_p_matrix[idx_left];
                    double value_top = d_p_matrix[idx_top];
                    double value_bottom = d_p_matrix[idx_bottom];

                    uint8_t bpos = d_border_position[idx];

                    double zero_gradient_val = (bpos == TOP) ? -value_top :
                                               (bpos == BOTTOM) ? -value_bottom :
                                               (bpos == LEFT) ? -value_left :
                                               (bpos == RIGHT) ? -value_right : 0;

                    double fixed_wall_val = (bpos == TOP) ? value_top :
                                            (bpos == BOTTOM) ? value_bottom :
                                            (bpos == LEFT) ? value_left :
                                            (bpos == RIGHT) ? value_right :
                                            (bpos == TOP_LEFT) ? (value_top + value_left) / 2 :
                                            (bpos == BOTTOM_RIGHT) ? (value_bottom + value_right) / 2 :
                                            (bpos == TOP_RIGHT) ? (value_top + value_right) / 2 :
                                            (bpos == BOTTOM_LEFT) ? (value_bottom + value_left) / 2 : 0;

                    d_p_matrix[idx] = (d_boundary_type[idx] == ZERO_GRADIENT) ? zero_gradient_val : fixed_wall_val;
            }
        }
    }
}

double gpu_psolve(double *d_p_matrix, double *d_p_matrix_new, const double *d_rs_matrix, const bool *d_fluid_mask,
                  const uint8_t *d_boundary_type, const uint8_t *d_border_position, const gridParams grid,
                  const int num_iterations) {

    double dx = grid.dx;
    double dy = grid.dy;
    int imax = grid.imax;
    int jmax = grid.jmax;
    int size_linear = (imax + 2) * (jmax + 2);
    double coeff = 1.0 / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
    double res = 0;
    size_t size_fluid_cells = grid.size_fluid_cells;

    dim3 threadsPerBlock(8, 8);
    dim3 numBlocks((imax + 2 + threadsPerBlock.x - 1) / threadsPerBlock.x, (jmax + 2 + threadsPerBlock.y - 1) / threadsPerBlock.y);

    for (int iter = 0; iter < num_iterations; iter++) {
        jacobiKernel<<<numBlocks, threadsPerBlock>>>(d_p_matrix_new, d_p_matrix, d_rs_matrix, d_fluid_mask, coeff, imax, jmax, dx, dy);
    }

    apply_pressure_bcs(d_p_matrix, d_fluid_mask, d_boundary_type, d_border_position, grid);
    #pragma acc wait

    // CALCULATE RESIDUAL
    #pragma acc parallel loop collapse(2) reduction(+:res) deviceptr(d_p_matrix, d_rs_matrix, d_fluid_mask, d_rs_matrix)
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            int idx = i + j * (imax + 2);
            int idx_left = idx - 1;
            int idx_right = idx + 1;
            int idx_top = idx + (imax + 2);
            int idx_bottom = idx - (imax + 2);

            double val = (d_p_matrix[idx_left] - 2.0 * d_p_matrix[idx] + d_p_matrix[idx_right]) / (dx * dx) +
                         (d_p_matrix[idx_bottom] - 2.0 * d_p_matrix[idx] + d_p_matrix[idx_top]) / (dy * dy) - d_rs_matrix[idx];
            res += val * val * d_fluid_mask[idx];
        }
    }

    return res / size_fluid_cells;

}
#endif