#include <cmath>
#include "Communication.hpp"
#include "PressureSolver.hpp"

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

double gpu_solve(double *p_matrix, double *p_matrix_old, const double *rs_matrix, const bool *fluid_mask,
                 const uint8_t *boundary_type, const gridParams grid, const int num_iterations) {

    double dx = grid.dx;
    double dy = grid.dy;
    int imax = grid.imax;
    int jmax = grid.jmax;
    int size_linear = (imax + 2) * (jmax + 2);
    double coeff = 1.0 / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));

    //present means data present on GPU
    // JACOBI ITERATION

    for (int iter = 0; iter < num_iterations; iter++) {
        #pragma acc parallel loop collapse(2) present(p_matrix[0 : size_linear], p_matrix_old[0 : size_linear], rs_matrix[0 : size_linear], fluid_mask[0 : size_linear], boundary_type[0 : size_linear])
        for (int i = 1; i <= imax; i++) {
            for (int j = 1; j <= jmax; j++) {
                int idx = i + j * (imax + 2);
                int idx_left = idx - 1;
                int idx_right = idx + 1;
                int idx_top = idx + (imax + 2);
                int idx_bottom = idx - (imax + 2);

                p_matrix[idx] =
                    coeff * ((p_matrix_old[idx_left] + p_matrix_old[idx_right]) / (dx * dx) +
                             (p_matrix_old[idx_top] + p_matrix_old[idx_bottom]) / (dy * dy) - rs_matrix[idx]);
                p_matrix[idx] *= fluid_mask[idx]; // set to zero if not fluid cell
            }
        }

        // swap p_matrix and p_matrix_old pointers
        double *temp = p_matrix;
        p_matrix = p_matrix_old;
        p_matrix_old = temp;

        // Update GPU with new pointer values
        #pragma acc update device(p_matrix, p_matrix_old)
    }


    double res = 0.0;
    int size_fluid_cells = 0;

    // CALCULATE RESIDUAL
    #pragma acc parallel loop collapse(2) reduction(+:res,size_fluid_cells) present(p_matrix[0:size_linear], p_matrix_old[0:size_linear], rs_matrix[0:size_linear], fluid_mask[0:size_linear], boundary_type[0:size_linear])
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            int idx = i + j * (imax + 2);
            int idx_left = idx - 1;
            int idx_right = idx + 1;
            int idx_top = idx + (imax + 2);
            int idx_bottom = idx - (imax + 2);

            double val = (p_matrix[idx_left] - 2.0 * p_matrix[idx] + p_matrix[idx_right]) / (dx * dx) +
                         (p_matrix[idx_bottom] - 2.0 * p_matrix[idx] + p_matrix[idx_top]) / (dy * dy) -
                         rs_matrix[idx];
            res += val * val * fluid_mask[idx];
            size_fluid_cells += fluid_mask[idx];
        }
    }

    res = std::sqrt(res / size_fluid_cells);

    //TODO apply boundary conditions


    return res;

}

