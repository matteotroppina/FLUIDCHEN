#pragma once

#include <utility>

#include "Boundary.hpp"
#include "Fields.hpp"
#include "Grid.hpp"

/**
 * @brief Abstract class for pressure Poisson equation solver
 *
 */
class PressureSolver {
  public:
    PressureSolver() = default;
    virtual ~PressureSolver() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) = 0;
};

/**
 * @brief Successive Over-Relaxation algorithm for solution of pressure Poisson
 * equation
 *
 */
class SOR : public PressureSolver {
  public:
    SOR() = default;

    /**
     * @brief Constructor of SOR solver
     *
     * @param[in] relaxation factor
     */
    SOR(double omega);

    virtual ~SOR() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    double _omega;
};

// Jacobi for GPU
class JAC_GPU : public PressureSolver {
  public:
    JAC_GPU() = default;

    virtual ~JAC_GPU() = default;

    //TODO : Implement an enum for different type of boundary conditions
    virtual double solve(std::vector<double> p_matrix, const std::vector<std::tuple<int, int>> &fluid_cell_indices,
                         const std::vector<std::tuple<int, int>> &boundary_indices, const std::vector<int> &boundary_type);
};
