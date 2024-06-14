#pragma once

#include <utility>

#include "Boundary.hpp"
#include "Fields.hpp"
#include "Grid.hpp"


class ViscositySolver {
  public:
    ViscositySolver() = default;
    virtual ~ViscositySolver() = default;

    /**
     * @brief Solve viscosity equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) = 0;
};

/**
 * @brief K-Epsilon model for turbulence viscosity
 *
 */
class K_EPS_model : public ViscositySolver {
  public:
    K_EPS_model() = default;

    /**
     * @brief Constructor of K-Epsilon model
     *
     * @param[in] 
     */
    K_EPS_model();

    virtual ~K_EPS_model() = default;

    /**
     * @brief Solve the viscosity equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:

    // nu : viscosity, nuT: turbulent eddy viscosity, nuT = C0 * (k^2 / epsilon)
    double _C0{0.09}; // in paper reffered as C_nu
    double _C1{1.44};
    double _C2{1.92};
    double _sk{1.0}; // in paper reffered as sigma k
    double _se{1.3}; // in paper reffered as sigma e
    double _decouple{0.0}; //decouple turbulent equations

};
