#pragma once

#include <utility>

#include "Boundary.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
#include <assert.h>
#include <math.h>


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
    virtual void solve(Fields &field, Grid &grid) = 0;
};

/**
 * @brief K-Epsilon model for turbulence viscosity
 *
 */
class K_EPS_model : public ViscositySolver {
  public:
    K_EPS_model() = default;

    // /**
    //  * @brief Constructor of K-Epsilon model
    //  *
    //  * @param[in] 
    //  */
    // K_EPS_model();

    virtual ~K_EPS_model() = default;

    /**
     * @brief Solve the viscosity equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual void solve(Fields &field, Grid &grid);

  private:

    // nu : viscosity, nuT: turbulent eddy viscosity, nuT = C0 * (k^2 / epsilon)

    // http://www.dicat.unige.it/guerrero/turbulence2021/slides/lecture6/6closure_models_RANS_part7_1.pdf
    // and
    // Numerical Investigation of Turbulent Flows Using k-epsilon By Reza Barati

    // Launder-Sharma model constants
    double _C_nu{0.09}; // in paper reffered as C_nu
    double _C1{1.44};
    double _C2{1.92};
    double _damp1{1.0};
    double _sk{1.0};  // in paper reffered as sigma k
    double _se{1.3};  // in paper reffered as sigma e

};
