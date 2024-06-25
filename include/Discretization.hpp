#pragma once

#include "Datastructures.hpp"

/**
 * @brief Static discretization methods to modify the fields
 *
 */
class Discretization {
  public:
    Discretization() = default;

    /**
     * @brief Constructor to set the discretization parameters
     *
     * @param[in] cell size in x direction
     * @param[in] cell size in y direction
     * @param[in] upwinding coefficient
     */
    Discretization(double dx, double dy, double gamma);

    /**
     * @brief Convection in x direction using donor-cell scheme
     *
     * @param[in] x-velocity field
     * @param[in] y-velocity field
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j);

    /**
     * @brief Convection in y direction using donor-cell scheme
     *
     * @param[in] x-velocity field
     * @param[in] y-velocity field
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j);

    /**
     * @brief Laplacian term discretization using central difference
     *
     * @param[in] data to be discretized
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double laplacian(const Matrix<double> &A, int i, int j);

    /**
     * @brief Terms of laplacian needed for SOR, i.e. excluding unknown value at
     * (i,j)
     *
     * @param[in] data to be discretized
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double sor_helper(const Matrix<double> &P, int i, int j);

    /**
     * @brief Compute interpolated value in the middle between two grid points via linear interpolation.
     *
     * @param[in] A data to be interpolated
     * @param[in] i index of first value used for interpolation
     * @param[in] j index of first value used for interpolation
     * @param[in] i_offset defines index of the second value used for interpolation as i+i_offset
     * @param[in] j_offset defines index of the second value used for interpolation as j+j_offset
     * @param[out] result
     *
     */
    static double interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset);

    /**
     * @brief Convection term for temperature
     *
     * @param[in] temperature field
     * @param[in] x-velocity field
     * @param[in] y-velocity field
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double convection_t(const Matrix<double> &T, const Matrix<double> &U, const Matrix<double> &V, int i, int j);

    /**
     * @brief Convection term for turbulent kinetic energy and dissipation rate
     *
     * @param[in] turbulent kinetic energy or dissipation rate field
     * @param[in] x-velocity field
     * @param[in] y-velocity field
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double convection_KEPS(const Matrix<double> &K, const Matrix<double> &U, const Matrix<double> &V, int i, int j);

    static double Discretization::laplacian_KEPS(const Matrix<double> &K, const Matrix<double> &nuT, const double nu, const double _sk, int i,int j);

    /**
     * @brief Strain rate tensor
     *
     * @param[in] x-velocity field
     * @param[in] y-velocity field
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double strain_rate(const Matrix<double> &U, const Matrix<double> &V, int i, int j);

  private:
    static double _dx;
    static double _dy;
    static double _gamma; //donor cells
};
