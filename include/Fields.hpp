#pragma once

#include "Datastructures.hpp"
#include "Discretization.hpp"
#include "Grid.hpp"
#include <assert.h>
#include <math.h>

/**
 * @brief Class of container and modifier for the physical fields
 *
 */
class Fields {
  public:
    Fields() = default;

    /**
     * @brief Constructor for the fields
     *
     * @param[in] kinematic viscosity
     * @param[in] initial timestep size
     * @param[in] adaptive timestep coefficient
     * @param[in] number of cells in x direction
     * @param[in] number of cells in y direction
     * @param[in] initial x-velocity
     * @param[in] initial y-velocity
     * @param[in] initial pressure
     *
     */
    Fields(double _nu, double _dt, double _tau, int size_x, int size_y, double UI, double VI, double PI, double alpha, double beta, double GX, double GY, double TI, double KI, double EI);

    void printMatrix(Grid &grid);
    void printCellTypes(Grid &grid);
    void printBorders(Grid &grid);


    /**
     * @brief Calculates the convective and diffusive fluxes in x and y
     * direction based on explicit discretization of the momentum equations
     *
     * @param[in] grid in which the fluxes are calculated
     *
     */
    void calculate_fluxes(Grid &grid, bool turbulence_started);

    /**
     * @brief Right hand side calculations using the fluxes for the pressure
     * Poisson equation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_rs(Grid &grid);

    /**
     * @brief Velocity calculation using pressure values
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_velocities(Grid &grid);

     /**
     * @brief temperature calculation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_temperature(Grid &grid);

    /**
     * @brief Adaptive step size calculation using x-velocity condition,
     * y-velocity condition and CFL condition
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_dt(Grid &grid, bool turbulence_started);

    /**
     * @brief turbulent kinetic energy and dissipation rate calculation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_nuT(Grid &grid, const double &C0);

    void calculate_yplus(Grid &grid);
    void calculate_delta_y(Grid &grid);

    /**
     * @brief damping factors calculation for Low-Reynolds formulation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_damping(Grid &grid);

    /// x-velocity index based access and modify
    double &u(int i, int j);

    /// y-velocity index based access and modify
    double &v(int i, int j);

    /// pressure index based access and modify
    double &p(int i, int j);

    /// RHS index based access and modify
    double &rs(int i, int j);

    /// x-momentum flux index based access and modify
    double &f(int i, int j);

    /// y-momentum flux index based access and modify
    double &g(int i, int j);


    /// temperature index based access and modify
    double &T(int i, int j);

    // turbulent kinetic energy index based access and modify
    double &K(int i, int j);

    /// dissipation rate index based access and modify
    double &E(int i, int j);

    /// turbulent viscosity index based access and modify
    double &nuT(int i, int j);

    /// turbulent viscosity index based access and modify
    double &nuT_i(int i, int j);

    /// turbulent viscosity index based access and modify
    double &nuT_j(int i, int j);

    double &yplus(int i, int j);
    double &delta_y(int i, int j);

    /// turbulent reynolds index based access and modify
    double &ReT(int i, int j);

    /// damping factor 1 index based access and modify
    double &damp1(int i, int j);

    /// damping factor 2 index based access and modify
    double &damp2(int i, int j);

    /// damping factor mu index based access and modify
    double &dampmu(int i, int j);

    /// get timestep size
    double dt() const;

    /// pressure matrix access and modify
    Matrix<double> &p_matrix();

    /// velocity u matrix access and modify
    Matrix<double> &u_matrix();

    /// velocity v matrix access and modify
    Matrix<double> &v_matrix();

    /// f matrix access and modify
    Matrix<double> &f_matrix();

    /// g pressure matrix access and modify
    Matrix<double> &g_matrix();

    /// RHS matrix access and modify
    Matrix<double> &rs_matrix();

    /// temperature matrix access and modify
    Matrix<double> &t_matrix();

    /// turbulent kinetic energy matrix access and modify
    Matrix<double> &k_matrix();

    /// dissipation rate matrix access and modify 
    Matrix<double> &e_matrix();

    /// turbulent viscosity matrix access and modify
    Matrix<double> &nuT_matrix();

    /// turbulent viscosity i matrix access and modify
    Matrix<double> &nuT_i_matrix();

    /// turbulent viscosity j matrix access and modify
    Matrix<double> &nuT_j_matrix();

    Matrix<double> &yplus_matrix();

    Matrix<double> &delta_y_matrix();

    Matrix<double> &ReT_matrix();

    Matrix<double> &damp1_matrix();

    Matrix<double> &damp2_matrix();

    Matrix<double> &dampmu_matrix();

    /// get kinematic viscosity
    double &nu();

    /// get diameter for pipe case
    double &size_y();

  private:
    /// x-velocity matrix
    Matrix<double> _U;
    /// y-velocity matrix
    Matrix<double> _V;
    /// pressure matrix
    Matrix<double> _P;
    /// x-momentum flux matrix
    Matrix<double> _F;
    /// y-momentum flux matrix
    Matrix<double> _G;
    /// right hand side matrix
    Matrix<double> _RS;
    /// temperature matrix
    Matrix<double> _T;
    /// turbulent kinetic energy matrix
    Matrix<double> _K;
    /// dissipation rate matrix
    Matrix<double> _E;
    /// turbulent viscosity matrix
    Matrix<double> _nuT;
    /// turbulent viscosity matrix i
    Matrix<double> _nuT_i;
    /// turbulent viscosity matrix j
    Matrix<double> _nuT_j;
    /// turbulent reynolds
    Matrix<double> _ReT;

    Matrix<double> _damp1;
    Matrix<double> _damp2;
    Matrix<double> _dampmu;

    Matrix<double> _yplus;
    Matrix<double> _delta_y;

    /// kinematic viscosity
    double _nu;
    /// timestep size
    double _dt;
    /// adaptive timestep coefficient
    double _tau;
    /// thermal diffusivity
    double _alpha;
    double _beta;
    double _size_y;

    /// gravitional acceleration in x direction
    double _gx{0.0};
    /// gravitional acceleration in y direction
    double _gy{0.0};

    double _Cmu{0.09};
};
