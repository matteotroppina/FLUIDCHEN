#include <cmath>

#include "Communication.hpp"
#include "ViscositySolver.hpp"


// Default constructor
// K_EPS_model::K_EPS_model() {}

void K_EPS_model::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    
    double dx = grid.dx();
    double dy = grid.dy();
    double dt = field.dt();
    double nu = field.nu();

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        
        double k1 = Discretization::convection_KEPS(field.k_matrix(),field.u_matrix(),field.v_matrix(),i,j);
        double k2 = Discretization::laplacian_KEPS(field.k_matrix(), field.nuT_i_matrix(), field.nuT_j_matrix(), nu, _sk, i ,j);
        double k3 = (nu + field.nuT(i,j)) * Discretization::strain_rate(field.u_matrix(), field.v_matrix(), i, j);
        double k4 = field.E(i, j);
        
        double e1 = Discretization::convection_KEPS(field.e_matrix(),field.u_matrix(),field.v_matrix(),i,j);
        double e2 = Discretization::laplacian_KEPS(field.k_matrix(), field.nuT_i_matrix(), field.nuT_j_matrix(), nu, _se, i ,j);
        double e3 = _C1 * (field.E(i, j) / field.K(i,j)) * k3;
        double e4 = _C2 * std::pow(field.E(i, j),2)/field.K(i,j);

        field.K(i, j) = field.K(i, j) + dt * (-k1 + k2 + k3 - k4);
        field.E(i, j) = field.E(i, j) + dt * (-e1 + e2 + e3 - e4);

        // std::cout << "k1: " << k1 << std::endl;
        // std::cout << "k2: " << k2 << std::endl;
        // std::cout << "k3: " << k3 << std::endl;
        // std::cout << "k4: " << k4 << std::endl;
        // std::cout << "e1: " << e1 << std::endl;
        // std::cout << "e2: " << e2 << std::endl;
        // std::cout << "e3: " << e3 << std::endl;
        // std::cout << "e4: " << e4 << std::endl;
        // std::cout << "K: " << field.K(i, j) << std::endl;
        // std::cout << "E: " << field.E(i, j) << std::endl;

        // double R_T = std::pow(field.K(i,j),2) / (field.E(i,j) * nu);
        // _C0 = 0.09 * std::exp(-2.50/(1 + (R_T/50.0)));
        // _C2 = 1.92 * ( 1 - 0.3 * std::exp(-std::pow(R_T, 2)) );     //updates for next iteration?

    }
}