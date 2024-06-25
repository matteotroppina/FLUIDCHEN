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
        double k2 = Discretization::laplacian_KEPS(field.k_matrix(),field.nuT_matrix(), nu, _sk, i ,j);
        double k3 = (nu + field.nuT(i,j)) * Discretization::strain_rate(field.u_matrix(), field.v_matrix(), i, j);
        double k4 = field.E(i, j);
        
        double e1 = Discretization::convection_KEPS(field.e_matrix(),field.u_matrix(),field.v_matrix(),i,j);
        double e2 = Discretization::laplacian_KEPS(field.k_matrix(),field.nuT_matrix(), nu, _se, i ,j);
        double e3 = (_C1 * field.E(i, j)/field.K(i,j) ) * k3;
        double e4 = _C2 * std::pow(field.E(i, j),2)/field.K(i,j);

        field.K(i, j) = field.K(i, j) + dt * (-k1 + k2 + k3 - k4);
        field.E(i, j) = field.E(i, j) + dt * (-e1 + e2 + e3 - e4);

        // field.nuT(i,j) = _C0 * field.K(i,j) * field.K(i,j) / field.E(i,j);

        //std::cout << "working" << i << "," << j << std::endl;

        // TODO --> implement stopping criteria
        //      --> update gamma


        /* TODO --> 
         1. choose default mixing length l0
         2  What is nu_0? is it the constant viscosity?
         3. how do we choose t*?
         4. we have to understand how to start the model. At the beginning we simulate with constant k and eps initialized with the inital values
         then at t* the turbulence model is activated and we update k and eps as above

         INITIAL CONDITIONS
         at t < t*
         k0 = (nu0 / l0)^2
         eps_0 = Cµ * k0^(3/2) / l0
        */
    }
}