#include <cmath>

#include "Communication.hpp"
#include "ViscositySolver.hpp"


// Default constructor
// K_EPS_model::K_EPS_model() {}

double K_EPS_model::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    
    double dx = grid.dx();
    double dy = grid.dy();
    double dt = field.dt();

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        _decopule = field.E(i,j) / field.K(i,j);

        double k1 = (field.nuT(i,j) / 2) * Discretization::strain_rate(field.u_matrix(), field.v_matrix(), i, j);
        double k2 = field.nuT(i,j) / _sk *  Discretization::laplacian(field.k_matrix(), i, j);
        double k3 = Discretization::convection_KEPS(field.k_matrix(), field.u_matrix(), field.v_matrix(), i, j);
        double k4 = _decouple * field.K(i, j);
        
        double e1 = _C1 *_decouple * (field.nuT(i,j) / 2) * Discretization::strain_rate(field.u_matrix(), field.v_matrix(), i, j);
        double e2 = field.nuT(i,j) / _se * Discretization::laplacian(field.e_matrix(), i, j);
        double e3 = Discretization::convection_KEPS(field.e_matrix(), field.u_matrix(), field.v_matrix(), i, j);
        double e4 = _C2 * _decouple * field.E(i, j);

        field.K(i, j) = field.K(i, j) + dt * (k1 + k2 - (k3 + k4));
        field.E(i, j) = field.E(i, j) + dt * (e1 + e2 - (e3 + e4));

        field.nuT(i,j) = _C0 * field.K(i,j) * field.K(i,j) / field.E(i,j);


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