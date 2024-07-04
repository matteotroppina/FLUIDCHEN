#include <algorithm>
#include <iostream>
#include <iomanip>

#include "Communication.hpp"
#include "Fields.hpp"

Fields::Fields(double nu, double dt, double tau, int size_x, int size_y, double length_x, double length_y, double UI,
               double VI, double PI, double alpha, double beta, double GX, double GY, double TI, double KI, double EI)
    : _nu(nu), _dt(dt), _tau(tau), _alpha(alpha),  _beta(beta), _gx(GX), _gy(GY), _length_x(length_x), _length_y(length_y) {
    _U  = Matrix<double>(size_x + 2, size_y + 2, UI);
    _V  = Matrix<double>(size_x + 2, size_y + 2, VI);
    _P  = Matrix<double>(size_x + 2, size_y + 2, PI);
    _T  = Matrix<double>(size_x + 2, size_y + 2, TI);
    _F  = Matrix<double>(size_x + 2, size_y + 2, 0.0);
    _G  = Matrix<double>(size_x + 2, size_y + 2, 0.0);
    _RS = Matrix<double>(size_x + 2, size_y + 2, 0.0);

    //turbulence model
    double nuTI = _C_nu * std::pow(KI,2)/EI; //initial turbulent viscosity
    _K     = Matrix<double>(size_x + 2, size_y + 2, KI);
    _E     = Matrix<double>(size_x + 2, size_y + 2, EI);
    _nuT   = Matrix<double>(size_x + 2, size_y + 2, nuTI);

    // Low-Reynolds formulation matrices
    _ReT    = Matrix<double>(size_x + 2, size_y + 2, 0.0);
    _damp2  = Matrix<double>(size_x + 2, size_y + 2, 1.0);
    _dampmu = Matrix<double>(size_x + 2, size_y + 2, 1.0);
    _L_k = Matrix<double>(size_x + 2, size_y + 2, 0.0);
    _L_e = Matrix<double>(size_x + 2, size_y + 2, 0.0);
    _yplus = Matrix<double>(size_x + 2, size_y + 2, 0.0);
    _dist_y = Matrix<double>(size_x + 2, size_y + 2, 1e10);
    _dist_x = Matrix<double>(size_x + 2, size_y + 2, 1e10);
}



void Fields::calculate_fluxes(Grid &grid, bool turbulence_started) {

    for (int i = 1; i <= grid.itermax_x() - 1; i++) {
        for (int j = 1; j <= grid.size_y(); j++) {

            double nu_Tx;
            if(turbulence_started){ nu_Tx = Discretization::interpolate(_nuT, i, j, 1, 0); }
            else{ nu_Tx = 0.0; }

            _F(i, j) =  _U(i, j) +
                        _dt * ((_nu + nu_Tx) * (Discretization::laplacian(_U, i, j)) - Discretization::convection_u(_U, _V, i, j)); // - _beta * _dt/2 * (_T(i,j)+_T(i+1,j)) * _gx;
        }
    }

    for (int i = 1; i <= grid.size_x(); i++) {
        for (int j = 1; j <= grid.itermax_y() - 1; j++) {

            double nu_Ty;
            if(turbulence_started){ nu_Ty = Discretization::interpolate(_nuT, i, j, 0, 1); }
            else{ nu_Ty = 0.0; }

            _G(i, j) =  _V(i, j) +
                        _dt * ((_nu + nu_Ty)* (Discretization::laplacian(_V, i, j)) - Discretization::convection_v(_U, _V, i, j)); //- _beta * _dt/2 * (_T(i,j)+_T(i,j+1)) * _gy;

        }
    }
}

void Fields::calculate_rs(Grid &grid) {

    for (int i = 1; i <= grid.size_x(); i++) {
        for (int j = 1; j <= grid.size_y(); j++) {
            _RS(i, j) = 1 / _dt * ((_F(i, j) - _F(i - 1, j)) / grid.dx() + (_G(i, j) - _G(i, j - 1)) / grid.dy());
        }
    }
}

#ifdef __CUDACC__
void Fields::calculate_velocities(Grid &grid) {

    int itermax_x = grid.itermax_x();
    int itermax_y = grid.itermax_y();
    int size_y = grid.size_y();
    int size_x = grid.size_x();
    double dx = grid.dx();
    double dy = grid.dy();
    double dt = this->dt();

    double * u_matrix = _U.raw_pointer();
    double * v_matrix = _V.raw_pointer();
    double * p_matrix = _P.raw_pointer();
    double * F_matrix = _F.raw_pointer();
    double * G_matrix = _G.raw_pointer();

    double * d_u_matrix;
    cudaMalloc(&d_u_matrix, sizeof(double) * (size_x + 2) * (size_y + 2));
    cudaMemcpy(d_u_matrix, u_matrix, sizeof(double) * (size_x + 2) * (size_y + 2), cudaMemcpyHostToDevice);
    double * d_v_matrix;
    cudaMalloc(&d_v_matrix, sizeof(double) * (size_x + 2) * (size_y + 2));
    cudaMemcpy(d_v_matrix, v_matrix, sizeof(double) * (size_x + 2) * (size_y + 2), cudaMemcpyHostToDevice);
    double * d_p_matrix;
    cudaMalloc(&d_p_matrix, sizeof(double) * (size_x + 2) * (size_y + 2));
    cudaMemcpy(d_p_matrix, p_matrix, sizeof(double) * (size_x + 2) * (size_y + 2), cudaMemcpyHostToDevice);
    double * d_F_matrix;
    cudaMalloc(&d_F_matrix, sizeof(double) * (size_x + 2) * (size_y + 2));
    cudaMemcpy(d_F_matrix, F_matrix, sizeof(double) * (size_x + 2) * (size_y + 2), cudaMemcpyHostToDevice);
    double * d_G_matrix;
    cudaMalloc(&d_G_matrix, sizeof(double) * (size_x + 2) * (size_y + 2));
    cudaMemcpy(d_G_matrix, G_matrix, sizeof(double) * (size_x + 2) * (size_y + 2), cudaMemcpyHostToDevice);


    #pragma acc parallel loop collapse(2) deviceptr(d_u_matrix, d_p_matrix, d_F_matrix)
    for (int i = 1; i <= itermax_x - 1; i++) {
        for (int j = 1; j <= size_y; j++) {
            int idx = i + j * (size_x + 2);
            int idx_right = idx + 1;
            d_u_matrix[idx] = d_F_matrix[idx] - dt / dx * (d_p_matrix[idx_right] - d_p_matrix[idx]);
        }
    }

    #pragma acc parallel loop collapse(2) deviceptr(d_v_matrix, d_p_matrix, d_G_matrix)
    for (int i = 1; i <= size_x; i++) {
        for (int j = 1; j <= itermax_y - 1; j++) {
            int idx = i + j * (size_x + 2);
            int idx_top = idx + (size_x + 2);
            d_v_matrix[idx] = d_G_matrix[idx] - dt / dy * (d_p_matrix[idx_top] - d_p_matrix[idx]);
        }
    }

    cudaMemcpy(u_matrix, d_u_matrix, sizeof(double) * (size_x + 2) * (size_y + 2), cudaMemcpyDeviceToHost);
    cudaMemcpy(v_matrix, d_v_matrix, sizeof(double) * (size_x + 2) * (size_y + 2), cudaMemcpyDeviceToHost);
    cudaFree(d_u_matrix);
    cudaFree(d_v_matrix);
    cudaFree(d_p_matrix);
    cudaFree(d_F_matrix);
    cudaFree(d_G_matrix);
}
#else
void Fields::calculate_velocities(Grid &grid) {

    for (int i = 1; i <= grid.itermax_x() - 1; i++) {
        for (int j = 1; j <= grid.size_y(); j++) {
            _U(i, j) = _F(i, j) - _dt / grid.dx() * (_P(i + 1, j) - _P(i, j));
        }
    }

    for (int i = 1; i <= grid.size_x(); i++) {
        for (int j = 1; j <= grid.itermax_y() - 1; j++) {
            _V(i, j) = _G(i, j) - _dt / grid.dy() * (_P(i, j + 1) - _P(i, j));
        }
    }
}
#endif

void Fields::calculate_temperature(Grid &grid) {
    for (int i = 1; i <= grid.size_x(); i++) {
        for (int j = 1; j <= grid.size_y(); j++){
            if (grid.cell(i,j).type() == cell_type::FLUID){
                _T(i, j) = _T(i,j) + _dt * ((_alpha * Discretization::laplacian(_T,i,j)) - Discretization::convection_t(_T,_U,_V,i,j));
            }
        }
    }
}

void Fields::calculate_walldist(Grid &grid) {
    const double PI = 3.14159265358979311599796346854;
    double dx = grid.dx();
    double dy = grid.dy();

    for (auto fluid_cell : grid.fluid_cells()) {
        int i = fluid_cell->i();
        int j = fluid_cell->j();
        double least_dist_x = 1e10;
        double least_dist_y = 1e10;
        double least_dist = 1e10;
        double angle = PI / 4;
        for (auto wall_cell : grid.fixed_wall_cells()) {
            int i_wall = wall_cell->i();
            int j_wall = wall_cell->j();

            double x_component = (i - i_wall) * dx;
            double y_component = (j - j_wall) * dy;
            double distance = std::sqrt(x_component * x_component + y_component * y_component);
            if (distance < least_dist) {
                least_dist = distance;
                angle = std::atan2(y_component, x_component);
                least_dist_x = std::abs(least_dist * std::cos(angle));
                least_dist_y = std::abs(least_dist * std::sin(angle));
            }
        }

        if (least_dist_x > grid.dx()*0.95){ // floating point error
            _dist_x(i, j) = least_dist_x;
        }
        if (least_dist_y > grid.dy()*0.95){
            _dist_y(i, j) = least_dist_y;
        }

        double abs_angle = std::abs(angle);
        // angle is less than 45 degrees
        if ((abs_angle < PI / 4) || (abs_angle > 3 * PI / 4)) { // wall is along x-axis
            _dist_x(i, j) = least_dist_x;
            _dist_y(i, j) = 1e10;
        } else if ((abs_angle > PI / 4 && abs_angle < PI / 2) ||
                   (abs_angle > PI / 2 && abs_angle < 3 * PI / 4)) { // wall is along y-axis
            _dist_y(i, j) = least_dist_y;
            _dist_x(i, j) = 1e10;
        }
    }
}

void Fields::calculate_yplus(Grid &grid) {
    //skin friction coefficient
    double C_f = 0.058 * pow(_nu, 0.2);
    for (int i = 1; i <= grid.size_x(); i++) {
        for (int j = 1; j <= grid.size_y(); j++) {
            if (grid.cell(i,j).type() == cell_type::FLUID){
                double u_tau = C_f * std::pow(_K(i,j), 0.5); // friction velocity
                double y; // distance from wall
                if (_dist_x(i,j) < _dist_y(i,j)){
                    y = _dist_x(i,j);
                } else {
                    y = _dist_y(i,j);
                }
                _yplus(i,j) = y * u_tau / _nu;
            }
        }
    }
}

void Fields::calculate_damping(Grid &grid){
        for (int i = 1; i <= grid.size_x(); i++) {
        for (int j = 1; j <= grid.size_y(); j++){
            if (grid.cell(i,j).type() == cell_type::FLUID){

                // Low-Reynolds formulation up to the viscous sublayer
                if (_yplus(i, j) < 30) {
                    // Launder-Sharma coefficients
                    ReT(i, j) = yplus(i, j) * std::pow(K(i, j), 2) / (_nu * E(i, j));
                    damp2(i, j) = 1 - 0.3 * std::exp(-std::pow(ReT(i, j), 2));
                    dampmu(i, j) = std::exp(-3.40 / std::pow((1 + 0.02*ReT(i, j)), 2));
                    L_e(i, j) = 2 * _nu * nuT(i,j) * std::pow((Discretization::laplacian(u_matrix(), i, j) + Discretization::laplacian(v_matrix(), i, j)), 2);

                    //dissipation rate of K at the wall
                    if (_dist_x(i, j) < _dist_y(i, j)){
                        //divergence of sqrt(k) at the wall
                        L_k(i, j) = 2*_nu * std::pow(std::pow(_K(i,j), 0.5) / _dist_x(i,j), 2); // assuming k = 0 at wall
                    } else {
                        L_k(i, j) = 2*_nu * std::pow(std::pow(_K(i,j), 0.5) / _dist_y(i, j), 2); // assuming k = 0 at wall
                    }
                } else {
                    damp2(i, j) = 1;
                    dampmu(i, j) = 1;
                    L_k(i, j) = 0;
                    L_e(i, j) = 0;
                }

            }
        }
    }
}

void Fields::calculate_dt(Grid &grid, bool turbulence_started) {
    double dx_2 = grid.dx() * grid.dx();
    double dy_2 = grid.dy() * grid.dy();

    double u_max = _U.max_abs_value();
    double v_max = _V.max_abs_value();

    double coefficient = (dx_2 * dy_2) / (dx_2 + dy_2);
    double conv_cond = coefficient / (2 * _nu);
    double cfl_x = grid.dx() / u_max;
    double cfl_y = grid.dy() / v_max;
    double new_cond = (1 / (2 * _alpha)) * (1/ (1/ dx_2 + 1/dy_2));

    _dt = std::min({conv_cond, cfl_x, cfl_y, new_cond});

//    if (turbulence_started){
//        double k_max = _K.max_abs_value();
//        double eps_max = _E.max_abs_value();
//        double k_cond = 1 / (2 * k_max * (1 / dx_2 + 1 / dy_2));
//        double eps_cond = 1 / (2 * eps_max * (1 / dx_2 + 1 / dy_2));
//        _dt = std::min(_dt, 25 * k_cond);
//        _dt = std::min(_dt, 25 * eps_cond);
//    }

    _dt = _tau * _dt;

    _dt = Communication::reduce_min(_dt);
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }
double &Fields::T(int i, int j) { return _T(i, j); }
double &Fields::K(int i, int j) {return _K(i,j);}
double &Fields::E(int i, int j) {return _E(i,j);}
double &Fields::nuT(int i, int j) {return _nuT(i,j);}
double &Fields::nu(){return _nu;}
double &Fields::yplus(int i, int j) {return _yplus(i,j);}
double &Fields::dist_y(int i, int j) {return _dist_y(i,j);}
double &Fields::dist_x(int i, int j) {return _dist_x(i,j);}
double &Fields::ReT(int i, int j) {return _ReT(i,j);}
double &Fields::damp2(int i, int j) {return _damp2(i,j);}
double &Fields::dampmu(int i, int j) {return _dampmu(i,j);}
double &Fields::L_k(int i, int j) {return _L_k(i,j);}
double &Fields::L_e(int i, int j) {return _L_e(i,j);}


Matrix<double> &Fields::p_matrix() { return _P; }
Matrix<double> &Fields::u_matrix() { return _U; }
Matrix<double> &Fields::v_matrix() { return _V; }
Matrix<double> &Fields::f_matrix() { return _F; }
Matrix<double> &Fields::g_matrix() { return _G; }
Matrix<double> &Fields::rs_matrix() { return _RS; }
Matrix<double> &Fields::t_matrix() { return _T; }
Matrix<double> &Fields::k_matrix() { return _K; }
Matrix<double> &Fields::e_matrix() { return _E; }
Matrix<double> &Fields::nuT_matrix() { return _nuT; }
Matrix<double> &Fields::yplus_matrix() { return _yplus; }
Matrix<double> &Fields::dist_y_matrix() { return _dist_y; }
Matrix<double> &Fields::dist_x_matrix() { return _dist_x; }
Matrix<double> &Fields::ReT_matrix() { return _ReT; }
Matrix<double> &Fields::damp2_matrix() { return _damp2; }
Matrix<double> &Fields::dampmu_matrix() { return _dampmu; }
Matrix<double> &Fields::L_k_matrix() { return _L_k; }
Matrix<double> &Fields::L_e_matrix() { return _L_e; }

double Fields::dt() const { return _dt; }
double Fields::length_x() const {return _length_x;}
double Fields::length_y() const {return _length_y;}
