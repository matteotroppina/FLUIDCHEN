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

void Fields::printMatrix(Grid &grid) {
    
    std::cout << std::fixed;
    std::cout << std::setprecision(7); // digits after decimal point

    std::cout << "K matrix" << std::endl;
    for (auto j = grid.size_y() + 1; j >= 0; j--) {
        for (auto i = 0; i <= grid.size_x() + 1; i++) {
            std::cout << _K(i, j) << ", ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    std::cout << "E matrix" << std::endl;
    for (auto j = grid.size_y() + 1; j >= 0; j--) {
        for (auto i = 0; i <= grid.size_x() + 1; i++) {
            std::cout << _E(i, j) << ", ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    // std::cout << "P matrix" << std::endl;
    // for (auto j = grid.size_y() + 1; j >= 0; j--) {
    //     for (auto i = 0; i <= grid.size_x() + 1; i++) {
    //         std::cout << _P(i, j) << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << std::endl;

    // std::cout << "U matrix" << std::endl;
    // for (auto j = grid.size_y() + 1; j >= 0; j--) {
    //     for (auto i = 0; i <= grid.size_x() + 1; i++) {
    //         std::cout << _U(i, j) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout << "V matrix" << std::endl;
    // for (auto j = grid.size_y() + 1; j >= 0; j--) {
    //     for (auto i = 0; i <= grid.size_x() + 1; i++) {
    //         std::cout << _V(i, j) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout << "T matrix" << std::endl;
    // for (auto j = grid.size_y() + 1; j >= 0; j--) {
    //     for (auto i = 0; i <= grid.size_x() + 1; i++) {
    //         std::cout << _T(i, j) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout << std::setprecision(4); // digits after decimal point
    
}

void Fields::printCellTypes(Grid &grid){
    

    std::map<cell_type, char> cellTypeToChar = {
        {cell_type::FLUID, GeometryIDs::fluid},
        {cell_type::FIXED_WALL, GeometryIDs::fixed_wall},
        {cell_type::FIXED_VELOCITY, GeometryIDs::fixed_velocity},
        {cell_type::ZERO_GRADIENT, GeometryIDs::zero_gradient},
        {cell_type::MOVING_WALL, GeometryIDs::moving_wall},
        {cell_type::INNER_OBSTACLE, ' ' - '0'}};


    std::cout << std::fixed;
    std::cout << std::setprecision(1);

    std::cout <<std::endl << "Cell types" << std::endl;
    for (auto j = grid.size_y() + 1; j >= 0; j--) {
        for (auto i = 0; i <= grid.size_x() + 1; i++) {
            char cell_id = cellTypeToChar[grid.cell(i, j).type()];
            cell_id += '0';
            std::cout << cell_id  << ", ";
        }
        std::cout << std::endl;
    }

    
}

void Fields::printBorders(Grid &grid) {
    

    std::cout << std::endl << "Borders" << std::endl;

    for (auto j = grid.size_y() + 1; j >= 0; j--) {
        for (auto i = 0; i <= grid.size_x() + 1; i++) {
            if (grid.cell(i, j).is_border(border_position::LEFT)) {
                std::cout << "L, ";
            } else if (grid.cell(i, j).is_border(border_position::RIGHT)) {
                    std::cout << "R, ";
            } else if (grid.cell(i, j).is_border(border_position::TOP)) {
                    std::cout << "T, ";
            } else if (grid.cell(i, j).is_border(border_position::BOTTOM)) {
                    std::cout << "B, ";
            } else {
                std::cout << ".  ";
            }
        }
        std::cout << std::endl;
    }

    
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
                        L_k(i, j) = 2*_nu * std::pow((std::pow(_K(i,j), 0.5) - std::pow(_K(i-1,j), 0.5)) / grid.dx(), 2);
                    } else {
                        L_k(i, j) = 2*_nu * std::pow((std::pow(_K(i,j), 0.5) - std::pow(_K(i,j-1), 0.5)) / grid.dy(), 2);
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

    if (turbulence_started){
        double k_max = _K.max_abs_value();
        double eps_max = _E.max_abs_value();
        double k_cond = 1 / (2 * k_max * (1 / dx_2 + 1 / dy_2));
        double eps_cond = 1 / (2 * eps_max * (1 / dx_2 + 1 / dy_2));
        _dt = std::min(_dt, k_cond);
        _dt = std::min(_dt, eps_cond);
    }

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
