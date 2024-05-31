#include <algorithm>
#include <iostream>
#include <iomanip>

#include "Communication.hpp"
#include "Fields.hpp"

Fields::Fields(double nu, double dt, double tau, int size_x, int size_y, double UI, double VI, double PI, double alpha, double beta, double GX, double GY, double TI)
    : _nu(nu), _dt(dt), _tau(tau), _alpha(alpha),  _beta(beta), _gx(GX), _gy(GY) {
    _U = Matrix<double>(size_x + 2, size_y + 2, UI);
    _V = Matrix<double>(size_x + 2, size_y + 2, VI);
    _P = Matrix<double>(size_x + 2, size_y + 2, PI);
    _T = Matrix<double>(size_x + 2, size_y + 2, TI);
    _F = Matrix<double>(size_x + 2, size_y + 2, 0.0);
    _G = Matrix<double>(size_x + 2, size_y + 2, 0.0);
    _RS = Matrix<double>(size_x + 2, size_y + 2, 0.0);
}

void Fields::printMatrix(Grid &grid) {
    std::cout << std::fixed;
    std::cout << std::setprecision(1); // digits after decimal point

    std::cout << "P matrix" << std::endl;
    for (auto j = grid.size_y() + 1; j >= 0; j--) {
        for (auto i = 0; i <= grid.size_x() + 1; i++) {
            std::cout << _P(i, j) << ", ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    std::cout << "U matrix" << std::endl;
    for (auto j = grid.size_y() + 1; j >= 0; j--) {
        for (auto i = 0; i <= grid.size_x() + 1; i++) {
            std::cout << _U(i, j) << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "V matrix" << std::endl;
    for (auto j = grid.size_y() + 1; j >= 0; j--) {
        for (auto i = 0; i <= grid.size_x() + 1; i++) {
            std::cout << _V(i, j) << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "T matrix" << std::endl;
    for (auto j = grid.size_y() + 1; j >= 0; j--) {
        for (auto i = 0; i <= grid.size_x() + 1; i++) {
            std::cout << _T(i, j) << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << std::setprecision(4); // digits after decimal point
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

void Fields::calculate_fluxes(Grid &grid) {
    for (int i = 1; i <= grid.size_x() - 1; i++) {
        for (int j = 1; j <= grid.size_y(); j++) {
            _F(i, j) = _U(i, j) +
                       _dt * (_nu * (Discretization::laplacian(_U, i, j)) - Discretization::convection_u(_U, _V, i, j)) - _beta * _dt/2 * (_T(i,j)+_T(i+1,j)) * _gx;
        }
    }

    for (int i = 1; i <= grid.size_x(); i++) {
        for (int j = 1; j <= grid.size_y() - 1; j++) {
            _G(i, j) = _V(i, j) +
                       _dt * (_nu * (Discretization::laplacian(_V, i, j)) - Discretization::convection_v(_U, _V, i, j)) - _beta * _dt/2 * (_T(i,j)+_T(i,j+1)) * _gy;
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
    for (int i = 1; i <= grid.size_x() - 1; i++) {
        for (int j = 1; j <= grid.size_y(); j++) {
            _U(i, j) = _F(i, j) - _dt / grid.dx() * (_P(i + 1, j) - _P(i, j));
        }
    }

    for (int i = 1; i <= grid.size_x(); i++) {
        for (int j = 1; j <= grid.size_y() - 1; j++) {
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

void Fields::calculate_dt(Grid &grid) {
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

    _dt = _tau * _dt;
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }
double &Fields::T(int i, int j) { return _T(i, j); }


Matrix<double> &Fields::p_matrix() { return _P; }
Matrix<double> &Fields::u_matrix() { return _U; }
Matrix<double> &Fields::v_matrix() { return _V; }
Matrix<double> &Fields::f_matrix() { return _F; }
Matrix<double> &Fields::g_matrix() { return _G; }
Matrix<double> &Fields::rs_matrix() { return _RS; }
Matrix<double> &Fields::t_matrix() { return _T; }




double Fields::dt() const { return _dt; }
