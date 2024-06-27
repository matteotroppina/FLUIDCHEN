#include "Discretization.hpp"

#include <cmath>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
}

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double du2_dx = 1 / _dx * (std::pow(interpolate(U, i, j, 1, 0), 2) - std::pow(interpolate(U, i, j, -1, 0), 2)) +
                    _gamma / _dx *
                        ((std::abs(interpolate(U, i, j, 1, 0)) * (U(i, j) - U(i + 1, j)) / 2) -
                         std::abs(interpolate(U, i, j, -1, 0)) * (U(i - 1, j) - U(i, j)) / 2);
    double dv2_dy = 1 / _dy *
                        (((interpolate(V, i, j, 1, 0)) * (interpolate(U, i, j, 0, 1))) -
                         ((interpolate(V, i, j - 1, 1, 0)) * (interpolate(U, i, j, 0, -1)))) +
                    _gamma / _dy *
                        ((std::abs(interpolate(V, i, j, 1, 0)) * (U(i, j) - U(i, j + 1)) / 2) -
                         (std::abs(interpolate(V, i, j - 1, 1, 0)) * (U(i, j - 1) - U(i, j)) / 2));

    return du2_dx + dv2_dy;;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double dv2_dy = 1/ _dy * (std::pow(interpolate(V,i,j,0,1), 2) - std::pow(interpolate(V,i,j - 1,0,1), 2)) +
                    _gamma/_dy * ((std::abs(interpolate(V,i,j,0,1)) * (V(i, j) - V(i, j + 1)) / 2) -
                                    std::abs(interpolate(V, i, j - 1, 0, 1)) * (V(i, j - 1)- V(i, j)) / 2) ;
    double du2_dx = 1/ _dx * (((interpolate(U, i, j, 0, 1)) * (interpolate(V, i, j, 1, 0)))  - ( (interpolate(U,i-1,j,0,1)) * (interpolate(V,i-1,j,1,0)))) +
                    _gamma / _dx * ( (std::abs(interpolate(U,i,j,0,1))*(V(i, j) - V(i + 1, j)) / 2) -
                                    (std::abs(interpolate(U,i - 1,j,0,1)) * (V(i - 1, j) - V(i, j)) / 2 ));

    return dv2_dy + du2_dx;;
}

double Discretization::convection_t(const Matrix<double> &T, const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double duT_dx = 1/_dx * (U(i,j)* interpolate(T,i,j,1,0) - U(i-1,j)* interpolate(T,i,j,-1,0))
                    + _gamma/_dx * (std::abs(U(i,j))* (T(i,j)-T(i+1,j))/2 - std::abs(U(i-1,j))*(T(i-1,j)-T(i,j))/2);
    double dvT_dy = 1/_dy * (V(i,j) * interpolate(T,i,j,0,1) - V(i,j-1)* interpolate(T,i,j,0,-1))
                    + _gamma/_dy * (std::abs(V(i,j))* (T(i,j)-T(i,j+1))/2 - std::abs(V(i,j-1))*(T(i,j-1)-T(i,j))/2);
    return duT_dx + dvT_dy;
}

double Discretization::laplacian(const Matrix<double> &A, int i, int j) {
    return (A(i+1,j) - 2*A(i,j) + A(i-1,j))/std::pow(_dx,2) + (A(i,j+1) - 2*A(i,j) + A(i,j-1))/std::pow(_dy,2);
}


double Discretization::mixed_derivative(const Matrix<double> &A, int i, int j) {
    return (A(i+1,j+1) - A(i-1,j+1) - A(i+1,j-1) + A(i-1,j-1))/(4*_dx*_dy);
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {
    return (A(i, j) + A(i+i_offset, j+j_offset))/2;
}

double Discretization::convection_KEPS(const Matrix<double> &K, const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double duK_dx = 1/_dx * (U(i,j)* interpolate(K,i,j,1,0) - U(i-1,j)* interpolate(K,i-1,j,1,0))
                    + _gamma/_dx * (std::abs(U(i,j))* (K(i,j)-K(i+1,j))/2 - std::abs(U(i-1,j))*(K(i-1,j)-K(i,j))/2.0);
    double dvK_dy = 1/_dy * (V(i,j) * interpolate(K,i,j,0,1) - V(i,j-1)* interpolate(K,i,j,0,-1))
                    + _gamma/_dy * (std::abs(V(i,j))* (K(i,j)-K(i,j+1))/2 - std::abs(V(i,j-1))*(K(i,j-1)-K(i,j))/2.0);
    return duK_dx + dvK_dy;
}

double Discretization::laplacian_KEPS(const Matrix<double> &K, const Matrix<double> &nuT_i, const Matrix<double> &nuT_j, const double nu, const double _sk,int i,int j){

    double laplacian_x = ( (K(i+1,j) - K(i,j)) * (nu + nuT_i(i,j)/_sk) - (nu + nuT_i(i-1, j)/_sk) * (K(i,j)+ K(i-1,j)) ) / std::pow(_dx, 2);
    double laplacian_y = ( (K(i,j+1) - K(i,j)) * (nu + nuT_j(i,j)/_sk) - (nu + nuT_j(i, j-1)/_sk) * (K(i,j)+ K(i,j-1)) ) / std::pow(_dy, 2); 
    return laplacian_x + laplacian_y;
}


double Discretization::strain_rate(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // double S_xx = std::pow((U(i,j) - U(i-1,j))/std::pow(_dx,2),2);
    // double S_yy = std::pow((V(i,j) - V(i,j-1))/std::pow(_dy,2),2);
    // double su = (interpolate(U,i,j+1,-1,0) - interpolate(U,i,j-1,-1,0)) / (2 * _dy);
    // double sv = (interpolate(V,i+1,j,0,-1) - interpolate(V,i-1,j,0,-1)) / (2 * _dx);;
    // double S_xy = std::pow((su + sv),2);

    // return (2 * S_xx + S_xy + 2 * S_yy);

    double S_xx = std::pow((U(i,j) - U(i-1,j))/std::pow(_dx,2),2);
    double S_yy = std::pow((V(i,j) - V(i,j-1))/std::pow(_dy,2),2);
    double su = (U(i, j + 1) + U(i - 1, j + 1) - U(i, j - 1) - U(i - 1, j - 1)) * (0.25 / _dy);
    double sv = (V(i + 1, j) + V(i + 1, j - 1) - V(i - 1, j) - V(i - 1, j - 1)) * (0.25 / _dx);
    double Sxy = std::pow((su+sv), 2);

    return (2 * S_xx + Sxy + 2 * S_yy);
}

double Discretization::laplacian_x(const Matrix<double> &A, int i, int j) {
    return (A(i+1,j) - 2*A(i,j) + A(i-1,j))/std::pow(_dx,2);
}

double Discretization::laplacian_y(const Matrix<double> &A, int i, int j) {
    return (A(i,j+1) - 2*A(i,j) + A(i,j-1))/std::pow(_dy,2);
}