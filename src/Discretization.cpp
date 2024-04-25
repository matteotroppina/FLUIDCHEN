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
    return 1/_dx * (pow(interpolate(U, i, j, 1, 0), 2) - pow(interpolate(U, i, j, -1, 0), 2))
           + _gamma/_dx * (abs(interpolate(U, i, j, 1, 0))*(U(i,j) - U(i+1, j))/2
                             - abs(interpolate(U, i, j, -1, 0))*(U(i-1, j) - U(i,j))/2);
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {}

double Discretization::laplacian(const Matrix<double> &A, int i, int j) {}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {
    return (A(i, j) + A(i+i_offset, j+j_offset))/2;
}