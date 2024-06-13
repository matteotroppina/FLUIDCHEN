#include <cmath>
#include "Communication.hpp"
#include "ViscositySolver.hpp"

// nu : viscosity, nuT: turbulent eddy viscosity, nuT = C0 * (k^2 / epsilon)
double C0 = 0.09; // in paper reffered as C_nu
double C1 = 1.44;
double C2 = 1.92;
double sk = 1.0; // in paper reffered as sigma k
double se = 1.3; // in paper reffered as sigma e
