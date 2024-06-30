#include "Boundary.hpp"

#include <cmath>

Boundary::Boundary(std::vector<Cell *> cells) : _cells(cells) {}
void Boundary::applyFlux(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        // B_NW cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
            field.f(i - 1, j) = field.u(i - 1, j);
            field.g(i, j) = field.v(i, j);
            continue;
        }
        // B_SE cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
            field.f(i, j) = field.u(i, j);
            field.g(i, j - 1) = field.v(i, j - 1);
            continue;
        }
        // B_NE cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
            field.f(i, j) = field.u(i, j);
            field.g(i, j) = field.v(i, j);
            continue;
        }
        // B_SW cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
            field.f(i, j) = field.u(i - 1, j);
            field.g(i, j - 1) = field.v(i, j - 1);
            continue;
        }

        if (cell->is_border(border_position::RIGHT)) {
            field.f(i, j) = field.u(i, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.f(i - 1, j) = field.u(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.g(i, j) = field.v(i, j);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.g(i, j - 1) = field.v(i, j - 1);
        }
    }
}

void Boundary::applyVelocity(Fields &field){ (void)field; }
void Boundary::applyPressure(Fields &field){ (void)field; }
void Boundary::applyTemperature(Fields &field) { (void)field; }

void Boundary::applyTurbulence(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        // B_NW cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
            field.K(i,j) = ( field.K(i-1,j) + field.K(i,j+1) ) / 2.0;
            field.E(i,j) = ( field.E(i-1,j) + field.E(i,j+1) ) / 2.0;
            
            field.nuT_i(i, j) = C0 * std::pow(field.K(i, j), 2) / field.E(i, j);
            continue;
        }
        // B_SE cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
            field.K(i,j) = ( field.K(i+1,j) + field.K(i,j+1) ) / 2.0;
            field.E(i,j) = ( field.E(i+1,j) + field.E(i,j+1) ) / 2.0;
            
            field.nuT_i(i, j) = C0 * std::pow(field.K(i, j), 2) / field.E(i, j);
            continue;
        }
        // B_NE cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
            field.K(i,j) = ( field.K(i+1,j) + field.K(i,j+1) ) / 2.0;
            field.E(i,j) = ( field.E(i+1,j) + field.E(i,j+1) ) / 2.0;
            
            field.nuT_i(i, j) = C0 * std::pow(field.K(i, j), 2) / field.E(i, j);
            continue;
        }
        // B_SW cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
            field.K(i,j) = ( field.K(i-1,j) + field.K(i,j-1) ) / 2.0;
            field.E(i,j) = ( field.E(i-1,j) + field.E(i,j-1) ) / 2.0;
            
            field.nuT_i(i, j) = C0 * std::pow(field.K(i, j), 2) / field.E(i, j);
            continue;
        }

        if (cell->is_border(border_position::RIGHT)) {
            field.K(i,j) = field.K(i+1,j);
            field.E(i,j) = field.E(i+1,j);
            
            field.nuT_i(i, j) = C0 * std::pow(field.K(i, j), 2) / field.E(i, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.K(i,j) = field.K(i-1,j);
            field.E(i,j) = field.E(i-1,j);
            
            field.nuT_i(i, j) = C0 * std::pow(field.K(i, j), 2) / field.E(i, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.K(i,j) = field.K(i,j+1);
            field.E(i,j) = field.E(i,j+1);
            
            field.nuT_i(i, j) = C0 * std::pow(field.K(i, j), 2) / field.E(i, j);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.K(i,j) = field.K(i,j-1);
            field.E(i,j) = field.E(i,j-1);
            
            field.nuT_i(i, j) = C0 * std::pow(field.K(i, j), 2) / field.E(i, j);
        }
    }
}


InnerObstacle::InnerObstacle(std::vector<Cell *> cells) : Boundary(cells) {}

void InnerObstacle::applyVelocity(Fields &field) {
    for (auto cell: _cells) {
        int i = cell->i();
        int j = cell->j();

        field.u(i,j) = 0;
        field.v(i,j) = 0;
    }
} // do nothing

// void InnerObstacle::applyFlux(Fields &field) {}
// void InnerObstacle::applyPressure(Fields &field) {}


FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, double wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::applyVelocity(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        // B_NW cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
            field.v(i, j) = 0;
            field.u(i, j) = -field.u(i, j+1);
            field.v(i, j-1) = - field.v(i-1, j-1);

            continue;
        }
        // B_SE cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j-1) = 0;
            field.u(i-1, j) = -field.u(i-1, j-1);
            field.v(i, j) = -field.v(i+1, j);

            continue;
        }
        // B_NE cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j) = 0;
            field.u(i-1, j) = -field.u(i-1, j+1);
            field.v(i, j-1) = - field.v(i+1, j-1);

            continue;
        }
        // B_SW cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
            field.u(i-1, j) = 0;
            field.v(i, j-1) = 0;
            field.u(i, j) = -field.u(i, j-1);
            field.v(i,j) = -field.v(i-1, j);

            continue;
        }

        // B_E cell
        if (cell->is_border(border_position::RIGHT)) {
            field.v(i, j) = -field.v(i + 1, j);
            field.u(i, j) = 0;
        }
        // B_W cell
        if (cell->is_border(border_position::LEFT)) {
            field.v(i, j) = -field.v(i - 1, j);
            field.u(i - 1, j) = 0;
        }
        // B_N cell
        if (cell->is_border(border_position::TOP)) {
            field.u(i, j) = -field.u(i, j + 1);
            field.v(i, j) = 0;
        }
        // B_S cell
        if (cell->is_border(border_position::BOTTOM)) {
            field.u(i, j) = -field.u(i, j - 1);
            field.v(i, j - 1) = 0;
        }

        // forbidden cells with two opposite borders or three boundaries
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::TOP)) {
            std::cerr << "there are forbidden cells with two opposite borders or three boundaries \n";
            std::cerr << "i: " << i << " j: " << j << std::endl;
            exit(1);
        }
        if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::RIGHT)) {
            std::cerr << "there are forbidden cells with two opposite borders or three boundaries \n";
            std::cerr << "i: " << i << " j: " << j << std::endl;
            exit(1);
        }
        // forbidden cells with obstacles only consisting of 1 cell
        if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::RIGHT) &&
            cell->is_border(border_position::TOP) && cell->is_border(border_position::BOTTOM)) {
            std::cerr << "there are forbidden cells with four boundaries \n";
            std::cerr << "i: " << i << " j: " << j << std::endl;
            exit(1);
        }
    }
}

void FixedWallBoundary::applyPressure(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j); // i = 0
        }

        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j); // i = imax + 1
        }

        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1); // j = 0
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1); // j = jmax + 1
        }

        // B_NW cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
            field.p(i, j) = (field.p(i, j + 1) + field.p(i - 1, j)) / 2;
        }
        // B_SE cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = (field.p(i + 1, j) + field.p(i, j - 1)) / 2;
        }
        // B_SW cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
            field.p(i, j) = (field.p(i - 1, j) + field.p(i, j - 1)) / 2;
        }
        // B_NE cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = (field.p(i, j + 1) + field.p(i + 1, j)) / 2;
        }
    }
}

void FixedWallBoundary::applyTemperature(Fields &field) {

    if (_wall_temperature == -1) { // Neumann
        for (auto cell : _cells) {
            int i = cell->i();
            int j = cell->j();

            if (cell->is_border(border_position::RIGHT)) {
                field.T(i, j) = field.T(i + 1, j);
            }

            if (cell->is_border(border_position::LEFT)) {
                field.T(i, j) = field.T(i - 1, j);
            }

            if (cell->is_border(border_position::TOP)) {
                field.T(i, j) = field.T(i, j + 1);
            }

            if (cell->is_border(border_position::BOTTOM)) {
                field.T(i, j) = field.T(i, j - 1);
            }
        }
    } else { // Dirichlet

        for (auto cell : _cells) {
            int i = cell->i();
            int j = cell->j();

            if (cell->is_border(border_position::BOTTOM)) {
                field.T(i, j) = 2 * _wall_temperature - field.T(i, j - 1);
            }
            if (cell->is_border(border_position::TOP)) {
                field.T(i, j) = 2 * _wall_temperature - field.T(i, j + 1);
            }

            if (cell->is_border(border_position::RIGHT)) {
                field.T(i, j) = 2 * _wall_temperature - field.T(i + 1, j);
            }

            if (cell->is_border(border_position::LEFT)) {
                field.T(i, j) = 2 * _wall_temperature - field.T(i - 1, j);
            }
        }
    }
}


// INFLOW
FixedVelocityBoundary::FixedVelocityBoundary(std::vector<Cell *> cells, double inflow_u_velocity,
                                             double inflow_v_velocity)
    : Boundary(cells) {
    _inflow_u_velocity.insert(std::pair<int, double>(GeometryIDs::fixed_velocity, inflow_u_velocity));
    _inflow_v_velocity.insert(std::pair<int, double>(GeometryIDs::fixed_velocity, inflow_v_velocity));
}

FixedVelocityBoundary::FixedVelocityBoundary(std::vector<Cell *> cells, std::map<int, double> inflow_u_velocity,
                                             std::map<int, double> inflow_v_velocity,
                                             std::map<int, double> wall_temperature)
    : Boundary(cells), _inflow_u_velocity(inflow_u_velocity), _inflow_v_velocity(inflow_v_velocity),
      _wall_temperature(wall_temperature) {}

void FixedVelocityBoundary::applyVelocity(Fields &field) {

    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::BOTTOM)) {
            field.u(i, j) = 2 * _inflow_u_velocity[GeometryIDs::fixed_velocity] - field.u(i, j - 1);
            field.v(i, j - 1) = _inflow_v_velocity[GeometryIDs::fixed_velocity];
        }
        if (cell->is_border(border_position::TOP)) {
            field.u(i, j) = 2 * _inflow_u_velocity[GeometryIDs::fixed_velocity] - field.u(i, j + 1);
            field.v(i, j) = _inflow_v_velocity[GeometryIDs::fixed_velocity];
        }

        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = _inflow_u_velocity[GeometryIDs::fixed_velocity];
            field.v(i, j) = 2 * _inflow_v_velocity[GeometryIDs::fixed_velocity] - field.v(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = _inflow_u_velocity[GeometryIDs::fixed_velocity];
            field.v(i, j) = 2 * _inflow_v_velocity[GeometryIDs::fixed_velocity] - field.v(i - 1, j);
        }
    }
}

void FixedVelocityBoundary::applyPressure(Fields &field) {
    // Neumann condition
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1);
        }
    }
}

// do we have to prescribe something for k and eps?
void FixedVelocityBoundary::applyTurbulence(Fields &field) {

    double d_pipe = 2.0; //TO DO: get pipe diameter (length y, physical)
    double l = 0.07 * d_pipe;
    double I = 0.1;

    for (auto cell : _cells) {

        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::BOTTOM)) {
            double u = (field.u(i, j - 1) + field.u(i, j)) / 2;
            double v = field.v(i, j - 1);
            double mean_uv = std::sqrt(std::pow(u,2) + std::pow(v,2));
            double k_boundary = 3/2 * std::pow(mean_uv*I,2);
            double eps_boundary = std::pow(C0,0.75) * std::pow(k_boundary, 1.5) / l; 

            field.K(i, j) = 2*k_boundary - field.K(i, j - 1);
            field.E(i, j) = 2*eps_boundary - field.E(i, j-1);
        }

        if (cell->is_border(border_position::TOP)) {
            double u = (field.u(i, j) + field.u(i, j + 1)) / 2;
            double v = field.v(i, j);
            double mean_uv = std::sqrt(std::pow(u,2) + std::pow(v,2));
            double k_boundary = 3/2 * std::pow(mean_uv*I,2);
            double eps_boundary = std::pow(C0,0.75) * std::pow(k_boundary, 1.5) / l; 

            field.K(i, j) = 2*k_boundary - field.K(i, j + 1);
            field.E(i, j) = 2*eps_boundary - field.E(i, j+1);
            
        }

        if (cell->is_border(border_position::RIGHT)) {
            double u = field.u(i, j); // = _inflow_u_velocity[GeometryIDs::fixed_velocity]
            double v = (field.v(i, j) + field.v(i + 1, j)) / 2;
            double mean_uv = std::sqrt(std::pow(u,2) + std::pow(v,2));
            double k_boundary = 3/2 * std::pow(mean_uv*I,2);
            double eps_boundary = std::pow(C0,0.75) * std::pow(k_boundary, 1.5) / l; 

            field.K(i, j) = 2*k_boundary - field.K(i + 1, j);
            field.E(i, j) = 2*eps_boundary - field.E(i + 1, j);
            field.nuT(i,j) = C0 * std::pow(field.K(i, j), 2) / field.E(i, j);

        }

        if (cell->is_border(border_position::LEFT)) {
            double u = field.u(i - 1, j);
            double v = (field.v(i, j) + field.v(i - 1, j)) / 2;
            double mean_uv = std::sqrt(std::pow(u,2) + std::pow(v,2));
            double k_boundary = 3/2 * std::pow(mean_uv*I,2);
            double eps_boundary = std::pow(C0,0.75) * std::pow(k_boundary, 1.5) / l; 

            field.K(i, j) = 2*k_boundary - field.K(i - 1, j);
            field.E(i, j) = 2*eps_boundary - field.E(i - 1, j);
        }
    }
}

//OUTFLOW
ZeroGradientBoundary::ZeroGradientBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

ZeroGradientBoundary::ZeroGradientBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void ZeroGradientBoundary::applyVelocity(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = field.u(i + 1, j);
            field.v(i, j) = field.v(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.u(i, j) = field.u(i - 2, j);
            field.u(i - 1, j) = field.u(i, j);
            field.v(i, j) = field.v(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.v(i, j) = field.v(i, j + 1);
            field.u(i, j) = field.u(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.v(i, j) = field.v(i, j - 2);
            field.v(i, j -1) = field.v(i, j);
            field.u(i, j) = field.u(i, j -1);
        }
    }
}

void ZeroGradientBoundary::applyPressure(Fields &field) {
    // Dirichlet condition pressure on boundary = 0
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = -field.p(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = -field.p(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = -field.p(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = -field.p(i, j - 1);
        }
    }
}

/*
    TODO --> 
    do we want to make the code more general and handle inflow/outflow for all possible combinations? (change the code above)

    INFLOW
    K = c_bc * |u|^2   (euclidean norm)    c_bc [0.003, 0.01]
    eps = Cµ * k^(3/2) / l0   

    OUTFLOW
    grad(k) = 0
    grad(eps) = 0
    do we have to change the boundaries conditions for velocity   n*[grad(u) + grad(u)^T] = 0
*/
void ZeroGradientBoundary::applyTurbulence(Fields &field) {
        for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.K(i, j) = field.K(i + 1, j);
            field.E(i, j) = field.E(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.K(i, j) = field.K(i - 1, j);
            field.E(i, j) = field.E(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.K(i, j) = field.K(i, j + 1);
            field.E(i, j) = field.E(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.K(i, j) = field.K(i, j - 1);
            field.E(i, j) = field.E(i, j - 1);
        }
    }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : Boundary(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::applyVelocity(Fields &field) {

    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::BOTTOM)) {
            field.u(i, j) = 2 * _wall_velocity[GeometryIDs::moving_wall] - field.u(i, j - 1);
            field.v(i, j - 1) = 0;
        }
        if (cell->is_border(border_position::TOP)) {
            field.u(i, j) = 2 * _wall_velocity[GeometryIDs::moving_wall] - field.u(i, j + 1);
            field.v(i, j) = 0;
        }

        if (cell->is_border(border_position::RIGHT)) {
            field.v(i, j) = 2 * _wall_velocity[GeometryIDs::moving_wall] - field.v(i + 1, j);
            field.u(i, j) = 0;
        }

        if (cell->is_border(border_position::LEFT)) {
            field.v(i, j) = 2 * _wall_velocity[GeometryIDs::moving_wall] - field.v(i - 1, j);
            field.u(i - 1, j) = 0;

        }
    }
}

void MovingWallBoundary::applyPressure(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1);
        }
    }
}
