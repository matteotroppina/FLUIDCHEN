#include "Boundary.hpp"

Boundary::Boundary(std::vector<Cell *> cells) : _cells(cells) {}

void Boundary::applyFlux(Fields &field) {


}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::applyVelocity(Fields &field) {
//    std::cout << "printing fixedwallboundary cells" << "\n";
//    for (auto cell : _cells){
//        std::cout << cell->i() << ", " << cell->j() << "\n";
//    }
    int imax = 0;
    int jmax = 0;
    for (auto cell : _cells) {
        if (cell->i() > imax){
            imax = cell->i();
        }
        if (cell->j() > jmax){
            jmax = cell->j();
        }
    }
    imax = imax-1; // because we have the ghost cells
    jmax = jmax-1;


    for (int j = 1; j < jmax; j++){
        field.v(0, j) = -field.v(1, j);
        field.u(0, j) = 0;
    }

    for (int j = 1; j < jmax; j++){
        field.v(imax, j) = -field.v(imax+1, j);
        field.u(imax, j) = 0;
    }

    for (int i = 1; i < imax; i++){
        field.v(i, 0) = 0;
        field.u(i, 0) = -field.u(i, 1);
    }

//    for (int i = 1; i < imax; i++){
//        field.v(i, jmax) = 0;
//        field.u(i, jmax) = -field.u(i, jmax+1);
//    }
}

void FixedWallBoundary::applyPressure(Fields &field) {}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : Boundary(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::applyVelocity(Fields &field) {
    int imax = 0;
    int jmax = 0;
    for (auto cell : _cells) {
        if (cell->i() > imax){
            imax = cell->i();
        }
        if (cell->j() > jmax){
            jmax = cell->j();
        }
    }
    imax = imax-1; // because we have the ghost cells
    jmax = jmax-1;


//    for (int j = 1; j < jmax; j++){
//        field.v(0, j) = -field.v(1, j);
//        field.u(0, j) = 0;
//    }
//
//    for (int j = 1; j < jmax; j++){
//        field.v(imax, j) = -field.v(imax+1, j);
//        field.u(imax, j) = 0;
//    }
//
//    for (int i = 1; i < imax; i++){
//        field.v(i, 0) = 0;
//        field.u(i, 0) = -field.u(i, 1);
//    }

    for (int i = 1; i < imax; i++){
        field.v(i, jmax) = 0;
        field.u(i, jmax) = -field.u(i, jmax+1);
    }
}

void MovingWallBoundary::applyPressure(Fields &field) {}
