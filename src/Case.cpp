#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace filesystem = std::filesystem;

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

#include "Case.hpp"
#include "Enums.hpp"

Case::Case(std::string file_name, int argn, char **args) {
    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    std::ifstream file(file_name);
    double nu{};      /* viscosity   */
    double UI{};      /* velocity x-direction */
    double VI{};      /* velocity y-direction */
    double PI{};      /* pressure */
    double GX{};      /* gravitation x-direction */
    double GY{};      /* gravitation y-direction */
    double xlength{}; /* length of the domain x-dir.*/
    double ylength{}; /* length of the domain y-dir.*/
    double dt{};      /* time step */
    int imax{};       /* number of cells x-direction*/
    int jmax{};       /* number of cells y-direction*/
    double gamma{};   /* uppwind differencing factor*/
    double omg{};     /* relaxation factor */
    double tau{};     /* safety factor for time step*/
    int itermax{};    /* max. number of iterations for pressure per time step */
    double eps{};     /* accuracy bound for pressure*/
    double UIN{};     /*inlet velocity*/
    double VIN{};
    int num_of_walls{};

    if (file.is_open()) {

        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
            if (var[0] == '#') { /* ignore comment line*/
                file.ignore(MAX_LINE_LENGTH, '\n');
            } else {
                if (var == "xlength") file >> xlength;
                if (var == "ylength") file >> ylength;
                if (var == "nu") file >> nu;
                if (var == "t_end") file >> _t_end;
                if (var == "dt") file >> dt;
                if (var == "omg") file >> omg;
                if (var == "eps") file >> eps;
                if (var == "tau") file >> tau;
                if (var == "gamma") file >> gamma;
                if (var == "dt_value") file >> _output_freq;
                if (var == "UI") file >> UI;
                if (var == "VI") file >> VI;
                if (var == "GX") file >> GX;
                if (var == "GY") file >> GY;
                if (var == "PI") file >> PI;
                if (var == "itermax") file >> itermax;
                if (var == "imax") file >> imax;
                if (var == "jmax") file >> jmax;
                if (var == "UIN") file >> UIN;
                if (var == "VIN") file >> VIN;
                // read geometry file name from .dat file and directly assign it to private member fo Case
                if (var == "geo_file") file >> _geom_name;
                if (var == "num_of_walls") file >> num_of_walls;
            }
        }
    }
    file.close();

    if (_geom_name.compare("NONE") == 0) {
        std::map<int, double> wall_vel;
        wall_vel.insert(std::pair<int, double>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    }

    // Set file names for geometry file and output directory
    set_file_names(file_name);

    // Build up the domain
    Domain domain;
    domain.dx = xlength / static_cast<double>(imax);
    domain.dy = ylength / static_cast<double>(jmax);
    domain.domain_imax = imax;
    domain.domain_jmax = jmax;

    build_domain(domain, imax, jmax);

    _grid = Grid(_geom_name, domain);
    _field = Fields(nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI);

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
    _max_iter = itermax;
    _tolerance = eps;

    if (_geom_name.compare("NONE") == 0) { // Construct boundaries for lid driven cavity
        if (not _grid.moving_wall_cells().empty()) {
            _boundaries.push_back(
                std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
        }
        if (not _grid.fixed_wall_cells().empty()) {
            _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells()));
        }
    } else { // general case
        if (not _grid.inner_obstacle_cells().empty()) {
            _boundaries.push_back(std::make_unique<InnerObstacle>(_grid.inner_obstacle_cells()));
        }
        if (not _grid.moving_wall_cells().empty()) {
            //            _boundaries.push_back(std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), ??)); //
            //            TODO: set wall velocity according to input file
        }
        if (not _grid.fixed_wall_cells().empty()) {
            _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells()));
        }
        if (not _grid.fixed_velocity_cells().empty()) {
            _boundaries.push_back(std::make_unique<FixedVelocityBoundary>(_grid.fixed_velocity_cells(), UIN, VIN));
        }
        if (not _grid.zero_gradient_cells().empty()) {
            _boundaries.push_back(std::make_unique<ZeroGradientBoundary>(_grid.zero_gradient_cells()));
        }
    }
}

void Case::set_file_names(std::string file_name) {
    std::string temp_dir;
    bool case_name_flag = true;
    bool prefix_flag = false;

    for (int i = file_name.size() - 1; i > -1; --i) {
        if (file_name[i] == '/') {
            case_name_flag = false;
            prefix_flag = true;
        }
        if (case_name_flag) {
            _case_name.push_back(file_name[i]);
        }
        if (prefix_flag) {
            _prefix.push_back(file_name[i]);
        }
    }

    for (int i = file_name.size() - _case_name.size() - 1; i > -1; --i) {
        temp_dir.push_back(file_name[i]);
    }

    std::reverse(_case_name.begin(), _case_name.end());
    std::reverse(_prefix.begin(), _prefix.end());
    std::reverse(temp_dir.begin(), temp_dir.end());

    _case_name.erase(_case_name.size() - 4);
    _dict_name = temp_dir;
    _dict_name.append(_case_name);
    _dict_name.append("_Output");

    if (_geom_name.compare("NONE") != 0) {
        _geom_name = _prefix + _geom_name;
    }

    // Create output directory
    filesystem::path folder(_dict_name);
    try {
        filesystem::create_directory(folder);
    } catch (const std::exception &e) {
        std::cerr << "Output directory could not be created." << std::endl;
        std::cerr << "Make sure that you have write permissions to the "
                     "corresponding location"
                  << std::endl;
    }
}

/**
 * This function is the main simulation loop. In the simulation loop, following steps are required
 * Calculate and apply velocity boundary conditions for all the boundaries in _boundaries container
 * using applyVelocity() member function of Boundary class
 * Calculate fluxes (F and G) using calculate_fluxes() member function of Fields class.
 * Flux consists of diffusion and convection part, which are located in Discretization class
 * Apply Flux boundary conditions using applyFlux()
 * Calculate right-hand-side of PPE using calculate_rs() member function of Fields class
 * - Iterate the pressure poisson equation until the residual becomes smaller than the desired tolerance
 *   or the maximum number of the iterations are performed using solve() member function of PressureSolver
 * - Update pressure boundary conditions after each iteration of the SOR solver
 * - Calculate the velocities u and v using calculate_velocities() member function of Fields class
 * - calculate the maximal timestep size for the next iteration using calculate_dt() member function of Fields class
 * - Write vtk files using output_vtk() function
 *
 * Please note that some classes such as PressureSolver, Boundary are abstract classes which means they only provide the
 * interface and/or common functions. You need to define functions with individual functionalities in inherited
 * classes such as MovingWallBoundary class.
 * For information about the classes and functions, you can check the header files.
 */
void Case::simulate() {

    double t = 0.0;
    double dt = _field.dt();
    int timestep = 0;
    double output_counter = 0.0;

    double residual = 1;
    int iter = 0;
    std::vector<int> iter_vec;

    _field.printCellTypes(_grid);
    _field.printBorders(_grid);

    while (t < _t_end) {

        _field.calculate_dt(_grid);
        dt = _field.dt();

        for (auto &b : _boundaries) {
            b->applyVelocity(_field);
        }

        _field.calculate_fluxes(_grid);

        for (auto &b : _boundaries) {
            b->applyFlux(_field);
        }

        _field.calculate_rs(_grid);

        residual = 1;
        iter = 0;
        while (iter < _max_iter and residual > _tolerance) {
            for (auto &b : _boundaries) {
                b->applyPressure(_field);
            }
            residual = _pressure_solver->solve(_field, _grid, _boundaries);
            iter += 1;
        }

        iter_vec.push_back(iter);

        if (iter == _max_iter) {
            std::cerr << "Exceeded max iterations" << std::endl;
        }

        _field.calculate_velocities(_grid);

        timestep += 1;
        output_counter += dt;
        t += dt;

        if (output_counter + dt/2 >= _output_freq or timestep == 1) {
            for (auto &b : _boundaries) {
                b->applyVelocity(_field); // for debugging
            }

            std::cout << "time: " << t << " timestep: " << timestep << " residual: " << residual << std::endl;
            std::cout << "min/max p: " << _field.p_matrix().min_value() << " / " << _field.p_matrix().max_value()
                      << std::endl;

            double max_p = _field.p_matrix().max_abs_value();
            if (max_p > 1e6 or max_p != max_p) { // check larger than or nan
                std::cerr << "Divergence detected" << std::endl;
                break;
            }
            output_vtk(timestep, 0);
            output_counter = output_counter - _output_freq;

            _field.printMatrix(_grid); // remove this line to run normally
        }
    }
     output_csv(iter_vec);
}

void Case::output_csv(const std::vector<int> &vec) {
    std::string filename = _dict_name + "/iterations.csv";

    std::ofstream file(filename);
    if (file.is_open()) {
        for (size_t i = 0; i < vec.size(); ++i) {
            file << vec[i];
            if (i != vec.size() - 1) {
                file << ",";
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void Case::output_vtk(int timestep, int my_rank) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    double dx = _grid.dx();
    double dy = _grid.dy();

    double x = _grid.domain().iminb * dx;
    double y = _grid.domain().jminb * dy;

    { y += dy; }
    { x += dx; }

    double z = 0;
    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        x = _grid.domain().iminb * dx;
        { x += dx; }
        for (int row = 0; row < _grid.domain().size_x + 1; row++) {
            points->InsertNextPoint(x, y, z);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(_grid.domain().size_x + 1, _grid.domain().size_y + 1, 1);
    structuredGrid->SetPoints(points);

    // Set blank cells
    for (int col = 0; col < _grid.domain().size_y; col++) {
        for (int row = 0; row < _grid.domain().size_x; row++) {
            if (_grid.cell(row+1, col+1).type() == cell_type::FLUID) {
                continue;
            }
            structuredGrid->BlankCell(row + col * (_grid.domain().size_x));
        }
    }

    // Pressure Array
    vtkSmartPointer<vtkDoubleArray> Pressure = vtkSmartPointer<vtkDoubleArray>::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Velocity Array for cell data
    vtkSmartPointer<vtkDoubleArray> Velocity = vtkSmartPointer<vtkDoubleArray>::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print pressure, velocity and temperature from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double pressure = _field.p(i, j);
            Pressure->InsertNextTuple(&pressure);
            vel[0] = (_field.u(i - 1, j) + _field.u(i, j)) * 0.5;
            vel[1] = (_field.v(i, j - 1) + _field.v(i, j)) * 0.5;
            Velocity->InsertNextTuple(vel);
        }
    }

    // Velocity Array for point data
    vtkSmartPointer<vtkDoubleArray> VelocityPoints = vtkSmartPointer<vtkDoubleArray>::New();
    VelocityPoints->SetName("velocity");
    VelocityPoints->SetNumberOfComponents(3);



    // Print Velocity from bottom to top
    for (int j = 0; j < _grid.domain().size_y + 1; j++) {
        for (int i = 0; i < _grid.domain().size_x + 1; i++) {
            vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
            vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
            VelocityPoints->InsertNextTuple(vel);
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured Grid
    structuredGrid->GetCellData()->AddArray(Velocity);
    structuredGrid->GetPointData()->AddArray(VelocityPoints);

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname =
        _dict_name + '/' + _case_name + "_" + std::to_string(my_rank) + "." + std::to_string(timestep) + ".vtk";

    //    std::cout << "saving output to : " << outputname << std::endl;

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain) {
    domain.iminb = 0;
    domain.jminb = 0;
    domain.imaxb = imax_domain + 2;
    domain.jmaxb = jmax_domain + 2;
    domain.size_x = imax_domain;
    domain.size_y = jmax_domain;
}
