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
    double TI{};    /*initial temperature*/
    double alpha{}; /*thermal diffusivity*/
    double beta{};  /*coefficient of thermal expansion*/
    double wall_temp_3{};
    double wall_temp_4{};
    double wall_temp_5{};
    double KI{}; /*initial value for turbulent kinetic energy*/
    double EI{}; /*initial value for the dissipation rate*/
        
        // if(_turbulence && timestep % 10 < 1e-1){
        //         _field.printMatrix(_grid);
        //         cc++;
        // }
    int num_of_walls{};

    // initialized to sequential execution
    int iproc{1};
    int jproc{1};

    if (argn > 2){
        iproc = *args[2] - 48;  //convert to ASCII
        jproc = *args[3] - 48;
    }

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
                if (var == "TI") file >> TI;
                if (var == "alpha") file >> alpha;
                if (var == "beta") file >> beta;
                if (var == "turbulence") file >> _turbulence;
                if (var == "t_init") file >> _t_init; 
                if (var == "KI") file >> KI;
                if (var == "EI") file >> EI;
                // read geometry file name from .dat file and directly assign it to private member fo Case
                if (var == "geo_file") file >> _geom_name;
                if (var == "num_of_walls") file >> num_of_walls;
                if (var == "wall_temp_3") file >> wall_temp_3;
                if (var == "wall_temp_4") file >> wall_temp_4;
                if (var == "wall_temp_5") file >> wall_temp_5;
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

    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank_global == 0){
        std::cout << "\n(2/4) BUILDING DOMAINS...\n " << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    build_domain(domain, imax, jmax, iproc, jproc);

    MPI_Barrier(MPI_COMM_WORLD);

    // double l0 = 1;
    // double k0 = std::pow(nu / l0, 2);
    // double eps0 = _C0 * std::pow(k0, 1.5) / l0;

    // k0 = 0.1;
    // eps0 = 0.1;

    // std::cout << "k0 = " << k0 << std::endl;
    // std::cout << "eps0 = " << eps0 << std::endl;

    _grid = Grid(_geom_name, domain);
    _field = Fields(nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI, alpha, beta, GX, GY, TI, KI, EI);

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
    _viscosity_solver = std::make_unique<K_EPS_model>();
    _max_iter = itermax;
    _tolerance = eps;

    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank_global == 0){
        std::cout << "\n(3/4) READING GEOMETRY...\n" << std::endl;
        if(_geom_name.compare("NONE") != 0){
            std::cout << "Reading geometry file: " << _geom_name << std::endl;
        }else{
            std::cout << "No geometry file provided. Constructing boundaries for lid driven cavity" << std::endl;
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    if (_geom_name.compare("NONE") == 0) { // Construct boundaries for lid driven cavity

        if (not _grid.moving_wall_cells().empty()) {
            _boundaries.push_back(std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
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
        if (not _grid.fixed_velocity_cells().empty()) {
            _boundaries.push_back(std::make_unique<FixedVelocityBoundary>(_grid.fixed_velocity_cells(), UIN, VIN));
        }
        if (not _grid.zero_gradient_cells().empty()) {
            _boundaries.push_back(std::make_unique<ZeroGradientBoundary>(_grid.zero_gradient_cells()));
        }
        if (not _grid.fixed_wall_cells().empty()) {
            _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells(), wall_temp_3));
        }
        if (not _grid.hot_wall_cells().empty()) {
            _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.hot_wall_cells(), wall_temp_4));
        }
        if (not _grid.cold_wall_cells().empty()) {
            _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.cold_wall_cells(), wall_temp_5));
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

void Case::simulate() {

    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank_global == 0){
        std::cout << "\n(4/4) FLUIDCHEN SIMULATION...\n" << std::endl;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    double t = 0.0;
    double dt = _field.dt();
    int timestep = 0;
    double output_counter = 0.0;

    double residual = 1;
    int iter = 0;
    std::vector<int> iter_vec;
    bool turbulence_started = false;

    while (t < _t_end) {

        _field.calculate_dt(_grid, turbulence_started);
        dt = _field.dt();
        
        for (auto &b : _boundaries) {
            b->applyVelocity(_field);
            b->applyTemperature(_field);
            
        }


        _field.calculate_temperature(_grid);
        Communication::communicate(_field.t_matrix());


        _field.calculate_fluxes(_grid, turbulence_started);
        Communication::communicate(_field.f_matrix());
        Communication::communicate(_field.g_matrix());


        for (auto &b : _boundaries) {
            b->applyFlux(_field);
        }

        _field.calculate_rs(_grid);

        residual = 1;
        iter = 0;
        while (iter < _max_iter and residual > _tolerance) {
            residual = _pressure_solver->solve(_field, _grid);
            Communication::communicate(_field.p_matrix());

            for (auto &b : _boundaries) {
                b->applyPressure(_field);
            }
            iter += 1;

            residual = Communication::reduce_sum(residual);
        }

        iter_vec.push_back(iter);

        _field.calculate_velocities(_grid);
        Communication::communicate(_field.u_matrix());
        Communication::communicate(_field.v_matrix());

        if(_turbulence && t > _t_init){

            turbulence_started = true;

            _viscosity_solver->solve(_field, _grid);
            _field.calculate_nuT(_grid, _C0);


            for(auto &b : _boundaries){
                b->applyTurbulence(_field);
            }

            Communication::communicate(_field.k_matrix());
            Communication::communicate(_field.e_matrix());
            Communication::communicate(_field.nuT_matrix());
            Communication::communicate(_field.nuT_i_matrix());
            Communication::communicate(_field.nuT_j_matrix());
            
        }



        // TO DO: here turbulence loop, only enter if a certain t value is reached? at the end: replace nu with nu+nuT from viscosity solver
        
        // if(_turbulence && timestep % 10 < 1e-1){
        //         _field.printMatrix(_grid);
        // }
        

        timestep += 1;
        output_counter += dt;
        t += dt;

        if (output_counter >= _output_freq or timestep == 1) {

            output_counter = 0;

            // double max_p = _field.p_matrix().max_abs_value();
            // if (max_p > 1e6 or max_p != max_p or residual != residual) { // check larger than or nan
            //     std::cerr << "Divergence detected" << std::endl;
            //     break;
            // }

            output_vtk(timestep, my_rank_global);
            if (my_rank_global == 0) {
                std::cout << "\n[" << static_cast<int>((t / _t_end) * 100) << "%"
                          << " completed] " << "Writing Output at t = " << t << "s" << std::endl;
                std::cout << std::left << "[ " << "Timestep: " << timestep << "\t\tSOR Iterations: " << iter << "\tSOR Residual: " << residual << " ]"<< std::flush;
                if (iter == _max_iter) {
                    std::cout << "\t\t ---> Exceeded max iterations";
                }
                std::cout << "\n------------------------------------------------------------------------------------" << std::flush;
                // std::cout << "min/max p: " << _field.p_matrix().min_value() << " / " << _field.p_matrix().max_value()
                //           << std::endl;
            }

        }

        // output for performance analysis - comment the output above
        
        // if (output_counter >= _output_freq && my_rank_global == 0) {
        //     std::cout << "\n[" << static_cast<int>((t / _t_end) * 100) << "%" << " completed] " << std::endl; 
        //     output_counter = 0;
        // }

    }
    // output_csv(iter_vec);

    if (my_rank_global == 0) {
        std::cout << "\n\n[100% completed] Simulation completed successfully!\n" << std::endl;
    }
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
            if (_grid.cell(row + 1, col + 1).type() == cell_type::FLUID) {
                continue;
            }
            structuredGrid->BlankCell(row + col * (_grid.domain().size_x));
        }
    }

    // Pressure Array
    vtkSmartPointer<vtkDoubleArray> Pressure = vtkSmartPointer<vtkDoubleArray>::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Temperature Array
    vtkSmartPointer<vtkDoubleArray> Temperature = vtkSmartPointer<vtkDoubleArray>::New();
    Temperature->SetName("Temperature");
    Temperature->SetNumberOfComponents(1);

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
            double temperature = _field.T(i, j);
            Pressure->InsertNextTuple(&pressure);
            Temperature->InsertNextTuple(&temperature);
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
    structuredGrid->GetCellData()->AddArray(Temperature);

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

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain, int iproc, int jproc) {

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Building domain for process: " << my_rank_global << std::endl;

    int i = my_coords_global[0];
    int j = my_coords_global[1];

    int size_x = imax_domain / iproc;
    int size_y = jmax_domain / jproc;

    domain.size_x = size_x;
    domain.size_y = size_y;

    domain.itermax_x = size_x;
    domain.itermax_y = size_y;

    domain.iminb = i * size_x;
    domain.jminb = j * size_y;
    domain.imaxb = (i+1) * size_x + 2;
    domain.jmaxb = (j+1) * size_y + 2;

    std::array<int, 4> neighbours = Communication::get_neighbours();

    if (neighbours[RIGHT] != MPI_PROC_NULL){ // if there is a right neighbour
        domain.itermax_x = size_x + 1;
        // std::cout << "rank: " << my_rank_global << " has right neighbour" << std::endl;
    }

    if (neighbours[UP] != MPI_PROC_NULL){ // if there is an upper neighbour
        domain.itermax_y = size_y + 1;
        // std::cout << "rank: " << my_rank_global << " has upper neighbour" << std::endl;
    }

}
