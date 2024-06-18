#include <iostream>
#include <string>

#include "Case.hpp"
#include "Communication.hpp"
#include <chrono>

int main(int argn, char **args) {

    auto start = std::chrono::steady_clock::now();

    if (argn > 1) {
        std::string file_name{args[1]}; // input file name

        Communication::init_parallel(argn, args);
        Case problem(file_name, argn, args);
        problem.simulate();

    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
    Communication::finalize();

    if(my_rank_global == 0){
        auto end = std::chrono::steady_clock::now();
        std::cout << "Simulation Runtime:" << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s\n\n";
    }
}
