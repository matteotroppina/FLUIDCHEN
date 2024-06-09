#include <iostream>
#include <string>

#include "Case.hpp"
#include "Communication.hpp"

int main(int argn, char **args) {

    if (argn > 1) {
        std::string file_name{args[1]}; // input file name

        // save logs to file for each process
        Communication::init_parallel(argn, args);
        std::string outputlogs = "output_" + std::to_string(my_rank_global);
        std::freopen(outputlogs.c_str(), "w", stdout);


        Case problem(file_name, argn, args);
        problem.simulate();

    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
    Communication::finalize();
}
