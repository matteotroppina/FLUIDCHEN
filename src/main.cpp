#include <iostream>
#include <string>

#include "Case.hpp"
#include "Communication.hpp"

int main(int argn, char **args) {

    if (argn > 1) {
        std::string file_name{args[1]};
        int iproc{*args[2]};
        int jproc{*args[3]};
        iproc = iproc - 48; // convert ascii to integer
        jproc = jproc - 48;

        std::cout << "iproc: " << iproc << std::endl;
        std::cout << "jproc: " << jproc << std::endl;

        init_parallel(iproc, jproc);

        Case problem(file_name, argn, args);
        problem.simulate();
    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
}
