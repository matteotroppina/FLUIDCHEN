clear  # clear terminal
rm ../example_cases/ShearFlow/ShearFlowUpscaled2_Output/* #remove old output files
#make && mpirun -np 1  ./fluidchen ../example_cases/ShearFlow/ShearFlowUpscaled2.dat 1 1
make && mpirun -np 10  ./fluidchen ../example_cases/ShearFlow/ShearFlowUpscaled2.dat 5 2