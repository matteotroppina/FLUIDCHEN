clear # clear terminal
rm ../example_cases/ChannelWithObstacle/ChannelWithObstacle_Output/* #remove old output files
make && mpirun -np 5 ./fluidchen ../example_cases/ChannelWithObstacle/ChannelWithObstacle.dat 5 1