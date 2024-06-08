clear # clear terminal
cd build
rm ../example_cases/ChannelWithObstacle/ChannelWithObstacle_Output/*
make && mpirun -np 5 ./fluidchen ../example_cases/ChannelWithObstacle/ChannelWithObstacle.dat 5 1
# mpirun -np 2 ./fluidchen ../example_cases/ShearFlow/ShearFlow.dat 2 1