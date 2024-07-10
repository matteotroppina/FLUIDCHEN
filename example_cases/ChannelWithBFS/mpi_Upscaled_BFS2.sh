clear # clear terminal
rm ../example_cases/ChannelWithBFS/ChannelWithBFSUpscaled2_Output/* #remove old output files
make && mpirun --oversubscribe -np 1 ./fluidchen ../example_cases/ChannelWithBFS/ChannelWithBFSUpscaled2.dat 1 1
