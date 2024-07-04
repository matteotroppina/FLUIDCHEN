clear # clear terminal
rm ../example_cases/ChannelWithBFS/ChannelWithBFSUpscaled_Output/* #remove old output files
make && mpirun --oversubscribe -np 1 ./fluidchen ../example_cases/ChannelWithBFS/ChannelWithBFSUpscaled.dat 1 1