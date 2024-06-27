clear # clear terminal
rm -rf ChannelWithBFS_Output/* #remove old output files
make && mpirun --oversubscribe -np 6 ./fluidchen ChannelWithBFS.dat 3 2