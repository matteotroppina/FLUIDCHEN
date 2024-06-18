# Installation instructions for GPU

We are using OpenACC for GPU acceleration. OpenACC is a directive-based programming model for programming accelerators. 
OpenACC requires a compiler that supports OpenACC. 
Since we are using NVIDIA GPUs, we will use the NVIDIA HPC SDK compiler: [Official NVIDIA HPC SDK installation guide](https://docs.nvidia.com/hpc-sdk//hpc-sdk-install-guide/index.html).

Make sure to check your Nvidia driver version and CUDA version before installing the HPC SDK. 
Since I have an RTX 2080Ti and an older driver with support for CUDA 12.0, I will be using the SDK version 23.3. All releases can be found on
[NVIDIA HPC SDK Releases](https://developer.nvidia.com/nvidia-hpc-sdk-releases).

You will need to unpack the tarball and run the installation script (with sudo). The default installation directory will be `/opt/nvidia/hpc_sdk`.

After a successful install, you will need to add the compilers and OpenMPI to your PATH. If you're running on a cluster, follow the instructions in the documentation.

In bash, sh, or ksh, use these commands:
You can add them to your `.bashrc` or `.bash_profile` file to make them available every time you open a terminal.

`MANPATH=$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/compilers/man; export MANPATH`
`PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/compilers/bin:$PATH; export PATH`

Once the 64-bit compilers are available, you can make the OpenMPI
commands and man pages accessible using these commands

`export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/comm_libs/mpi/bin:$PATH`
`export MANPATH=$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/comm_libs/mpi/man`

Finally, invoke the `nvaccelinfo` command to see that your GPU and drivers are properly installed and available. 
At the end of the output, you will see something like: `Default Target:                cc75`. 
Use the compiler flag `-gpu=cc75` to specifically target a certain class of GPU. If you do not specify the GPU architecture, 
the compiler will target the default GPU that is installed on your system.

To compile the code with NVC++ compiler and OpenACC support, use the following cmake command:
```Shell
cmake -DCMAKE_BUILD_TYPE=GPU ..
```

# Debugging
For debugging, read the following reply from an Nvidia forum ([@MatColgrove](https://forums.developer.nvidia.com/t/how-to-debug-illegal-address-during-kernel-execution-q/135831)):
> You can use cuda-gdb to debug the device code. Just compile with “-g” and you’ll be able to step into the generated CUDA kernels. [...]
> 
> 1. Compile the OpenACC regions to target the CPUs (-ta=multicore). You can then run this code in pgdbg to track down any issues with the parallel code.
> 2. Compile with CUDA Unified Memory (-ta=tesla:managed). This will put all dynamic memory allocations in an address space accessible from both the host and device. If the program works, then it’s an indication that the code is accessing a host address.
> 3. Set the environment variable PGI_ACC_NOTIFY=1 (or the more verbose PGI_ACC_DEBUG=1). This will print out all the kernel calls made to the device and help track down which one is failing.

If the cuda-gdb is not working, try exporting the following environment variable:
```Shell
export CUDBG_USE_LEGACY_DEBUGGER=1
```
for timing information, you can use the following environment variable:
```Shell
export NV_ACC_TIME=1
```

In my case I couldn't run without mpirun, so even the cuda debugger needs to be run with mpirun:
```Shell
mpirun -np 1 cuda-gdb fluidchen
````
Then, inside the cuda-gdb, you can run the program with the input file:
```Shell
(cuda-gdb) r ../example_cases/LidDrivenCavity/LidDrivenCavity.dat
```


