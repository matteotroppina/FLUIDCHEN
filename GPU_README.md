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


