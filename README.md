# LaplacianSmoothing


I forked this repository from [Alec Jacobson](https://github.com/alecjacobson/geometry-processing-smoothing) and then I implemented the functions in order to perform the Laplacian Smoothing as learned in that repository notes. 

## Prerequisite installation

> #### ¹ Mac Users
>
> You will need to install Xcode if you haven't already. 
>
> #### ² Linux Users
>
> Many linux distributions do not include gcc and the basic development tools
> in their default installation. On Ubuntu, you need to install the following
> packages:
>
>     sudo apt-get install git
>     sudo apt-get install build-essential
>     sudo apt-get install cmake
>     sudo apt-get install libx11-dev
>     sudo apt-get install mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev
>     sudo apt-get install libxrandr-dev
>     sudo apt-get install libxi-dev
>     sudo apt-get install libxmu-dev
>     sudo apt-get install libblas-dev
>
>
> #### ³ Windows Users
>
> libigl only supports the Microsoft Visual Studio
> 2015 compiler in 64bit mode. It will not work with a 32bit build and it will
> not work with older versions of visual studio.


## Compilation and execution

>     mkdir build
>     cd build/
>     cmake ..
>     make -j 2 
>     ./smoothing [parameters]
