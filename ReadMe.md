Project Overview
================================================================================

This project is a fork of the [NTPoly](https://github.com/william-dawson/NTPoly)
code aimed at providing drivers and data for the
[SPEC MPI Accelerator Benchmark Suite](https://www.spec.org/hpg/search/). Note
that this ReadMe only contains the information relevant for benchmarking.

Set Up Guide
--------------------------------------------------------------------------------
Installing NTPoly-SPEC suite requires the following software:

* A Fortran Compiler.
* An MPI Installation (MPI-3 Standard+).
* CMake (Version 3.2+).

The following optional software can greatly enhance the NTPoly experience:

* BLAS: for multiplying dense matrices, if they emerge in the calculation.

NTPoly-SPEC uses CMake as a build system. First, take a look in the Targets
directory. You'll find a list of `.cmake` files which have example configurations
on popular systems. You should copy one of these files, and create your own
mymachine.cmake file. Then, cd into the Build directory, and type:
> cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/mymachine.cmake -DFORTRAN_ONLY=YES ..

There are a few options you can pass to CMake to modify the build. You can set
`-DCMAKE_BUILD_TYPE=Debug` for debugging purposes. You can set the install
directory using the standard `-DCMAKE_INSTALL_PREFIX=/path/to/dir`. The flag
`-DFORTRAN_ONLY=YES` is set because we only need the Fortran bindings for
benchmarking. You can set the flag `-DNOIALLGATHER=Yes` if your MPI implementation
does not have nonblocking collectives.

After that you can build using:
> make

Driver Usage
--------------------------------------------------------------------------------
Drivers are to put into the `Driver` directory. Currently there is only one
driver, `InverseDriver`, which computes the inverse of a sparse matrix.
The following input parameters should be passed to this driver:

> --input followed by the name of the input matrix file.

> --process_rows the process cube dimension in the row direction.

> --process_columns the process cube dimension in the column direction.

> --process_slices the process cube dimension in the slice direction.

> --threshold for flushing small values to zero.

> --loop_times how many times to convert the inverse.

Note that the product of the process rows, columns, and slices must equal
the total number of processes. Slices should be kept low, and only increased
when strong scaling degrades. For optimal performance, keep the number of
rows and columns equal, or with one equal to two times the other.

For benchmarking purposes, the threshold should be set to `1e-5` or
`1e-6`. A smaller values will lead to denser matrices and longer compute
times. Input files are in the `Benchmarks` directory.

The final parameter `loop_times` controls how many times the inverse is
computed. In real world applications, we frequently need to compute the
same function on many different matrices with similar structure. By increasing
the loop_times, you can simulate this type of workload.

Example run:
> mpirun -np 4 ./bin/InverseDriver \
> --input ../Benchmarks/004.xyz-input631G.in.mtx  \
> --process_rows 2 --process_columns 2 --process_slices 1 \
> --threshold 1e-6 --loop_times 1

Citation
--------------------------------------------------------------------------------
A description of the techniques used in NTPoly can be found in the following
Computer Physics Communications paper:

> Dawson, William, and Takahito Nakajima. "Massively parallel sparse matrix
> function calculations with NTPoly." Computer Physics Communications (2017).

Please cite this paper in accordance to the practices in your field.
