# Numerical Methods and High Performance Computing Assignment
The WENO Schemes has not been successfully executed. There are errors generated in the F1 scheme, creating far different results to the analytical.

For efficiency, the codes for parallel code are simplified, thus generates less time overall compared to the serial codes and create more efficiency in transfering data from personal to the HPC facility. This is done by deleting the norms calculation, dividing into different files for every schemes, etc. 


F1.f90 and F0.f90 is for the nonparallel calculation. The other fortran files (which contains _parallel) is for parallel computing.

Input.dat is used for parallel input, and input_noparallel.dat is for no parallel input.
