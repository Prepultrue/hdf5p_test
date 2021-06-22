# HDF5 MPI test

Quick test that ranks can see others writes

```shell
$ module load hdf5/1.10.7p
$ module load openmpi/4.0.3
$ h5pcc main.c 
$ mpiexec -n 2 a.out 
Writing
Rank 1: Selecting offset=[128, 0] count=[128, 5]
Rank 0: Selecting offset=[0, 0] count=[128, 5]

Reading
Rank 0: Selecting offset=[128, 0] count=[128, 5]
Rank 1: Selecting offset=[0, 0] count=[128, 5]
```
