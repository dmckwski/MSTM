gfortran -O2 -fallow-argument-mismatch -c -o mstm-intrinsics.obj mstm-intrinsics.f90 -IC:\lib\mpi\include -LC:\lib\mpi\lib -lmsmpi
gfortran -O2 -fallow-argument-mismatch -c -o mpidefs-parallel.obj mpidefs-parallel.f90 -IC:\lib\mpi\include -LC:\lib\mpi\lib -lmsmpi
gfortran -O2 -fallow-argument-mismatch -c -o mstm.obj mstm.f90 -IC:\lib\mpi\include -LC:\lib\mpi\lib -lmsmpi
gfortran -O2 -fallow-argument-mismatch -o mstm.exe mstm-intrinsics.obj mpidefs-parallel.obj mstm.obj -IC:\lib\mpi\include -LC:\lib\mpi\lib -lmsmpi
