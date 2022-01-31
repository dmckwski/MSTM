gfortran -O2 -fallow-argument-mismatch -c -o mstm-intrinsics.obj mstm-intrinsics.f90
gfortran -O2 -fallow-argument-mismatch -c -o mpidefs-serial.obj mpidefs-serial.f90
gfortran -O2 -fallow-argument-mismatch -c -o mstm.obj mstm.f90
gfortran -O2 -fallow-argument-mismatch -o mstm.exe mstm-intrinsics.obj mpidefs-serial.obj mstm.obj
