/***************
 * Basil Lin
 * ECE 6730
 * README.txt
 ***************/
 
To compile serial programs:
gcc make-matrix.c -o make-matrix
gcc print-matrix.c -o print-matrix
gcc mv-serial.c -o mv-serial

To compile parallel program:
mpicc mv-parallel.c -o mv-parallel

Example run:
./make-matrix -n 16 -l 0 -u 100 -m -o matrix.bin
./make-matrix -n 16 -l 0 -u 100 -v -o vector.bin
./mv-serial matrix.bin vector.bin serial.bin
./print-matrix serial.bin
qsub mv-serial.pbs
./print-matrix parallel.bin

Notes:
- mv-serial.pbs is currently set to use 4 processes
- mv-serial.pbs requires matrix.bin, vector.bin, and
  outputs to parallel.bin
- serial.bin and parallel.bin will not exactly match
  due to small rounding imprecisions
- large files must be run using mpirun or mpiexec as
  Palmetto will not queue tasks too large
