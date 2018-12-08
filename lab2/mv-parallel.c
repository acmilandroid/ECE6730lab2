/******************************
 * Basil Lin
 * ECE 6730
 * mv-parallel.c
 ******************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <sys/types.h>
#include <time.h>
#include "MyMPI.c"
#include "MyMPI.h"


int main (int argc, char **argv) {
	int i, j, k; 
	char name[10]; 
	int offset; 
	int root; 
    int n;  /* total data items in each dim */
    int zero = 0;
	int *rtmp, *ctmp; 
	MPI_Status stat; 
	offset = 0; 
    int SIZE, RANK; /* original rank and size */
    int rank; /* local rank in the 2D topology */
    int coords[2]; /* local coordinates in the 2D topology */
    int row_start;
    int row_end;
    int row_cnt;
    int col_start;
    int col_end;
    int col_cnt;
    double **matrix;
    double *vector;
    double *local_vector;
    double **storage;
    double time;
    double max_time;
    double *out_vector_send;
    int nfinal;
    int garbage;

    /* size of topology grid in each direction */
	int gridsz[2] = {0, 0};  
    /* declare several communicators */
	MPI_Comm Comm_2D, Comm_row, Comm_col;
    /* arrays with info about the 2D topology */
    int period[2] = {0, 0}; /* do not wrap around cart dims */
/*————————————————————————————————————————————————————*/
    MPI_Init(&argc, &argv);

    /* must read n from file or user */
    FILE *fpt;
    if ((fpt = fopen(argv[1], "rb")) == NULL) {
    	printf("Could not open file %s. Exiting\n", argv[1]);
    	exit(0);
    }
    fread(&n, sizeof(int), 1, fpt);
    fclose(fpt);
 
    MPI_Comm_size(MPI_COMM_WORLD, &SIZE);
 	MPI_Comm_rank(MPI_COMM_WORLD, &RANK);

	/* get cart topology size */
    MPI_Dims_create(SIZE, 2, gridsz);
    /* create 2D cart topology */
	MPI_Cart_create(MPI_COMM_WORLD, 2, gridsz, period, 0, &Comm_2D);

	MPI_Comm_rank(Comm_2D, &rank);           /* get 2D rank … */
	MPI_Cart_coords(Comm_2D, rank, 2, coords);  /* then 2D coords */

	MPI_Comm_split(MPI_COMM_WORLD, coords[0], coords[1], &Comm_row); 
	MPI_Comm_split(MPI_COMM_WORLD, coords[1], coords[0], &Comm_col);
	

	/* find row start and end index, then same for column */
	row_start = BLOCK_LOW(coords[0],  gridsz[0], n);
	row_end   = BLOCK_HIGH(coords[0], gridsz[0], n);
	row_cnt   = BLOCK_SIZE(coords[0], gridsz[0], n); 
 
	col_start = BLOCK_LOW(coords[1],  gridsz[1], n);
	col_end   = BLOCK_HIGH(coords[1], gridsz[1], n); 
	col_cnt   = BLOCK_SIZE(coords[1], gridsz[1], n); 
	
	//read and print matrix
	read_checkerboard_matrix(argv[1], (void ***) &matrix, (void **) &storage, MPI_DOUBLE, &n, &n, Comm_2D);
	
	//read and print vector
	read_block_vector(argv[2], (void **) &vector, MPI_DOUBLE, &garbage, Comm_row);

	out_vector_send = (double *)my_malloc(rank, row_cnt*sizeof(double));
	MPI_Barrier(MPI_COMM_WORLD); //syncronizes processes
	
	time = -MPI_Wtime(); //start timing
	
	//multiply
	for (i = 0; i < row_cnt; i++) {
		out_vector_send[i] = 0.0;
		for (j = 0; j < col_cnt; j++) {
			out_vector_send[i] += matrix[i][j] * vector[j];
		}
	}
	
	double *out_vector_rcv = (double *)my_malloc(rank, row_cnt*sizeof(double));
	
	MPI_Reduce(out_vector_send, out_vector_rcv, col_cnt, MPI_DOUBLE, MPI_SUM, 0, Comm_row); //partial sum reduction	
	MPI_Barrier(MPI_COMM_WORLD); //syncronize processes
	
	time += MPI_Wtime(); //end timing
	MPI_Allreduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); //send time to first process

	//print results
	if (!rank) {
		printf("N = %d, Processes = %d, Time = %0.7f seconds\n", n, SIZE, max_time);
	}
	if(coords[1] == 0) {
		print_block_vector(argv[3], (void *) out_vector_rcv, MPI_DOUBLE, n, Comm_col);
	}
	
	MPI_Finalize();
	return 0;
}

