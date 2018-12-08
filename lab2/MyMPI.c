/*	 Basil Lin
 *	 ECE 6730
 *
 *   MyMPI.c -- A library of matrix/vector
 *   input/output/redistribution functions
 *
 *   Programmed by Michael J. Quinn
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "MyMPI.h"


/***************** MISCELLANEOUS FUNCTIONS *****************/

/*
 *   Given MPI_Datatype 't', function 'get_size' returns the
 *   size of a single datum of that data type.
 */

int get_size (MPI_Datatype t) {
   if (t == MPI_BYTE) return sizeof(char);
   if (t == MPI_DOUBLE) return sizeof(double);
   if (t == MPI_FLOAT) return sizeof(float);
   if (t == MPI_INT) return sizeof(int);
   printf ("Error: Unrecognized argument to 'get_size'\n");
   fflush (stdout);
   MPI_Abort (MPI_COMM_WORLD, TYPE_ERROR);
}


/*
 *   Function 'my_malloc' is called when a process wants
 *   to allocate some space from the heap. If the memory
 *   allocation fails, the process prints an error message
 *   and then aborts execution of the program.
 */

void *my_malloc (
   int id,     /* IN - Process rank */
   int bytes)  /* IN - Bytes to allocate */
{
   void *buffer;
   if ((buffer = malloc ((size_t) bytes)) == NULL) {
      printf ("Error: Malloc failed for process %d\n", id);
      fflush (stdout);
      MPI_Abort (MPI_COMM_WORLD, MALLOC_ERROR);
   }
   return buffer;
}


/********************* INPUT FUNCTIONS *********************/

/*
 *   Function 'read_checkerboard_matrix' reads a matrix from
 *   a file. The first two elements of the file are integers
 *   whose values are the dimensions of the matrix ('m' rows
 *   and 'n' columns). What follows are 'm'*'n' values
 *   representing the matrix elements stored in row-major
 *   order.  This function allocates blocks of the matrix to
 *   the MPI processes.
 *
 *   The number of processes must be a square number.
 */
 
void read_checkerboard_matrix (
   char *s,              /* IN - File name */
   void ***subs,         /* OUT - 2D array */
   void **storage,       /* OUT - Array elements */
   MPI_Datatype dtype,   /* IN - Element type */
   int *m,               /* OUT - Array rows */
   int *n,               /* OUT - Array cols */
   MPI_Comm grid_comm)   /* IN - Communicator */
{
   void      *buffer;         /* File buffer */
   int        coords[2];      /* Coords of proc receiving
                                 next row of matrix */
   int        datum_size;     /* Bytes per elements */
   int        dest_id;        /* Rank of receiving proc */
   int        grid_coord[2];  /* Process coords */
   int        grid_id;        /* Process rank */
   int        grid_period[2]; /* Wraparound */
   int        grid_size[2];   /* Dimensions of grid */
   int        i, j, k;
   FILE      *infileptr;      /* Input file pointer */
   void      *laddr;          /* Used when proc 0 gets row */
   int        local_cols;     /* Matrix cols on this proc */
   int        local_rows;     /* Matrix rows on this proc */
   void     **lptr;           /* Pointer into 'subs' */
   int        p;              /* Number of processes */
   void      *raddr;          /* Address of first element
                                 to send */
   void      *rptr;           /* Pointer into 'storage' */
   MPI_Status status;         /* Results of read */

   MPI_Comm_rank (grid_comm, &grid_id);
   MPI_Comm_size (grid_comm, &p);
   datum_size = get_size (dtype);

   /* Process 0 opens file, gets number of rows and
      number of cols, and broadcasts this information
      to the other processes. */

   if (grid_id == 0) {
      infileptr = fopen (s, "r");
      if (infileptr == NULL) *m = 0;
      else {
         fread (m, sizeof(int), 1, infileptr);
         fread (n, sizeof(int), 1, infileptr);
      }
   }
   MPI_Bcast (m, 1, MPI_INT, 0, grid_comm);

   if (!(*m)) MPI_Abort (MPI_COMM_WORLD, OPEN_FILE_ERROR);

   MPI_Bcast (n, 1, MPI_INT, 0, grid_comm);

   /* Each process determines the size of the submatrix
      it is responsible for. */

   MPI_Cart_get (grid_comm, 2, grid_size, grid_period,
      grid_coord);
   local_rows = BLOCK_SIZE(grid_coord[0],grid_size[0],*m);
   local_cols = BLOCK_SIZE(grid_coord[1],grid_size[1],*n);

   /* Dynamically allocate two-dimensional matrix 'subs' */

   *storage = my_malloc (grid_id,
      local_rows * local_cols * datum_size);
   *subs = (void **) my_malloc (grid_id,local_rows*PTR_SIZE);
   lptr = (void *) *subs;
   rptr = (void *) *storage;
   for (i = 0; i < local_rows; i++) {
      *(lptr++) = (void *) rptr;
      rptr += local_cols * datum_size;
   }

   /* Grid process 0 reads in the matrix one row at a time
      and distributes each row among the MPI processes. */

   if (grid_id == 0)
      buffer = my_malloc (grid_id, *n * datum_size);

   /* For each row of processes in the process grid... */
   for (i = 0; i < grid_size[0]; i++) {
      coords[0] = i;

      /* For each matrix row controlled by this proc row...*/
      for (j = 0; j < BLOCK_SIZE(i,grid_size[0],*m); j++) {

         /* Read in a row of the matrix */

         if (grid_id == 0) {
            fread (buffer, datum_size, *n, infileptr);
         }

         /* Distribute it among process in the grid row */

         for (k = 0; k < grid_size[1]; k++) {
            coords[1] = k;

            /* Find address of first element to send */
            raddr = buffer +
               BLOCK_LOW(k,grid_size[1],*n) * datum_size;

            /* Determine the grid ID of the process getting
               the subrow */
            MPI_Cart_rank (grid_comm, coords, &dest_id);

            /* Process 0 is responsible for sending...*/
            if (grid_id == 0) {

               /* It is sending (copying) to itself */
               if (dest_id == 0) {
                  laddr = (*subs)[j];
                  memcpy (laddr, raddr,
                     local_cols * datum_size);

               /* It is sending to another process */
               } else {
                  MPI_Send (raddr,
                     BLOCK_SIZE(k,grid_size[1],*n), dtype,
                  dest_id, 0, grid_comm);
               }

            /* Process 'dest_id' is responsible for
               receiving... */
            } else if (grid_id == dest_id) {
               MPI_Recv ((*subs)[j], local_cols, dtype, 0,
                  0, grid_comm,&status);
            }
         }
      }
   }
   if (grid_id == 0) free (buffer);
}


/*
 *   Open a file containing a vector, read its contents,
 *   and distributed the elements by block among the
 *   processes in a communicator.
 */

void read_block_vector (
    char        *s,      /* IN - File name */
    void       **v,      /* OUT - Subvector */
    MPI_Datatype dtype,  /* IN - Element type */
    int         *n,      /* OUT - Vector length */
    MPI_Comm     comm)   /* IN - Communicator */
{
   int        datum_size;   /* Bytes per element */
   int        i;
   FILE      *infileptr;    /* Input file pointer */
   int        local_els;    /* Elements on this proc */
   MPI_Status status;       /* Result of receive */
   int        id;           /* Process rank */
   int        p;            /* Number of processes */
   int        x;            /* Result of read */
   int		  garbage;

   datum_size = get_size (dtype);
   MPI_Comm_size(comm, &p);
   MPI_Comm_rank(comm, &id);

   /* Process p-1 opens file, determines number of vector
      elements, and broadcasts this value to the other
      processes. */

   if (id == (p-1)) {
      infileptr = fopen (s, "r");
      if (infileptr == NULL) *n = 0;
      else {
      	fread (n, sizeof(int), 1, infileptr);
      	fread (&garbage, sizeof(int), 1, infileptr);
      }
   }
   MPI_Bcast (n, 1, MPI_INT, p-1, comm);
   if (! *n) {
      if (!id) {
         printf ("Input file '%s' cannot be opened\n", s);
         fflush (stdout);
      }
   }

   /* Block mapping of vector elements to processes */

   local_els = BLOCK_SIZE(id,p,*n);

   /* Dynamically allocate vector. */

   *v = my_malloc (id, local_els * datum_size);
   if (id == (p-1)) {
      for (i = 0; i < p-1; i++) {
         x = fread (*v, datum_size, BLOCK_SIZE(i,p,*n),
            infileptr);
         MPI_Send (*v, BLOCK_SIZE(i,p,*n), dtype, i, DATA_MSG,
            comm);
      }
      x = fread (*v, datum_size, BLOCK_SIZE(id,p,*n),
             infileptr);
      fclose (infileptr);
   } else {
      MPI_Recv (*v, BLOCK_SIZE(id,p,*n), dtype, p-1, DATA_MSG,
         comm, &status);
   }
}


/******************** OUTPUT FUNCTIONS ********************/

/*
 *   Print a vector that is block distributed among the
 *   processes in a communicator.
 */

void print_block_vector (
   char *arg,
   void        *v,       /* IN - Address of vector */
   MPI_Datatype dtype,   /* IN - Vector element type */
   int          n,       /* IN - Elements in vector */
   MPI_Comm     comm)    /* IN - Communicator */
{
   int        datum_size; /* Bytes per vector element */
   int        i;
   int        prompt;     /* Dummy variable */
   MPI_Status status;     /* Result of receive */
   void       *tmp;       /* Other process's subvector */
   int        id;         /* Process rank */
   int        p;          /* Number of processes */
   FILE *fpt;
   int one = 1;
   
   fpt = fopen(arg, "wb");
   fwrite(&n, sizeof(int), 1, fpt);
   fwrite(&one, sizeof(int), 1, fpt);

   MPI_Comm_size (comm, &p);
   MPI_Comm_rank (comm, &id);
   datum_size = get_size (dtype);

   if (!id) {
      fwrite((double *)v, sizeof(double), BLOCK_SIZE(id,p,n), fpt);
      if (p > 1) {
         tmp = my_malloc (id,BLOCK_SIZE(p-1,p,n)*datum_size);
         for (i = 1; i < p; i++) {
            MPI_Send (&prompt, 1, MPI_INT, i, PROMPT_MSG,
               comm);
            MPI_Recv (tmp, BLOCK_SIZE(i,p,n), dtype, i,
               RESPONSE_MSG, comm, &status);
            fwrite((double *)tmp, sizeof(double), BLOCK_SIZE(i,p,n), fpt);
         }
         free (tmp);
      }
   } else {
      MPI_Recv (&prompt, 1, MPI_INT, 0, PROMPT_MSG, comm,
         &status);
      MPI_Send (v, BLOCK_SIZE(id,p,n), dtype, 0,
         RESPONSE_MSG, comm);
   }
   fclose(fpt);
}
