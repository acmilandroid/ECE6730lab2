/*   Basil Lin
 *   ECE 6730
 *
 *	 MyMPI.h
 *
 *   Header file for a library of matrix/vector
 *   input/output/redistribution functions.
 *
 *   Programmed by Michael J. Quinn
 */

/************************* MACROS **************************/

#define DATA_MSG           0
#define PROMPT_MSG         1
#define RESPONSE_MSG       2

#define OPEN_FILE_ERROR    -1
#define MALLOC_ERROR       -2
#define TYPE_ERROR         -3

#define MIN(a,b)           ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) \
                     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j,p,n) (((p)*((j)+1)-1)/(n))
#define PTR_SIZE           (sizeof(void*))
#define CEILING(i,j)       (((i)+(j)-1)/(j))

/****************** INPUT FUNCTIONS ************************/

void read_checkerboard_matrix (char *, void ***, void **,
        MPI_Datatype, int *, int *, MPI_Comm);
void read_block_vector (char *, void **, MPI_Datatype,
        int *, MPI_Comm);

/****************** OUTPUT FUNCTIONS ***********************/
void print_block_vector (char *, void *, MPI_Datatype, int,
        MPI_Comm);
