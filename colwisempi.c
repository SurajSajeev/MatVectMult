#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef int dtype;
#define BLOCK_OWNER(j,p,n) (((p)*((j)+1)-1)/(n))
#define MIN(a,b)           ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))			// low_value= (in/p)
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)		// high value (i+1)n/p-1
#define BLOCK_SIZE(id,p,n) \
                     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)	// size of the interval
#define PTR_SIZE           (sizeof(void*))			// size of a pointer. remeber sizeof(int *) = sizeof(float *) = sizeof(void *) = .. (all pointer types have same size)

void writevector(char *filename,size_t * vector,int n){
    FILE *fptr=fopen(filename,"w");
    for(int i=0;i<n;i++)
    fprintf(fptr,"%lu\n",vector[i]);
    fclose(fptr);
}

void writevector2(char *filename,int * vector,int n){
    FILE *fptr=fopen(filename,"w");
    for(int i=0;i<n;i++)
    fprintf(fptr,"%d\n",vector[i]);
    fclose(fptr);
}
size_t **matrix_alloc_for_p2(size_t r,size_t c){
int len=0; 
    size_t *ptr, **arr; 
    int count = 0,i,j; 
    FILE *fptr;
    fptr = fopen("matrix.txt","w");
    len =  c * r; 
    arr = (size_t **)malloc(r*sizeof(size_t*)); 
    ptr=(size_t *)malloc(len*(sizeof(size_t)));
    // ptr is now pointing to the first element in of 2D array 
    // for loop to point rows pointer to appropriate location in 2D array 
    for(i = 0; i < r; i++) 
        arr[i] = (ptr + c * i); 
  
    for (i = 0; i < r; i++){ 
        for (j = 0; j < c; j++){ 
            arr[i][j] = count++; // OR *(*(arr+i)+j) = ++count 
        fprintf(fptr,"%lu ",arr[i][j]);
        }
    fprintf(fptr,"\n");
    }
    fclose(fptr);
    return arr; 
}


void *my_malloc ( int id, int bytes)
{	void *buffer;
	if ((buffer = malloc ((size_t) bytes)) == NULL)
	{	
		MPI_Abort (MPI_COMM_WORLD, -2);
	}
	return buffer;
}

void createmixedxferarrays ( int id, int p, int n, int *count, int *disp)
{	int i;
	count[0] = BLOCK_SIZE(0,p,n);
	disp[0] = 0;
	for (i = 1; i < p; i++)
	{	disp[i] = disp[i-1] + count[i-1];
		count[i] = BLOCK_SIZE(i,p,n);
	}
}

void createuniformparts ( int id, int p, int n, int *count, int *disp)
{	int i;
	count[0] = BLOCK_SIZE(id,p,n);
	disp[0] = 0;
	for (i = 1; i < p; i++)
	{	disp[i] = disp[i-1] + count[i-1];
		count[i] = BLOCK_SIZE(id,p,n);
	}
}




void read_block_vector ( dtype *v, int n, MPI_Comm comm)
{	int i, j;
	MPI_Status status;       /* Result of receive */
	int my_rank, p;
	MPI_Comm_size(comm, &p);
	MPI_Comm_rank(comm, &my_rank);
	MPI_Barrier(comm);
	/* Process p-1 opens file, determines number of vector elements, and broadcasts this value to the other processes. */
	if (my_rank == (p-1))
	{	for (i = 0; i < p-1; i++)
		{	for ( j=0; j < BLOCK_SIZE(i,p,n); j++ )
				v[j]= 1;
			MPI_Send (v, BLOCK_SIZE(i,p,n), MPI_INT, i, i, comm);
		}
		for ( j=0; j < BLOCK_SIZE(p-1,p,n); j++ )
			v[j]= 1;
	}
	else
		MPI_Recv (v, BLOCK_SIZE(my_rank,p,n), MPI_INT, p-1, my_rank, comm, &status);
	MPI_Barrier(comm);
}

void print_vector(size_t *vector,int n){
	for(int i=0;i<n;i++)
	printf("%lu ",vector[i]);
	printf("\n");
}



void read_col_striped_matrix ( dtype **subs, MPI_Datatype dtype, int m, int n, MPI_Comm comm, int local_cols)
{	int i, my_rank, p,j;
	MPI_Status status;
	MPI_Comm_size (comm, &p);
	MPI_Comm_rank (comm, &my_rank);
	int *send_count = my_malloc (my_rank, p * sizeof(int));
	int *send_disp = my_malloc (my_rank, p * sizeof(int));
	int storage[n];
	createmixedxferarrays (my_rank,p,n,send_count,send_disp);
	MPI_Barrier(comm);
	{	int i1,j;
        int temp=0;
		for ( j=0; j< m; j++ )
		{	if ( my_rank == p-1 )
			{	for ( i1=0; i1 < n; i1++)
					storage[i1]=temp++;
			}
			MPI_Scatterv (storage, send_count, send_disp, MPI_INT, subs[j], local_cols, MPI_INT, p-1, comm);
		}
		MPI_Barrier(comm);
	}
}
void replicate_block_vector_for_c ( size_t *ablock, int n, size_t *arep, MPI_Comm comm)
{	int my_rank, p;
	MPI_Comm_size (comm, &p);
	MPI_Comm_rank (comm, &my_rank);
	int *cnt = (int*) malloc (p*sizeof(int));
	int *disp = (int*) malloc (p*sizeof(int));
	createmixedxferarrays(my_rank, p, n, cnt, disp);				// deciding the order and initializing the 2 arrays
	MPI_Allgatherv (ablock, cnt[my_rank], MPI_LONG, arep, cnt, disp, MPI_LONG, comm);
	free (cnt);
	free (disp);
}

int main (int argc, char *argv[])
{	dtype **a;				/* matrix */
	dtype *b;				/* vector */
	size_t *c;				/* The product, a vector */
	size_t *c_part_out;			/* Partial sums, sent */
	size_t  *c_part_in;			/* Partial sums, received */
	int *cnt_out;				/* Elements sent to each proc */
	int *cnt_in;				/* Elements received per proc */
	int *disp_out;				/* Indices of sent elements */
	int *disp_in;				/* Indices of received elements */
	int i, j, my_rank, m, n, p, nb, local_cols, size_of_type;
	double max_seconds, seconds;
	dtype *rptr, *local_a;			/* This process's portion of 'a' */
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size (MPI_COMM_WORLD, &p);
	
	if ( my_rank == (p-1))
	{	m=atoi(argv[1]);
		n=atoi(argv[2]);
		nb=1000;
	}
	MPI_Bcast (&m, 1, MPI_INT, p-1, MPI_COMM_WORLD);
	MPI_Bcast (&n, 1, MPI_INT, p-1, MPI_COMM_WORLD);
	//MPI_Bcast (&nb, 1, MPI_INT, p-1, MPI_COMM_WORLD);
	matrix_alloc_for_p2(m,n);
	size_of_type = sizeof(int);
	local_cols = BLOCK_SIZE(my_rank,p,n);
	local_a = my_malloc (my_rank, local_cols * m * size_of_type);
	b = my_malloc (my_rank, local_cols * size_of_type);
	a = my_malloc (my_rank, m * PTR_SIZE);				// a will store pointer to a pointer type!!
	rptr = (dtype *)local_a;
	a[0]= rptr;
	for (i = 1; i < m; i++)
	{	rptr += (local_cols);					//size_of_type;		// go to next row
		a[i]= (dtype *) (rptr);
	}
	read_col_striped_matrix ( a, MPI_INT, m, n, MPI_COMM_WORLD, local_cols);
	read_block_vector (b, nb, MPI_COMM_WORLD);
	writevector2("vector.txt",b,n);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	c_part_out = (size_t *) my_malloc (my_rank, n * sizeof(size_t));
	MPI_Barrier (MPI_COMM_WORLD);
	seconds = -MPI_Wtime();
	for (i = 0; i < n; i++)
	{	c_part_out[i] = 0;
		for (j = 0; j < local_cols; j++){
			c_part_out[i] += a[i][j] * b[j];
    	}
	}
	cnt_out = my_malloc (my_rank, p * sizeof(int));
	disp_out = my_malloc (my_rank, p * sizeof(int));
	cnt_in = my_malloc (my_rank, p * sizeof(int));
	disp_in = my_malloc (my_rank, p * sizeof(int));
	createmixedxferarrays (my_rank, p, n, cnt_out, disp_out);
	createuniformparts ( my_rank, p, n, cnt_in, disp_in);
	c_part_in = (size_t*) my_malloc (my_rank, p*local_cols*sizeof(size_t));
	MPI_Alltoallv (c_part_out, cnt_out, disp_out, MPI_LONG, c_part_in, cnt_in, disp_in, MPI_LONG, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    c = (size_t*) my_malloc (my_rank, local_cols * sizeof(size_t));
	for (i = 0; i < local_cols; i++)
	{	c[i] = 0;
		for (j = 0; j < p; j++)
			c[i] += c_part_in[i + j*local_cols];
	}
	MPI_Barrier (MPI_COMM_WORLD);
	seconds += MPI_Wtime();
	MPI_Allreduce (&seconds, &max_seconds, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	size_t *c_final=(size_t*)malloc(sizeof(size_t)*n);
    replicate_block_vector_for_c(c,m,c_final,MPI_COMM_WORLD);
    //print_vector(c_final,n);
    
	writevector("result.txt",c_final,m);
	MPI_Finalize();
	return 0;
}