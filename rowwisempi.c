#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<omp.h>

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))			// low_value= (in/p)
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)		// high value (i+1)n/p-1
#define BLOCK_SIZE(id,p,n) \
                     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)

size_t* vector_alloc(size_t n){
size_t *buffer;
	if ((buffer = (size_t*)malloc ((size_t) n*sizeof(size_t))) == NULL)
	{	
		printf("No or less memory available");
	}
	for(int i=0;i<n;i++)
	buffer[i]=1;
	return (size_t*) buffer;
}
void vector_free(void *v){
free(v);
}
void matrix_free(void **m){
free((void *)m);
}
size_t *matrix_alloc(size_t r,size_t c){
int len=0; 
    size_t *ptr, **arr; 
    int count = 0,i,j; 
  
    len = sizeof(int *) * r + sizeof(int) * c * r; 
    arr = (size_t **)malloc(len*sizeof(size_t)); 
  
    // ptr is now pointing to the first element in of 2D array 
    ptr = (size_t *)(arr + r); 
  
    // for loop to point rows pointer to appropriate location in 2D array 
    for(i = 0; i < r; i++) 
        arr[i] = (ptr + c * i); 
  
    for (i = 0; i < r; i++) 
        for (j = 0; j < c; j++) 
            arr[i][j] = 0; // OR *(*(arr+i)+j) = ++count 
  
    
    return (size_t*)arr; 
}

size_t *matrix_alloc_for_p(size_t r,size_t c){
int len=0; 
    size_t *ptr, **arr; 
    int count = 0,i,j; 
  
    len = sizeof(size_t *) * r + sizeof(size_t) * c * r; 
    arr = (size_t **)malloc(len*sizeof(size_t)); 
  
    // ptr is now pointing to the first element in of 2D array 
    ptr = (size_t *)(arr + r); 
  
    // for loop to point rows pointer to appropriate location in 2D array 
    for(i = 0; i < r; i++) 
        arr[i] = (ptr + c * i); 
  
    for (i = 0; i < r; i++) 
        for (j = 0; j < c; j++) 
            arr[i][j] = count++; // OR *(*(arr+i)+j) = ++count 
  
    
    return (size_t*)(arr+c); 
}

void print_matrix(size_t *arr,int r,int c){
int i,j;
for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) 
            printf("%lu ", *(arr+i*c+j)); 
	printf("\n");
	}
}
void print_vector(size_t *vector,int n){
	for(int i=0;i<n;i++)
	printf("%lu ",vector[i]);
	printf("\n");
}
size_t *matvectmul(size_t *a,size_t *b,int m,int n){
//size_t * c=vector_alloc(m);
 for (int i=0;i<m;i++){
	 	//c[i]=0;
        for (int j=0;j<n;j++){
            //c[i]+=( *(a+i*n+j)*b[j]);
        }
}
return b;
}

//STARTING TO WRITE THE MPI CODE HERE


void createxarrays(int id,int p, size_t n, int *count, int *disp){
     int i;
	// order in the process number. coun = no. of elements being sent, disp is the local offset of ith process's data
	count[0] = BLOCK_SIZE(0,p,n);
	disp[0] = 0;
	for (i = 1; i < p; i++)
	{	disp[i] = disp[i-1] + count[i-1];
		count[i] = BLOCK_SIZE(i,p,n);
	}
}

void createmixedxarrays(int id,int p, size_t n, int **count, int **disp){
     int i;

    *count = (int*) malloc (p*sizeof(int));
   *disp = (int*) malloc (p*sizeof(int));
    (*count)[0] = BLOCK_SIZE(id,p,n);
   (*disp)[0] = 0;
   for (i = 1; i < p; i++) {
      (*disp)[i] = (*disp)[i-1] + (*count)[i-1];
      (*count)[i] = BLOCK_SIZE(i,p,n);
   }
}
 
size_t* read_row_stripped_matrix(int m,int n,MPI_Comm comm,size_t * storage){
    int i,i1,my_rank,p;
    int *rptr;
    MPI_Status status;
    
    MPI_Comm_size(comm,&p);
    MPI_Comm_rank(comm,&my_rank);
    int local_rows=BLOCK_SIZE(my_rank,p,m);
    if(my_rank==p-1){
        size_t *a=(size_t *)matrix_alloc_for_p(m,n);
        if(p==1)
        return a;
        for(i=0;i<p-1;i++){
            int temp=0;
            storage=matrix_alloc(BLOCK_SIZE(i,p,m),n);
            for(i1=BLOCK_LOW(i,p,m)*n;i1<(BLOCK_HIGH(i,p,m)+1)*n;i1++)
            {storage[temp++]=a[i1];
            //printf("%lu\n",a[i1]);
            }
            MPI_Send(storage,(BLOCK_SIZE(i,p,m))*n,MPI_LONG,i,i,comm);
            
            
        }
        int temp=0;
        for(i1=BLOCK_LOW(p-1,p,m)*n;i1<(BLOCK_HIGH(p-1,p,m)+1)*n;i1++)
            {storage[temp++]=a[i1];
            //printf("%lu\n",a[i1]);
            }
            
    }
    else
    {   storage=matrix_alloc(BLOCK_SIZE(i,p,m),n);
        
        MPI_Recv(storage,local_rows*n,MPI_LONG,p-1,my_rank,comm,&status);
        //printf("For process %d with %d rows the matrix is like:\n\n",my_rank,local_rows);
        //print_matrix(storage,local_rows,n);
    }
    return storage;
}
size_t* read_vector_and_replicate(int n,size_t *v,MPI_Comm comm){
  
  int i, my_rank, p;
	MPI_Comm_rank (comm, &my_rank);
	MPI_Comm_size (comm, &p);
	if (my_rank == (p-1)) 
	{	v=vector_alloc(n);
	}
    if(p==1)
    return v;
	MPI_Bcast (v, n, MPI_LONG, p-1, comm);
    
    return v;
}


void replicate_block_vector ( size_t *ablock, int n, size_t *arep, MPI_Comm comm)
{	int my_rank, p;
	MPI_Comm_size (comm, &p);
	MPI_Comm_rank (comm, &my_rank);
	int *cnt = (int*) malloc (p*sizeof(int));
	int *disp = (int*) malloc (p*sizeof(int));
	createxarrays(my_rank, p, n, cnt, disp);				// deciding the order and initializing the 2 arrays
	MPI_Allgatherv (ablock, cnt[my_rank], MPI_LONG, arep, cnt, disp, MPI_LONG, comm);
	free (cnt);
	free (disp);
}


int main(int argc, char *argv[]){
 size_t **a;
 size_t *b;
 size_t *local_a;
 size_t *c;
 size_t *c_block;
  double max_seconds, seconds;
  int i, j, my_rank, m, n, p, local_rows, its, nb, size_of_type;   // nb-Elements in vector
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &p);
  if ( my_rank == (p-1))
	{	m=atoi(argv[1]);
        n=atoi(argv[2]);							// nb must be equal to n
	}
MPI_Bcast (&m, 1, MPI_INT, p-1, MPI_COMM_WORLD);
	MPI_Bcast (&n, 1, MPI_INT, p-1, MPI_COMM_WORLD);
	MPI_Bcast (&nb, 1, MPI_INT, p-1, MPI_COMM_WORLD);
    local_rows=BLOCK_SIZE(my_rank,p,n);
    c_block=(size_t*)malloc(local_rows*sizeof(size_t));
    c=(size_t*)malloc(m*sizeof(size_t));
	double t=0.0;

    local_a=read_row_stripped_matrix(m,n,MPI_COMM_WORLD,local_a);
    b=read_vector_and_replicate(n,b,MPI_COMM_WORLD);
     t -= omp_get_wtime();
    fflush(stdout);
    for(int l=0;l<local_rows;l++){
      c_block[l]=0;  
      for(int k=0;k<n;k++)
      c_block[l]+=*(local_a+l*n+k)*b[k];
    }
    replicate_block_vector(c_block,m,c,MPI_COMM_WORLD);
     t += omp_get_wtime();
    if(my_rank==0)
    {//printf("For process %d the result of vector will be\n\n",my_rank);
    //print_vector(c,m);}
    printf("the time taken is %f",t);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    }
}