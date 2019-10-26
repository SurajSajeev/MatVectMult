#include<stdlib.h>
#include <omp.h>
#include <stdio.h>
#include<time.h>
#include<math.h>
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
size_t* vector_alloc_0(size_t n){
size_t *buffer;
	if ((buffer = (size_t*)malloc ((size_t) n*sizeof(size_t))) == NULL)
	{	
		printf("No or less memory available");
	}
	for(int i=0;i<n;i++)
	buffer[i]=0;
	return (size_t*) buffer;
}
size_t **matrix_alloc_for_p(size_t r,size_t c){
int len=0; 
    size_t *ptr, **arr; 
    int count = 0,i,j; 
  
    len = sizeof(size_t) * c * r; 
    ptr =malloc(len*sizeof(size_t)); 
    arr=(size_t**)malloc(r*(sizeof(size_t*)));
    // ptr is now pointing to the first element in of 2D array 
   
  
    // for loop to point rows pointer to appropriate location in 2D array 
    for(i = 0; i < r; i++) 
        arr[i] = (ptr + c * i); 
  
    for (i = 0; i < r; i++) 
        for (j = 0; j < c; j++) 
            arr[i][j] = count++; // OR *(*(arr+i)+j) = ++count 
  
    return (arr); 
}
int *matrix_alloc(int r,int c){
int len=0; 
    int *ptr, **arr; 
    int count = 0,i,j; 
  
    len = sizeof(int *) * r + sizeof(int) * c * r; 
    arr = (int **)malloc(len*sizeof(int)); 
  
    // ptr is now pointing to the first element in of 2D array 
    ptr = (int *)(arr + r); 
  
    // for loop to point rows pointer to appropriate location in 2D array 
    for(i = 0; i < r; i++) 
        arr[i] = (ptr + c * i); 
  
    for (i = 0; i < r; i++) 
        for (j = 0; j < c; j++) 
            arr[i][j] = 0; // OR *(*(arr+i)+j) = ++count 
  
    
    return (int*)arr; 
}
void print_vector(size_t *vector,int n){
	for(int i=0;i<n;i++)
	printf("%lu ",vector[i]);
	printf("\n");
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
void work(size_t** inputMatrix, size_t* inputVector, int p,int numRow,int numCol ){
   
    size_t *oV = (size_t*)malloc(numRow*sizeof(size_t));
    for(int i=0;i<numRow;i++)
        oV[i]=0;
    int rootp=(int)sqrt((double)p);
    int num_block_rows=numRow/rootp;
    int num_block_cols=numCol/rootp;
     

    #pragma omp parallel num_threads(p)
    {
    size_t * tres=vector_alloc_0(numRow);   
    #pragma omp for
    for(int l=0;l<p;l++){
      int rowbn=l/rootp;
      int colbn=l%rootp;
      for(int i=rowbn*num_block_rows;i<(rowbn+1)*num_block_rows;i++)
      { 
        tres[i]=0;
        for(int j=colbn*num_block_cols;j<(colbn+1)*num_block_cols;j++){
        tres[i]+=inputMatrix[i][j]*inputVector[j];
      }
      }
    }
    
    #pragma omp critical
    {
        for(int r=0;r<numRow;r++)
        {
           oV[r]+=tres[r]; 
    }
    free(tres);
    }
    
}


//print_vector(oV,numRow);
}
int main(int argc,char ** argv){
size_t numrows,numcols;
numrows=atoi(argv[1]);
numcols=atoi(argv[2]);
size_t **a=matrix_alloc_for_p(numrows,numcols);


size_t *inputVector = vector_alloc(numcols);

     double t = 0.0;
        t -= omp_get_wtime();
        work(a, inputVector,4,numrows ,numcols);
        t += omp_get_wtime();
        printf("PT for %d\tt:%lf\n",4, t);
        t = 0.0;
}