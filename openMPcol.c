#include<stdlib.h>
#include <omp.h>
#include <stdio.h>
#include<time.h>

int numRow = 10000;
int numCol = 10000;
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
size_t* vector_alloc_for_0(size_t n){
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

void print_vector(size_t *vector,int n){
	for(int i=0;i<n;i++)
	printf("%lu ",vector[i]);
	printf("\n");
}
void work(size_t** inputMatrix, size_t* inputVector, int p){
    size_t *oV = (size_t*)malloc(numRow*sizeof(size_t));
    for(int i=0;i<numRow;i++)
        oV[i]=0;
        
    #pragma omp parallel num_threads(p)
    {   size_t *pV = (size_t*)malloc(numRow*sizeof(size_t));
        for(int i=0;i<numRow;i++)
        pV[i]=0;
        int my_work = 0;
        for(int r = 0; r < numRow; r++){
            
            #pragma omp for
            for(int c = 0; c < numCol; c++)
                pV[r] += inputMatrix[r][c] * inputVector[c];
        }
    #pragma omp critical
    {
        for(int r=0;r<numRow;r++)
        {
            oV[r]+=pV[r];}
    }
    }
    //print_vector(oV,numRow);
}

int main(){
    struct timespec begin, end;
    double t=0.0;
    
    size_t **inputMatrix = matrix_alloc_for_p(numRow,numCol);
   
    
    size_t *inputVector = vector_alloc(numCol);
     t = 0.0;
    for(int i=1; i<=8; i++){
        t -= omp_get_wtime();
        work(inputMatrix, inputVector, i);
        t += omp_get_wtime();
        printf("PT for %d\tt:%lf\n",i, t);
        t = 0.0;
    }
   
}