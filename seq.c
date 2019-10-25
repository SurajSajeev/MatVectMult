#include<stdio.h>
#include<stdlib.h>
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
size_t **matrix_alloc(size_t r,size_t c){
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
            arr[i][j] = ++count; // OR *(*(arr+i)+j) = ++count 
  
    
    return (size_t**)arr; 
}
void print_matrix(size_t **arr,int r,int c){
int i,j;
for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) 
            printf("%lu ", arr[i][j]); 
	printf("\n");
	}
}
void print_vector(size_t *vector,int n){
	for(int i=0;i<n;i++)
	printf("%lu ",vector[i]);
	printf("\n");
}
size_t *matvectmul(size_t **a,size_t *b,int m,int n){
size_t * c=vector_alloc(m);
 for (int i=0;i<m;i++){
	 	c[i]=0;
        for (int j=0;j<n;j++){
            c[i]+=( a[i][j]*b[j]);
        }
}
return c;
}
int main(){
	size_t **a=matrix_alloc(10,10);
	size_t *b=vector_alloc(10);
	print_matrix(a,10,10);
	print_vector(b,10);
	size_t* c=matvectmul(a,b,10,10);
	print_vector(c,10);
	vector_free((void*)b);
	matrix_free((void **)a);
	
}