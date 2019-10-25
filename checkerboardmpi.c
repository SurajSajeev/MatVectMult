#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))			// low_value= (in/p)
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)		// high value (i+1)n/p-1
#define BLOCK_SIZE(id,p,n) \
                     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)

int* vector_alloc(int n){
int *buffer;
	if ((buffer = (int*)malloc ((int) n*sizeof(int))) == NULL)
	{	
		printf("No or less memory available");
	}
	for(int i=0;i<n;i++)
	buffer[i]=1;
	return (int*) buffer;
}
void vector_free(void *v){
free(v);
}
void matrix_free(void **m){
free((void *)m);
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

int *matrix_alloc_for_p(int r,int c){
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
            arr[i][j] = count++; // OR *(*(arr+i)+j) = ++count 
  
    
    return (int*)(arr+c); 
}

void print_matrix(int *arr,int r,int c){
int i,j;
for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) 
            printf("%d ", *(arr+i*c+j)); 
	printf("\n");
	}
}
void print_vector(int *vector,int n){
	for(int i=0;i<n;i++)
	printf("%d ",vector[i]);
	printf("\n");
}
int *matvectmul(int *a,int *b,int m,int n){
//int * c=vector_alloc(m);
 for (int i=0;i<m;i++){
	 	//c[i]=0;
        for (int j=0;j<n;j++){
            //c[i]+=( *(a+i*n+j)*b[j]);
        }
}
return b;
}


int main(int argc,char **argv){
MPI_Status status;
MPI_Comm row_comm,col_comm;
MPI_Group MPI_GROUP_WORLD;

int Numproc,MyRank,root_p,Root=0;
int irow,icol,iproc,jproc,index,Proc_Id;
int NoofRows,NoofCols,NoofRows_Bloc,NoofCols_Bloc;
int Bloc_MatrixSize,Bloc_VectorSize,VectorSize;
int Local_Index, Global_Row_Index, Global_Col_Index;
int **Matrix ,*Matrix_Array,*Bloc_Matrix,*Vector,*Bloc_Vector;
int *FinalResult, *MyResult,*FinalVector;

int *ranks,colsize,colrank,rowsize,rowrank,myrow;

MPI_Init(&argc,&argv);
MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP_WORLD);
MPI_Comm_rank(MPI_COMM_WORLD,&MyRank);
MPI_Comm_size(MPI_COMM_WORLD,&Numproc);
root_p=sqrt((double)Numproc);
if(Numproc!=(root_p*root_p))
{
    if(MyRank==0)
    printf("The number of processors should be a perfect square\n");
MPI_Finalize();
exit(-1);
}
NoofRows=20;
NoofCols=20;
if(MyRank==0){
if(NoofRows%root_p!=0 || NoofCols%root_p!=0){
    printf("The number of rows and columns sould be a multiple of root_p\n");
    MPI_Finalize();
    exit(-1);
}
Matrix_Array=matrix_alloc_for_p(NoofRows,NoofCols);
Matrix=malloc(NoofRows*sizeof(void*));
int *te=Matrix_Array;
Matrix[0]=te;
for(int i=1;i<NoofRows;i++){
te+=NoofCols;
Matrix[i]=te;
}
MPI_Bcast(&NoofRows,1,MPI_INT,0,MPI_COMM_WORLD);
MPI_Bcast(&NoofCols,1,MPI_INT,0,MPI_COMM_WORLD);
}
NoofRows_Bloc=NoofRows/root_p;
NoofCols_Bloc=NoofCols/root_p;
if(MyRank==0){
Matrix_Array=vector_alloc(NoofRows*NoofCols);
Local_Index=0;
for(iproc = 0; iproc < root_p; iproc++){
         for(jproc = 0; jproc < root_p; jproc++){
   	     Proc_Id = iproc * root_p + jproc;
	     for(irow = 0; irow < NoofRows_Bloc; irow++){
		 Global_Row_Index = iproc * NoofRows_Bloc + irow;
	         for (icol = 0; icol < NoofCols_Bloc; icol++){
		      Global_Col_Index = jproc * NoofCols_Bloc + icol;
	              Matrix_Array[Local_Index++] = 
			       Matrix[Global_Row_Index][Global_Col_Index];
                   //printf("%d\n",Matrix[Global_Row_Index][Global_Col_Index]);
	         }
	     }
	 }
       }
}
  Vector=vector_alloc(NoofCols);
  VectorSize=NoofRows;
  Bloc_VectorSize = VectorSize / root_p;
  Bloc_MatrixSize = NoofRows_Bloc * NoofCols_Bloc;
  Bloc_Matrix = (int *) malloc (Bloc_MatrixSize * sizeof(int));

  MPI_Scatter(Matrix_Array, Bloc_MatrixSize, MPI_INT, Bloc_Matrix, 
	      Bloc_MatrixSize, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); 
  
 /* Creating groups of procesors row wise */
  myrow=MyRank/root_p;
  MPI_Comm_split(MPI_COMM_WORLD,myrow,MyRank,&row_comm);
  MPI_Comm_size(row_comm,&rowsize);
  MPI_Comm_rank(row_comm,&rowrank);

  /* Creating groups of procesors column wise */
  myrow=MyRank%root_p;
  MPI_Comm_split(MPI_COMM_WORLD,myrow,MyRank,&col_comm);
  MPI_Comm_size(col_comm,&colsize);
  MPI_Comm_rank(col_comm,&colrank);
  Bloc_Vector = (int*) malloc(Bloc_VectorSize * sizeof(int));
  if(MyRank/root_p == 0){
     MPI_Scatter(Vector, Bloc_VectorSize, MPI_INT, Bloc_Vector, 
		 Bloc_VectorSize, MPI_INT, 0,row_comm);
  }

  MPI_Bcast(Bloc_Vector, Bloc_VectorSize, MPI_INT, 0, col_comm);

  MyResult   = (int*) malloc(NoofRows_Bloc * sizeof(int));
  index = 0;
  for(irow=0; irow < NoofRows_Bloc; irow++){
      MyResult[irow]=0;
      for(icol=0;icol< NoofCols_Bloc; icol++){
          MyResult[irow] += Bloc_Matrix[index++] * Bloc_Vector[icol];
      }
  }
    if(MyRank == Root) 
      FinalResult = (int *)malloc(NoofRows_Bloc*Numproc*sizeof(int));

   MPI_Gather (MyResult, NoofRows_Bloc, MPI_INT, FinalResult, 
	       NoofRows_Bloc, MPI_INT, 0, MPI_COMM_WORLD); 
   if(MyRank == 0){
      FinalVector = (int *) malloc(NoofRows * sizeof(int));
      index = 0;
      for(iproc=0; iproc<root_p; iproc++){
	  for(irow=0; irow<NoofRows_Bloc; irow++){
	      FinalVector[index]  = 0;
	      for(jproc=0; jproc<root_p; jproc++){
		  FinalVector[index] += 
		    FinalResult[iproc*root_p*NoofRows_Bloc + 
				jproc*NoofRows_Bloc +irow];

	      }
              printf("%d ",FinalVector[index]);
	      index++;
	  }
      }
   }


   /* Free the groups formed */
   MPI_Comm_free(&row_comm);
   MPI_Comm_free(&col_comm);


   MPI_Finalize(); 
}