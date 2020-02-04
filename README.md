# Matrix Vector Multiplication
1 Introduction
The multiplication of the matrix and vector is often used in wide variety of applications,as it
has inseparable significance in the field of applied mathematics, statistics, physics, economics,
and engineering.

So, the computation of large matrices and vectors are usually time consuming as it has to
perform M ∗N multiplications and M ∗(N −1) additions respectively and only one operations
can be performed at a time.
So parallelization must be done in order to improve the overall time required to compute the
overall time required to solve the problem.
2 Ways to approach the problem
The problem can be solved in following ways:
1. Sequential approach
2. Parallel approach

The parallel approaches can further be classified as

  * Row partitioning multiplication
  * Column partitioning multiplication
  * checkerboard partitioning multiplication
  
## 1.1  Sequential Implementation

It is pretty straightforward that the sequential code can be implemented from the pseudocode
as follows:
```
1. procedure MAT_MULT (A, B, C)
2. begin
3. for i := 0 to n - 1 do
4. for j := 0 to n - 1 do
5. begin
6. C[j] := C[j] + A[i, j] x B[j];
7. end
9. endfor;
10. end MAT_MULT
```
## 1.2 Parallel Implementation
#### 3.1 Row wise partitioned multiplication
in this method, we will divide the whole matrix into submatrix by dividing it row wise and
perform the further calculations on different processes where the number of proceses is equal
to number of processors P. also each process will get N/P number of rows for an N ∗ N
matrix,each process will result in (N/P) ∗ 1 sized vector.

#### 3.2 Column wise partitioned multiplication
in this method, we will divide the whole matrix into submatrix by dividing it column wise
and perform the further calculations on different processes where the number of proceses is
equal to number of processors P. also each process will get N/P number of column for an
N ∗ N matrix,each process will result in (N) ∗ 1 sized partial result vector from which all of
them are needed to be added.

#### 3.3 Checkerboard partitioning multiplication
In this case we will divide the whole matrix into submatrix by dividing it into grids which is
equal to number of processors P,assuming that P is a perfect square,now we will divide the
submatrix of N/√
P ∗ N/√
P, and now we will multiply with its relative subvector and store
partial result and do the remaining task of reducing for the result.
