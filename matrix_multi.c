/*
This program perform square matrix multiplication between two matrices (A & B).
A & B are generated with random value and will be store in text file.
This text file will be read as input file to generate matrix C.
Matrix C will be a matrix resultant and generated as output file.
*/

#include <stdio.h>
#include <time.h>
#include <sys/time.h> 
#include <stdlib.h>


//Generate matrix A
void GenInputMatrixA(int** matrix, int size){
	FILE *input;

	//creating file
	input = fopen("MatrixA.txt","w");

	if(input == NULL){
		printf("Unable to open file\n");
		printf("..::Failed creating file::..");
		return;
	}

	int row, col;
	//generate matrix element using random value
	for(row=0; row<size; row++){
		for (col = 0; col < size; col++)
		{
			matrix[row][col]=rand()%10;
		}
	}

	//print the element of matrix into a file
	for(row=0; row<size; row++){
		for (col = 0; col < size; col++)
		{
			fprintf(input, "%4d", matrix[row][col]);
		}
		fprintf(input, "\n");
	}

	fclose(input);

}

//Generate matrix B
void GenInputMatrixB(int** matrix, int size){
	FILE *input;

	//creating file
	input = fopen("MatrixB.txt","w");

	if(input == NULL){
		printf("Unable to open file\n");
		printf("..::Failed creating file::..");
		return;
	}

	int row, col;
	
	//generate matrix element using random value
	for(row=0; row<size; row++){
		for (col = 0; col < size; col++)
		{
			matrix[row][col]=rand()%10;
		}
	}

	//print the element of matrix into a file
	for(row=0; row<size; row++){
		for (col = 0; col < size; col++)
		{
			fprintf(input, "%4d", matrix[row][col]);
		}
		fprintf(input, "\n");
	}

	fclose(input);

}

//Read first input file
void ReadInputA(int** matrix, int size){
	FILE *input;
	//open and read the input file
	input = fopen("MatrixA.txt","r");
	
	int row, col;

	//scan and store the value of matrix in the input file
	for(row=0; row<size; row++){
		for (col = 0; col < size; col++)
		{
		   	if (!fscanf(input, "%4d", &matrix[row][col])) 
           	break;
		}
		
	}
}

//Read second input file
void ReadInputB(int** matrix, int size){
	FILE *input;
	//open and read the input file
	input = fopen("MatrixB.txt","r");
	
	int row, col;

	//scan and store the value of matrix in the input file
	for(row=0; row<size; row++){
		for (col = 0; col < size; col++)
		{
		   	if (!fscanf(input, "%4d", &matrix[row][col])) 
           	break;
		}
		
	}
}

void MatrixMultiplication(int** matA, int** matB, int** out, int size){
	int rowA, colA, rowB, colB, rowC, colC;
	FILE *output;
	//creating output file
	output = fopen("MatrixC.txt","w");

	rowC = rowB;
	colC = colB;

	//memory allocation for resultant of matrix
	out = (int**)malloc(size * sizeof(int*));
	for (colB = 0; colB < size; colB++)
	{
		out[colB] = (int*)malloc(size * sizeof(int));
	}

	//calculate matrix
	for (rowA = 0; rowA < size; rowA++){
		for (colB = 0; colB < size; colB++){
			out[rowA][colB] = 0;
			for (rowC = 0; rowC < size; rowC++){
				out[rowA][colB] += matA[rowA][rowC] * matB[rowC][colB]; 
			}
		}
	}
	
	//print the resultant of matrix into the output file
	for(rowC=0; rowC<size; rowC++){
		for (colC = 0; colC < size; colC++)
		{
		   	fprintf(output, "%4d", out[rowC][colC]);
		}
		fprintf(output, "\n");		
	}

	fclose(output);

	//free memory
	for(int rowC = 0; rowC < size; rowC++){
		free(out[rowC]);
	}
	free(out);

}

//free memory
void FreeMem(int** matrix, int size){
	for(int i = 0; i < size; i++){
		free(matrix[i]);
	}
	free(matrix);
}


int main(){
	//declare variable
	int size;
	int rowA, colA, rowB, colB;
	int **A, **B, **C;
	int rowC, colC;
	double start, tstart=0, tend=0, tcreat=0, tread=0, tmult=0;
	struct timeval tim;
	
	rowC = rowB;
	colC = colB;

	srand(time(NULL));
	
	//determine the ordo of matrices
	printf("Enter matrix size: ");
	scanf("%d", &size);

	gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
    tstart = start;

    //memory allocation for input matrix
	A = (int**)malloc(size * sizeof(int*));	
	for (colA = 0; colA < size; colA++)
	{
		A[colA] = (int*)malloc(size * sizeof(int));
	}

	B = (int**)malloc(size * sizeof(int*));
	for (colB = 0; colB < size; colB++)
	{
		B[colB] = (int*)malloc(size * sizeof(int));
	}
	
	//Generate input matrices
	gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
	GenInputMatrixA(A, size);
	GenInputMatrixB(B, size);
	gettimeofday(&tim, NULL);
    tcreat = tcreat + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
	
    //read input file
	gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
	ReadInputA(A, size);
	ReadInputB(B, size);
	gettimeofday(&tim, NULL);
    tread = tread + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
	
    //calculate the matrices
	gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
	MatrixMultiplication(A, B, C, size);
	gettimeofday(&tim, NULL);
    tmult = tmult + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
	
    //free memory
	FreeMem(A, size);
	FreeMem(B, size);

	gettimeofday(&tim, NULL);
    tend = tim.tv_sec+(tim.tv_usec/1000000.0);

	printf("Input file has been created: %.6lf\n", tcreat);
	printf("Input file has been readed: %.6lf\n", tread);
	printf("Output file created successfully: %.6lf \n", tmult);
	printf("Total time execution: %.6lf\n",tend - tstart);

	return 0;

}