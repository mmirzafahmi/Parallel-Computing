#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include "mpi.h"

int main(int argc, char **argv){
	//declare MPI variable
	int rank, size, start, finish;
	int i, j, k;
	int ordo = 7000;
	int A[ordo][ordo], B[ordo][ordo], C[ordo][ordo];
	double tstart, tend;
	FILE *MatrixA, *MatrixB, *MatrixC;
	//declare random seeder
	srand(time(NULL));
	//MPI initialization
	MPI_Init(&argc, &argv);
	MPI_Status stat;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	start = rank * ordo/size;
	finish = (rank+1) * ordo/size;

	
		
	if(rank == 0){
		tstart = MPI_Wtime();

		printf("Generate Matrix ...\n");
		

		MatrixA = fopen("MatrixA.txt","w");
		for(i = 0; i<ordo; i++){
			for(j = 0; j<ordo; j++){
				A[i][j] = rand()%10;
			}
		}
		
		for(i = 0; i<ordo; i++){
			for(j = 0; j<ordo; j++){
				fprintf(MatrixA, "%4d",A[i][j]);
			}
			fprintf(MatrixA, "\n");
		}
		
		MatrixB = fopen("MatrixB.txt","w");
		for(i = 0; i<ordo; i++){
			for(j = 0; j<ordo; j++){
				B[i][j] = rand()%10;
			}
		}
		
		for(i = 0; i<ordo; i++){
			for(j = 0; j<ordo; j++){
				fprintf(MatrixB, "%4d",B[i][j]);
			}
			fprintf(MatrixB, "\n");
		}

		fclose(MatrixA);
		fclose(MatrixB);
		
		printf("Generate Matrix Success...\n");


		
	}
		if(MPI_Bcast(B, ordo*ordo, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS){
			perror("Broadcast error");
			exit(1);
		}

		if(MPI_Scatter(A, ordo*rank, MPI_INT, A[start], ordo*rank, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS){
			perror("Scatter error");
			exit(1);
		} 

		for(i = start; i<finish; i++)
			for(j = 0; j<ordo; j++){
				C[i][j] = 0;
				for(k = 0; k<ordo; k++)
					C[i][j] = C[i][j] + A[i][k]*B[k][j];
				
			}
		

		MPI_Gather(C[start], ordo*rank, MPI_INT, C, ordo*rank, MPI_INT, 0, MPI_COMM_WORLD);	
	
			
	
	if(rank == 0){
		
		MatrixC = fopen("MatrixC.txt","w");
		for(i = 0; i<ordo; i++){
			for(j = 0; j<ordo; j++){
				fprintf(MatrixC, "%4d", C[i][j]);
			}	
			fprintf(MatrixC, "\n");
		}
		fclose(MatrixC);
		
		tend = MPI_Wtime();
		printf("Total Execution: %.6lf\n", tend - tstart);
		MPI_Finalize();
	
	}
	//tend = MPI_Wtime();
	//printf("Total execution: %.6lf\n", tend - tstart);
	//MPI_Finalize();
	return 0;

}
 
