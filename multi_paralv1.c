#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "mpi.h"
#define ORDO 2000
#define MASTER 0

//declare matrix variable
int A[ORDO][ORDO], B[ORDO][ORDO], C[ORDO][ORDO];

int main(int argc, char **argv){
	//declare MPI variable
	int rank, size, slave, sendrow, exrow, row, dest, source, balancer;
	int i, j, k, l, m;
	FILE *MatrixA, *MatrixB, *MatrixC;
	//declare random seeder
	srand(time(NULL));
	//MPI initialization
	MPI_Init(&argc, &argv);
	MPI_Status stat;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	slave = size - 1;	

	if(rank == MASTER){
		double start = MPI_Wtime();

		printf("Generate Matrix ...\n");
		MatrixA = fopen("MatrixA.txt","w");
		for(i = 0; i<ORDO; i++){
			for(j = 0; j<ORDO; j++){
				A[i][j] = rand()%10;
			}
		}

		for(i = 0; i<ORDO; i++){
			for(j = 0; j<ORDO; j++){
				fprintf(MatrixA, "%4d",A[i][j]);
			}
			fprintf(MatrixA, "\n");
		}

		MatrixB = fopen("MatrixB.txt","w");
		for(i = 0; i<ORDO; i++){
			for(j = 0; j<ORDO; j++){
				B[i][j] = rand()%10;
			}
		}

		for(i = 0; i<ORDO; i++){
			for(j = 0; j<ORDO; j++){
				fprintf(MatrixB, "%4d",B[i][j]);
			}
			fprintf(MatrixB, "\n");
		}

		fclose(MatrixA);
		fclose(MatrixB);

		MatrixC = fopen("MatrixC.txt","w");

		sendrow = ORDO/slave;
		exrow = ORDO%slave;
		balancer = 0;

		for(dest = 1; dest <= slave; dest++){
			row = (dest <= exrow) ? sendrow+1:sendrow;
			printf("Send %d row to task %d \n",row,dest);
			MPI_Send(&balancer, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
			MPI_Send(&row, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
			MPI_Send(&A[balancer][0], row*ORDO, MPI_INT, dest, 1, MPI_COMM_WORLD);
			MPI_Send(&B, ORDO*ORDO, MPI_INT, dest, 1, MPI_COMM_WORLD);			
			balancer += row;
		}

		for(source = 1; source <= slave; source++){
			MPI_Recv(&balancer, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &stat);
			MPI_Recv(&row, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &stat);
			MPI_Recv(&C[balancer][0], row*ORDO, MPI_INT, source, 2, MPI_COMM_WORLD, &stat);
			printf("Receive from task %d \n",source);
		}
		
		//printf("Print the result ....\n");
		for(l = 0; l<ORDO; l++){
			for(m = 0; m<ORDO; m++){
				fprintf(MatrixC, "%4d", C[l][m]);
			}	
			fprintf(MatrixC, "\n");
		}
		fclose(MatrixC);

		double finish = MPI_Wtime();
		printf("Total Execution: %.6lf\n", finish - start);
	}

	else{
		MPI_Recv(&balancer, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);
		MPI_Recv(&row, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);
		MPI_Recv(&A, row*ORDO, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);
		MPI_Recv(&B, ORDO*ORDO, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);

		for(i = 0; i<ORDO; i++){
			for(j = 0; j<ORDO; j++){
				C[i][j] = 0;
				for(k = 0; k<ORDO; k++){
					C[i][j] = C[i][j] + A[i][k]*B[k][j];
				}
			}
		}

		
		MPI_Send(&balancer, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
		MPI_Send(&row, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
		MPI_Send(&C, ORDO*row, MPI_INT, 0, 2, MPI_COMM_WORLD);		

	}	

	MPI_Finalize();

}
 
