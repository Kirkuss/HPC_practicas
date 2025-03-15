#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

// Function to add two matrices A and B
void addMatrices(int m, int n, double **A, double **B, double **C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

// Function to multiply two matrices A and B with loop ordering ijk
void multiplyMatrices_ijk(int m, int k, int n, double **A, double **B, double **C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0.0; 
            for (int l = 0; l < k; l++) {
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}

void initializeResultsMatrix(int m, int n, double **matrix) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = 0.0;
        }
    }
}

void printMatrix(int m, int n, double **matrix) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    int m1, k1, n1;
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4) {
        if (rank == 0) {
            printf("Usage: %s <m> <k> <n>\n", argv[0]);
        }
        MPI_Finalize();
        exit(1);
    } else {
        m1 = atoi(argv[1]);
        k1 = atoi(argv[2]);
        n1 = atoi(argv[3]);
    }

    double **A1 = (double **)malloc(m1 * sizeof(double *));
    double **B1 = (double **)malloc(k1 * sizeof(double *));
    double **C1 = (double **)malloc(m1 * sizeof(double *));
    double **D1 = (double **)malloc(m1 * sizeof(double *));
    double **E1 = (double **)malloc(m1 * sizeof(double *));
    
    for (int i = 0; i < m1; i++) {
        A1[i] = (double *)malloc(k1 * sizeof(double));
        C1[i] = (double *)malloc(n1 * sizeof(double));
        D1[i] = (double *)malloc(n1 * sizeof(double));
        E1[i] = (double *)malloc(n1 * sizeof(double));
    }
    for (int i = 0; i < k1; i++) {
        B1[i] = (double *)malloc(n1 * sizeof(double));
    }
    
    srand(time(NULL) + rank); 
    for (int i = 0; i < m1; i++) {
        for (int j = 0; j < k1; j++) {
            A1[i][j] = (double)rand() / RAND_MAX;
        }
    }
    for (int i = 0; i < k1; i++) {
        for (int j = 0; j < n1; j++) {
            B1[i][j] = (double)rand() / RAND_MAX;
        }
    }

    addMatrices(m1, n1, A1, B1, C1);
    if (rank == 0) {
        printf("Matrix C (A + B):\n");
        printMatrix(m1, n1, C1);
    }

    addMatrices(m1, n1, C1, B1, D1);
    if (rank == 0) {
        printf("Matrix D (C + B):\n");
        printMatrix(m1, n1, D1);
    }

    multiplyMatrices_ijk(m1, k1, n1, D1, C1, E1);
    if (rank == 0) {
        printf("Matrix E (D * C):\n");
        printMatrix(m1, n1, E1);
    }

    for (int i = 0; i < m1; i++) {
        free(A1[i]);
        free(C1[i]);
        free(D1[i]);
        free(E1[i]);
    }
    for (int i = 0; i < k1; i++) {
        free(B1[i]);
    }
    free(A1);
    free(B1);
    free(C1);
    free(D1);
    free(E1);

    MPI_Finalize();
    return 0;
}