#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

// Function to multiply two matrices A and B with loop ordering ijk
void multiplyMatrices_ijk(int m, int k, int n, double **A, double **B, double **C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int l = 0; l < k; l++) {
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}

// Function to multiply two matrices A and B with loop ordering ijk
void multiplyMatrices_jki(int m, int k, int n, double **A, double **B, double **C) {
    for (int j = 0; j < n; j++) {
        for (int l = 0; l < k; l++) {
            for (int i = 0; i < m; m++) {
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}

// Function to multiply two matrices A and B with loop ordering ijk
void multiplyMatrices_kji(int m, int k, int n, double **A, double **B, double **C) {
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}

void initializeResultsMatrix(int m, int n, double **matrix)  {
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = 0.0;
        }
    }
}

bool checkResults(int m, int n, double **matrix, double **check_matrix)  {
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (matrix[i][j] != check_matrix[i][j]){
              return false;
            }
        }
    }
    return true;
}

// Function to print a matrix
void printMatrix(int m, int n, double **matrix) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(int argc,char *argv[]) {
    int m1, k1, n1;

    if(argc != 4){
      exit(1);
    }else{
      m1 = atoi(argv[1]);
      k1 = atoi(argv[2]);
      n1 = atoi(argv[3]);
    }
    // Define the dimensions of matrices A, B, and C


    // Allocate memory for matrices A, B, and C
    double **A1 = (double **)malloc(m1 * sizeof(double *));
    double **B1 = (double **)malloc(k1 * sizeof(double *));
    double **C1 = (double **)malloc(m1 * sizeof(double *));
    
    for (int i = 0; i < m1; i++) {
        A1[i] = (double *)malloc(k1 * sizeof(double));
        C1[i] = (double *)malloc(n1 * sizeof(double));
    }
    for (int i = 0; i < k1; i++) {
        B1[i] = (double *)malloc(n1 * sizeof(double));
    }
    
    initializeResultsMatrix(m1, n1, C1);

    // Initialize matrices A and B with random values for the first case
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

    // Define variables to measure execution time
    clock_t start, end;
    double cpu_time_used;

    // Measure execution time for matrix multiplication with loop ordering ijk for the first case
    start = clock();
    multiplyMatrices_ijk(m1, k1, n1, A1, B1, C1);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Execution time for matrix multiplication (ijk) with matrix size 1000x1000: %f seconds\n", cpu_time_used);

    initializeResultsMatrix(m1, n1, C1);
    
    start = clock();
    multiplyMatrices_jki(m1, k1, n1, A1, B1, C1);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Execution time for matrix multiplication (jki) with matrix size 1000x1000: %f seconds\n", cpu_time_used);
    
    initializeResultsMatrix(m1, n1, C1);
    
    start = clock();
    multiplyMatrices_kji(m1, k1, n1, A1, B1, C1);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Execution time for matrix multiplication (kji) with matrix size 1000x1000: %f seconds\n", cpu_time_used);

    // Free allocated memory
    for (int i = 0; i < m1; i++) {
        free(A1[i]);
        free(C1[i]);
    }
    for (int i = 0; i < k1; i++) {
        free(B1[i]);
    }
    free(A1);
    free(B1);
    free(C1);

    return 0;
}
