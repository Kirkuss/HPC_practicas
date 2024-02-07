#include <stdio.h>
#include <time.h>

// Function to multiply two matrices A and B with loop ordering ijk
void multiplyMatrices_ijk(int m, int k, int n, double A[m][k], double B[k][n], double C[m][n]) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0.0;
            for (int l = 0; l < k; l++) {
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
}

// Function to print a matrix
void printMatrix(int m, int n, double matrix[m][n]) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main() {
    // Define the dimensions of matrices A, B, and C
    int m1 = 1000; // number of rows of A and C for the first case
    int k1 = 1000; // number of columns of A and rows of B
    int n1 = 1000; // number of columns of B and C

    // Define the matrices A, B, and C
    double A1[m1][k1];
    double B1[k1][n1];
    double C1[m1][n1] = {0};

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

    return 0;
}
