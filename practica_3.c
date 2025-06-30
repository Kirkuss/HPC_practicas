#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define M 500
#define K 500
#define N 500

void initialize_matrix(double* matrix, int rows, int cols) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i * cols + j] = (double)rand() / RAND_MAX;
        }
    }
}

void add_matrices(double* A, double* B, double* C, int rows, int cols) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            C[i * cols + j] = A[i * cols + j] + B[i * cols + j];
        }
    }
}

void multiply_matrices(double* A, double* B, double* C, int rows_a, int common, int cols_b) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows_a; i++) {
        for (int j = 0; j < cols_b; j++) {
            double sum = 0.0;
            for (int k = 0; k < common; k++) {
                sum += A[i * common + k] * B[k * cols_b + j];
            }
            C[i * cols_b + j] = sum;
        }
    }
}

void print_matrix(const char* name, double* matrix, int rows, int cols) {
    printf("Matrix %s:\n", name);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%0.2f ", matrix[i * cols + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char** argv) {
    if (argc != 4) {
        printf("Uso: %s <m> <k> <n>\n", argv[0]);
        return -1;
    }

    int m = atoi(argv[1]);
    int k = atoi(argv[2]);
    int n = atoi(argv[3]);

    double *A, *B, *C, *D, *E;

    A = (double*)malloc(m * k * sizeof(double));
    B = (double*)malloc(k * n * sizeof(double));
    C = (double*)malloc(m * n * sizeof(double));
    D = (double*)malloc(m * n * sizeof(double));
    E = (double*)malloc(m * n * sizeof(double));

    initialize_matrix(A, m, k);
    initialize_matrix(B, k, n);

    double start, end;

    // C = A + B
    start = omp_get_wtime();
    add_matrices(A, B, C, m, n);
    end = omp_get_wtime();
    printf("Tiempo para C = A + B: %f segunds\n", end - start);

    // D = C + B
    start = omp_get_wtime();
    add_matrices(C, B, D, m, n);
    end = omp_get_wtime();
    printf("Tiempo para D = C + B: %f segundos\n", end - start);

    // E = D * C
    start = omp_get_wtime();
    multiply_matrices(D, C, E, m, n, n);
    end = omp_get_wtime();
    printf("Tiempo para E = D * C: %f segundos\n", end - start);

    printf("\n");

    // print_matrix("A", A, M, K);
    // print_matrix("B", B, K, N);
    // print_matrix("C", C, M, N);
    // print_matrix("D", D, M, N);
    // print_matrix("E", E, M, N);

    free(A);
    free(B);
    free(C);
    free(D);
    free(E);

    return 0;
}
