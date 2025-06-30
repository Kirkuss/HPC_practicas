#include <omp.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


void add_matrices(double* A, double* B, double* C, int rows, int cols) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            C[i * cols + j] = A[i * cols + j] + B[i * cols + j];
        }
    }
}

void multiply_matrices(double* A, double* B, double* C, int rows_a, int rows_b, int cols) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows_a; i++) {
        for (int j = 0; j < cols; j++) {
            double sum = 0.0;
            for (int k = 0; k < rows_b; k++) {
                sum += A[i * rows_b + k] * B[k * cols + j];
            }
            C[i * cols + j] = sum;
        }
    }
}

void initialize_matrix(double* matrix, int rows, int cols) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i * cols + j] = (double)rand() / RAND_MAX;
        }
    }
}

void transpose_matrix(double* input, double* output, int rows, int cols) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            output[j * rows + i] = input[i * cols + j];
        }
    }
}

void print_matrix(const char* name, double* matrix, int rows, int cols) {
    printf("Matriz %s:\n", name);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%0.2f ", matrix[i * cols + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4) {
        if (rank == 0) {
            printf("Uso: %s <m> <k> <n>\n", argv[0]);
        }
        MPI_Finalize();
        return -1;
    }

    int m = atoi(argv[1]);
    int k = atoi(argv[2]);
    int n = atoi(argv[3]);

    double *A, *B, *C, *D, *E;
    double *local_A, *local_B, *local_C, *local_D, *local_E;

    double start_time_c, end_time_c, start_time_d, end_time_d, start_time_e, end_time_e;

    local_A = (double*)malloc((m / size) * k * sizeof(double));
    local_B = (double*)malloc((k / size) * n * sizeof(double));
    local_C = (double*)malloc(m * n * sizeof(double));
    local_D = (double*)malloc((m / size) * n * sizeof(double));
    local_E = (double*)malloc((m / size) * n * sizeof(double));

    if (rank == 0) {
        A = (double*)malloc(m * k * sizeof(double));
        B = (double*)malloc(k * n * sizeof(double));
        C = (double*)malloc(m * n * sizeof(double));
        D = (double*)malloc(m * n * sizeof(double));
        E = (double*)malloc(m * n * sizeof(double));

        // Inicializar matrices A y B con algunos valores
        initialize_matrix(A, m, k);
        initialize_matrix(B, k, n);
    }
    else
    {
        C = local_C;
    }

    MPI_Scatter(A, (m / size) * k, MPI_DOUBLE, local_A, (m / size) * k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(B, (k / size) * n, MPI_DOUBLE, local_B, (k / size) * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Medir tiempo para C = A + B
    start_time_c = MPI_Wtime();
    add_matrices(local_A, local_B, local_C, m / size, n);
    end_time_c = MPI_Wtime();

    // Medir tiempo para D = C + B
    start_time_d = MPI_Wtime();
    add_matrices(local_C, local_B, local_D, m / size, n);
    end_time_d = MPI_Wtime();

    MPI_Gather(local_C, (m / size) * n, MPI_DOUBLE, C, (m / size) * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(C, m * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Medir tiempo para E = D * C
    start_time_e = MPI_Wtime();
    multiply_matrices(local_D, local_C, local_E, m / size, n, n);
    end_time_e = MPI_Wtime();

    // Gather: Recolectar resultados
    MPI_Gather(local_A, (m / size) * k, MPI_DOUBLE, A, (m / size) * k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(local_D, (m / size) * n, MPI_DOUBLE, D, (m / size) * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(local_E, (m / size) * n, MPI_DOUBLE, E, (m / size) * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Sincronizar todos los procesos
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Tiempo para C = A + B: %f segundos\n", end_time_c - start_time_c);
        printf("Tiempo para D = C + B: %f segundos\n", end_time_d - start_time_d);
        printf("Tiempo para E = D * C: %f segundos\n", end_time_e - start_time_e);

        printf("\n");

        //print_matrix("A", A, m, k);
        //print_matrix("B", B, k, n);
        //print_matrix("C = A + B", C, m, n);
        //print_matrix("D = C + B", D, m, n);
        //print_matrix("E = D * C", E, m, n);

        free(A);
        free(B);
        free(C);
        free(D);
        free(E);
    }

    free(local_A);
    free(local_B);
    free(local_C);
    free(local_D);
    free(local_E);

    MPI_Finalize();
    return 0;
}
