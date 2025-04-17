#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


void initialize_matrix(double* matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i * cols + j] = (double)rand() / RAND_MAX; // Inicializa con números aleatorios entre 0 y 1
        }
    }
}
    

void add_matrices(double* A, double* B, double* C, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            C[i * cols + j] = A[i * cols + j] + B[i * cols + j];
        }
    }
}

void multiply_matrices(double* A, double* B, double* C, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            C[i * cols + j] = A[i * cols + j] * B[i * cols + j];
        }
    }
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
    double start_time, end_time;

    if (rank == 0) {
        A = (double*)malloc(m * k * sizeof(double));
        B = (double*)malloc(k * n * sizeof(double));
        C = (double*)malloc(m * n * sizeof(double));
        D = (double*)malloc(m * n * sizeof(double));
        E = (double*)malloc(m * n * sizeof(double));

        // Inicializar matrices A y B con algunos valores
        initialize_matrix(A, m, k);  // Por ejemplo, inicializar A con 1.0
        initialize_matrix(B, k, n);  // Por ejemplo, inicializar B con 2.0
    }

    // Scatter: Distribuir filas de A y B a todos los procesos
    double *local_A = (double*)malloc((m / size) * k * sizeof(double));
    double *local_B = (double*)malloc(k * n * sizeof(double));
    double *local_C = (double*)malloc((m / size) * n * sizeof(double));
    double *local_D = (double*)malloc((m / size) * n * sizeof(double));
    double *local_E = (double*)malloc((m / size) * n * sizeof(double));

    MPI_Scatter(A, (m / size) * k, MPI_DOUBLE, local_A, (m / size) * k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, k * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Medir tiempo para C = A + B
    start_time = MPI_Wtime();
    add_matrices(local_A, B, local_C, m / size, n);
    end_time = MPI_Wtime();
    if (rank == 0) {
        printf("Tiempo para C = A + B: %f segundos\n", end_time - start_time);
    }

    // Gather: Recolectar resultados de C en el proceso raíz
    MPI_Gather(local_C, (m / size) * n, MPI_DOUBLE, C, (m / size) * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Scatter: Distribuir filas de C a todos los procesos
    MPI_Scatter(C, (m / size) * n, MPI_DOUBLE, local_C, (m / size) * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Medir tiempo para D = C + B
    start_time = MPI_Wtime();
    add_matrices(local_C, B, local_D, m / size, n);
    end_time = MPI_Wtime();
    if (rank == 0) {
        printf("Tiempo para D = C + B: %f segundos\n", end_time - start_time);
    }

    // Gather: Recolectar resultados de D en el proceso raíz
    MPI_Gather(local_D, (m / size) * n, MPI_DOUBLE, D, (m / size) * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Scatter: Distribuir filas de D y C a todos los procesos
    MPI_Scatter(D, (m / size) * n, MPI_DOUBLE, local_D, (m / size) * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(C, (m / size) * n, MPI_DOUBLE, local_C, (m / size) * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Medir tiempo para E = D * C
    start_time = MPI_Wtime();
    multiply_matrices(local_D, local_C, local_E, m / size, n);
    end_time = MPI_Wtime();
    if (rank == 0) {
        printf("Tiempo para E = D * C: %f segundos\n", end_time - start_time);
    }

    // Gather: Recolectar resultados de E en el proceso raíz
    MPI_Gather(local_E, (m / size) * n, MPI_DOUBLE, E, (m / size) * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Sincronizar todos los procesos
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        // Imprimir matriz E
        // printf("Matriz E:\n");
        // for (int i = 0; i < m; i++) {
        //    for (int j = 0; j < n; j++) {
        //        printf("%f ", E[i * n + j]);
        //    }
        //    printf("\n");
        // }

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
