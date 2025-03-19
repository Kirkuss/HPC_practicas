#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

void addMatrices(int m, int n, double *A, double *B, double *C) {
    for (int i = 0; i < m * n; i++) {
        C[i] = A[i] + B[i];
    }
}

void multiplyMatrices(int m, int k, int n, double *A, double *B, double *C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i * n + j] = 0.0;
            for (int l = 0; l < k; l++) {
                C[i * n + j] += A[i * k + l] * B[l * n + j];
            }
        }
    }
}

void printMatrix(int m, int n, double *matrix) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", matrix[i * n + j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    int m, k, n;
    int rank, size;

    clock_t inicio, fin;
    double duration;

    m = atoi(argv[1]);
    k = atoi(argv[2]);
    n = atoi(argv[3]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4) {
        if (rank == 0) {
            printf("Usage: %s <m> <k> <n>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    int local_rows = m / size;
    int remainder = m % size;
    int *sendcounts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));
    
    double *A = malloc(m * k * sizeof(double));
    double *B = malloc(k * n * sizeof(double));
    double *final_E = malloc(m * n * sizeof(double));
    
    srand(time(NULL));
    for (int i = 0; i < m * k; i++) A[i] = (double)rand() / RAND_MAX;
    for (int i = 0; i < k * n; i++) B[i] = (double)rand() / RAND_MAX;
    
    // Configuración de la distribución de filas entre procesos
    int offset = 0;
    for (int i = 0; i < size; i++) {
        sendcounts[i] = (i < remainder) ? (local_rows + 1) * k : local_rows * k;
        displs[i] = offset;
        offset += sendcounts[i];
    }

    inicio = clock();

    // Difusión de la matriz B a todos los procesos
    MPI_Bcast(B, k * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Asignación de memoria para la porción de A que recibirá cada proceso
    double *local_A = malloc(sendcounts[rank] * sizeof(double));

    // Distribución de las filas de A entre los procesos
    MPI_Scatterv(A, sendcounts, displs, MPI_DOUBLE, local_A, sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int my_rows = sendcounts[rank] / k;  // Número de filas que cada proceso maneja

    // Matrices intermedias en cada proceso
    double *C = malloc(my_rows * n * sizeof(double));
    double *D = malloc(my_rows * n * sizeof(double));
    double *E = malloc(my_rows * n * sizeof(double));

    // C = A * B
    multiplyMatrices(my_rows, k, n, local_A, B, C);

    // D = C + B
    addMatrices(my_rows, n, C, B, D);

    // E = D * B
    multiplyMatrices(my_rows, k, n, D, B, E);

    // Configuración de la recepción de resultados en el proceso raíz
    int *recvcounts = malloc(size * sizeof(int));
    int *displsE = malloc(size * sizeof(int));

    offset = 0;
    for (int i = 0; i < size; i++) {
        recvcounts[i] = sendcounts[i] / k * n;
        displsE[i] = offset;
        offset += recvcounts[i];
    }

    // Recolección de los resultados en el proceso 0
    MPI_Gatherv(E, recvcounts[rank], MPI_DOUBLE, final_E, recvcounts, displsE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    fin = clock();
    duration = (double)(fin - inicio) / CLOCKS_PER_SEC;
    
    if (rank == 0) {
        printf("%d %d %d\n", m, n, k);
        printf("C=A*B, D=C+B, E=D*B: %2.1f segundos\n", duration);
        //printf("Matrix E (D * B):\n");
        //printMatrix(m, n, final_E);
    }
    
    free(sendcounts); free(displs); free(recvcounts); free(displsE);
    free(final_E); free(A); free(B); free(C); free(D); free(E); free(local_A);

    MPI_Finalize();

    return 0;
}
