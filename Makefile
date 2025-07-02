CC_MPI=mpicc
CC_OMP=gcc

mpi_g10: practica_2.c
        $(CC_MPI) practica_2.c -o mpi_grupo10.exe

omp_g10: practica_3.c
        $(CC_OMP) -fopenmp practica_3.c -o omp_grupo10.exe

hybrid_g10: practica_4.c
        $(CC_MPI) -fopenmp practica_4.c -o hybrid_grupo10.exe

clean:
        rm *.exe
