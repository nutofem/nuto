//
// Created by phuschke on 2/10/17.
//

#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[])
{

    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "Hello from process: \t" << rank << "/" << size << std::endl;

    MPI_Finalize();
}
