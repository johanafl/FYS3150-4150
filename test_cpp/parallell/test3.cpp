#include <mpi.h>
#include <stdio.h>

// Run with:  mpiexec -n 8 ./test.out 
// -n is number of processors, so this runs with 8 processors.

int main(int argc, char** argv) {
    // Initialize the MPI environment
    double total = 0.0;
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // printf("%d\n",world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // printf("%d\n",world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    // double a = 7.0*world_rank;    
    // double total2 = 0.0;

    // MPI_Reduce(&a,&total2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    // if (world_rank==0) total=total2;

    // std::cout << "Im rank = " << world_rank << " and total = " << total << std::endl;

    // // Finalize the MPI environment.
    MPI_Finalize();

    // std::cout <<  total << std::endl;
    // std::cout << "Im rank = " << world_rank << " and total = " << total << std::endl;
}
