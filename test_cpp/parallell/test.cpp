#include <omp.h>
#include <cstdio>
#include <iostream>

int main (int argc, char *argv[])
{
    int th_id, nthreads;

    #pragma omp parallel private(th_id) shared(nthreads)
    {
        th_id = omp_get_thread_num();
        printf("Hello World from thread %d\n", th_id);
        #pragma omp barrier
        if ( th_id == 0 ) 
        {
            nthreads = omp_get_num_threads();
            printf("There are %d threads\n",nthreads);
        }
    }

    // int id = omp_get_thread_num();
    // #pragma omp parallel private(id)
    // {
    //     std::cout << "My thread num" << id << std::endl; 
    // }

    #pragma omp parallel 
    {
        #pragma omp master
        {
            int id = omp_get_thread_num();
            std::cout << "My thread num" << id << std::endl; 
        } 
    }

    return 0;
}
