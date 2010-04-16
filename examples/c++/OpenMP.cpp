#include <omp.h>
#include <iostream>
#include <sstream>
#include <string>
int main()
{
    std::cout << "number of threads: " << omp_get_max_threads() << "(" << omp_get_num_procs() << ")" << std::endl;
    omp_set_num_threads(2);
#pragma omp parallel
    {
        // Code inside this region runs in parallel.
        std::string tmpString("Hello from thread: ");
        std::stringstream stringStream;
        stringStream << omp_get_thread_num();
        tmpString += stringStream.str() + "\n";
        std::cout << tmpString;
    }
    return 0;

}