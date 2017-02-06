#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <eigen3/Eigen/Core>
#include "math/EigenCompanion.h"

#include <ctime>
# ifdef _OPENMP
#include <omp.h>
# endif

void ComputeNoParallel(int rNumParallelLoops, int rNumRepetitionsInParallelLoop)
{
	for (long i=0; i<rNumParallelLoops; i++)
	{
		for (int j=0; j<rNumRepetitionsInParallelLoop; j++)
		{
			if ((int)(pow(sqrt((double)(i+j)),3))%500000000==-1)
				std::cout << "should not happen" << std::endl;
		}
	}
}

void ComputeStatic(int rNumParallelLoops, int rNumRepetitionsInParallelLoop)
{
#pragma omp parallel for schedule(static) default(shared)
	for (long i=0; i<rNumParallelLoops; i++)
	{
		for (int j=0; j<rNumRepetitionsInParallelLoop; j++)
		{
			if ((int)(pow(sqrt((double)(i+j)),3))%500000000==-1)
				std::cout << "should not happen" << std::endl;
		}
	}
}

void ComputeGuided(int rNumParallelLoops, int rNumRepetitionsInParallelLoop)
{
#pragma omp parallel for schedule(guided) default(shared)
	for (long i=0; i<rNumParallelLoops; i++)
	{
		for (int j=0; j<rNumRepetitionsInParallelLoop; j++)
		{
			if ((int)(pow(sqrt((double)(i+j)),3))%500000000==-1)
				std::cout << "should not happen" << std::endl;
		}
	}
}

void ComputeDynamic(int rNumParallelLoops, int rNumRepetitionsInParallelLoop)
{
#pragma omp parallel for schedule(dynamic) default(shared)
	for (long i=0; i<rNumParallelLoops; i++)
	{
		for (int j=0; j<rNumRepetitionsInParallelLoop; j++)
		{
			if ((int)(pow(sqrt((double)(i+j)),3))%500000000==-1)
				std::cout << "should not happen" << std::endl;
		}
	}
}

int main()
{
#ifdef _OPENMP
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

    Eigen::MatrixXd speedUp(omp_get_num_procs()*2-1,13);

	std::cout << "serial execution " << std::endl;
	int numInitParallelRegions(1);
	int numRepetitionsInParallelLoop(1000);
	int numParallelLoops(5000);

	double wstart = omp_get_wtime ( );
	for (int j=0; j<numInitParallelRegions; j++)
	{
		ComputeNoParallel(numParallelLoops,numRepetitionsInParallelLoop);
	}
	double wend = omp_get_wtime ( );
	double timeWithOneThread = wend-wstart;

    for (int effortInLoop=0; effortInLoop<2; effortInLoop++)
    {
		switch (effortInLoop)
		{
		case 0:
			numRepetitionsInParallelLoop=1000;
			numParallelLoops = 50000;
			break;
		case 1:
			numRepetitionsInParallelLoop=1;
			numParallelLoops = 50000000;
			break;
		default:
			exit(-1);
		}
    	std::cout << "static load balancing " << std::endl;
		for (int numThreads=1; numThreads<omp_get_num_procs()*2; numThreads++)
		{
				omp_set_num_threads(numThreads);
				std::clock_t start,end;
				start=clock();
				double wstart = omp_get_wtime ( );

				for (int j=0; j<numInitParallelRegions; j++)
				{
					ComputeStatic(numParallelLoops,numRepetitionsInParallelLoop);
				}
				end=clock();
				double wend = omp_get_wtime ( );
				if (numThreads==1)
					timeWithOneThread = wend-wstart;
				std::cout <<"[OpenMPTest with " << numThreads << "] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") speedup=" << (timeWithOneThread)/(wend-wstart) << "\n";
				speedUp(numThreads-1,0) = numThreads;
				speedUp(numThreads-1,6*effortInLoop+1) = wend-wstart;
				speedUp(numThreads-1,6*effortInLoop+2) = (timeWithOneThread)/(wend-wstart);
		}

		std::cout << "dynamic load balancing " << std::endl;
		for (int numThreads=1; numThreads<omp_get_num_procs()*2; numThreads++)
		{
				omp_set_num_threads(numThreads);
				std::clock_t start,end;
				start=clock();
				double wstart = omp_get_wtime ( );

				for (int j=0; j<numInitParallelRegions; j++)
				{
					ComputeDynamic(numParallelLoops,numRepetitionsInParallelLoop);
				}
				end=clock();
				double wend = omp_get_wtime ( );
				std::cout <<"[OpenMPTest with " << numThreads << "] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") speedup=" << (timeWithOneThread)/(wend-wstart) << "\n";
				speedUp(numThreads-1,6*effortInLoop+3) = wend-wstart;
				speedUp(numThreads-1,6*effortInLoop+4) = (timeWithOneThread)/(wend-wstart);
		}

		std::cout << "guided load balancing " << std::endl;
		for (int numThreads=1; numThreads<omp_get_num_procs()*2; numThreads++)
		{
				omp_set_num_threads(numThreads);
				std::clock_t start,end;
				start=clock();
				double wstart = omp_get_wtime ( );

				for (int j=0; j<numInitParallelRegions; j++)
				{
					ComputeGuided(numParallelLoops,numRepetitionsInParallelLoop);
				}
				end=clock();
				double wend = omp_get_wtime ( );
				std::cout <<"[OpenMPTest with " << numThreads << "] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") speedup=" << (timeWithOneThread)/(wend-wstart) << "\n";
				speedUp(numThreads-1,6*effortInLoop+5) = wend-wstart;
				speedUp(numThreads-1,6*effortInLoop+6) = (timeWithOneThread)/(wend-wstart);
		}
    }

    NuTo::EigenCompanion::WriteToFile(speedUp, "SpeedUp.txt"," ");
#endif
    return 0;
}
