// $Id$
#ifdef SHOW_TIME
    #include <ctime>
# ifdef _OPENMP
    #include <omp.h>
# endif
#endif


#include <iostream>
#include "math/FullMatrix.h"
//#include <boost/array.hpp>

int main()
{
    //vector size
   int size = 900000;
   int iter = 1000;
   //int size = 5;

	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Amat(size,1);
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Bmat(size,1);

	Eigen::VectorXd A(size);
	Eigen::VectorXd B(size);
	Eigen::VectorXd C(size);

	std::vector<double> g(size);
    std::vector<double> h(size);
    std::vector<double> k(size);

    double* p;
    double* q;
    double* r;

    p=new double[size];
    q=new double[size];
    r=new double[size];

    //boost::array<double,size> x; //size has to be a constant value, no variable

 	//initialize
    /* initialize random seed: */
    srand ( time(NULL) );
    for (int i=0;i<size;++i)
    {
    	Amat(i,0)=rand()/(static_cast<double>(RAND_MAX)+1.);
    	Bmat(i,0)=rand()/(static_cast<double>(RAND_MAX)+1.);

       	h.at(i)=Amat(i,0);
       	k.at(i)=Bmat(i,0);

       	q[i]=h.at(i);
       	r[i]=k.at(i);

    }
    A=Amat;
    B=Bmat;

    double result=0;
#ifdef SHOW_TIME
    std::clock_t start,end;
    double diff=0;
    start=clock();
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
#endif //SHOW_TIME
	for (int i=0;i<iter;++i)
	{
		//scalar product
		result += A.dot(B);
	}

#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    diff=wend-wstart;
    	std::cout<<"[NuTo::VectorOperations] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") - "<<(wend-wstart)/diff<<"\n";
#else
        std::cout<<"[NuTo::VectorOperations] " << (difftime(end,start)/CLOCKS_PER_SEC) << "sec  " << std::endl;
#endif
#endif
       std::cout<<"[NuTo::VectorOperations] "<< "result of scalar product with eigen vectors "<< result <<std::endl;

    result=0;
#ifdef SHOW_TIME
    start=clock();
#ifdef _OPENMP
    wstart = omp_get_wtime ( );
#endif
#endif //SHOW_TIME
	for (int num=0;num<iter;++num)
	{
		for (int i=0;i<size;++i)
		{
			//scalar product
			//result += h.at(i)*k.at(i);
			result += h[i]*k[i];
		}
	}

#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    wend = omp_get_wtime ( );
    	std::cout<<"[NuTo::VectorOperations] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") - "<<(wend-wstart)/diff<<"\n";
#else
        std::cout<<"[NuTo::VectorOperations] " << (difftime(end,start)/CLOCKS_PER_SEC)<< "sec  " << std::endl;
#endif
#endif
        std::cout<<"[NuTo::VectorOperations] "<< "result of scalar product with std::vector "<< result <<std::endl;

    result=0;
#ifdef SHOW_TIME
    start=clock();
#ifdef _OPENMP
    wstart = omp_get_wtime ( );
#endif
#endif //SHOW_TIME
	for (int num=0;num<iter;++num)
	{
		for (int i=0;i<size;++i)
		{
			//scalar product
			result += r[i]*q[i];
		}
	}
/*
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    wend = omp_get_wtime ( );
    	std::cout<<"[NuTo::VectorOperations] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") - "<<(wend-wstart)/diff<<"\n";
#else
       std::cout<<"[NuTo::VectorOperations] " << (difftime(end,start)/CLOCKS_PER_SEC)<< "sec  " << std::endl;
#endif
#endif
       std::cout<<"[NuTo::VectorOperations] "<< "result of scalar product with double field  "<< result <<std::endl;
*/
/*
#ifdef SHOW_TIME
    start=clock();
#ifdef _OPENMP
    wstart = omp_get_wtime ( );
#endif
#endif //SHOW_TIME
	Eigen::VectorXd D(size);

//	for (int i=0;i<iter;++i)
	//{
		Amat.mEigenMatrix=A;
		Bmat.mEigenMatrix=B;

		C = Amat.mEigenMatrix;
		D = Bmat.mEigenMatrix;
	//}

#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    wend = omp_get_wtime ( );
    	std::cout<<"[NuTo::VectorOperations:  Transforme matrix to eigen and back] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") - "<<(wend-wstart)/diff<<"\n";
#else
   std::cout<<"[NuTo::VectorOperations: Transforme matrix to eigen and back] " << difftime(end,start)/CLOCKS_PER_SEC << "sec  " << std::endl;
#endif
#endif
*/
#ifdef SHOW_TIME
    start=clock();
#ifdef _OPENMP
    wstart = omp_get_wtime ( );
#endif
#endif //SHOW_TIME
    //scalar factor
    double alpha =.3;
	for (int i=0;i<iter;++i)
	{
		C = A+alpha*B;
	}
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    wend = omp_get_wtime ( );
    diff=wend-wstart;
    	std::cout<<"[NuTo::VectorOperations: linear combination scalar factor] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") - "<<(wend-wstart)/diff<<"\n";
#else
    std::cout<<"[NuTo::VectorOperations: linear combination scalar factor] " << difftime(end,start)/CLOCKS_PER_SEC << "sec  " << std::endl;
#endif
#endif

#ifdef SHOW_TIME
    start=clock();
#ifdef _OPENMP
    wstart = omp_get_wtime ( );
#endif
#endif //SHOW_TIME
    //scalar factor
	for (int num=0;num<iter;++num)
	{
		for (int i=0;i<size;++i)
		{
			//g.at(i) = k.at(i)+alpha*h.at(i);
			g[i] = k[i]+alpha*h[i];

		}

	}
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    wend = omp_get_wtime ( );
    	std::cout<<"[NuTo::VectorOperations: linear combination scalar factor with std::vector [i]] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") - "<<(wend-wstart)/diff<<"\n";
#else
    std::cout<<"[NuTo::VectorOperations: linear combination scalar factor with std::vector] " << difftime(end,start)/CLOCKS_PER_SEC << "sec  " << std::endl;
#endif
#endif


#ifdef SHOW_TIME
    start=clock();
#ifdef _OPENMP
    wstart = omp_get_wtime ( );
#endif
#endif //SHOW_TIME
    //scalar factor
	for (int num=0;num<iter;++num)
	{
		for (int i=0;i<size;++i)
		{
			k[i] += alpha*h[i];

		}

	}
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    wend = omp_get_wtime ( );
    	std::cout<<"[NuTo::VectorOperations: k+=alpha h with std::vector [i]] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") - "<<(wend-wstart)/diff<<"\n";
#else
    std::cout<<"[NuTo::VectorOperations: linear combination scalar factor with std::vector] " << difftime(end,start)/CLOCKS_PER_SEC << "sec  " << std::endl;
#endif
#endif

#ifdef SHOW_TIME
    start=clock();
#ifdef _OPENMP
    wstart = omp_get_wtime ( );
#endif
#endif //SHOW_TIME
    //scalar factor
	for (int num=0;num<iter;++num)
	{
		for (int i=0;i<size;++i)
		{
			k[i] *= alpha;
			k[i]+=h[i];

		}

	}
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    wend = omp_get_wtime ( );
    	std::cout<<"[NuTo::VectorOperations: k*=alpa,k+=h with std::vector [i]] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") - "<<(wend-wstart)/diff<<"\n";
#else
    std::cout<<"[NuTo::VectorOperations: linear combination scalar factor with std::vector] " << difftime(end,start)/CLOCKS_PER_SEC << "sec  " << std::endl;
#endif
#endif

/*
#ifdef SHOW_TIME
    start=clock();
#ifdef _OPENMP
    wstart = omp_get_wtime ( );
#endif
#endif //SHOW_TIME
    //scalar factor
	for (int num=0;num<iter;++num)
	{
		for (int i=0;i<size;++i)
		{
			p[i] = q[i]+alpha*r[i];
		}

	}
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    wend = omp_get_wtime ( );
    	std::cout<<"[NuTo::VectorOperations: linear combination scalar factor with double] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") - "<<(wend-wstart)/diff<<"\n";
#else
    std::cout<<"[NuTo::VectorOperations: linear combination scalar factor with double] " << difftime(end,start)/CLOCKS_PER_SEC << "sec  " << std::endl;
#endif
#endif


#ifdef SHOW_TIME
    start=clock();
#ifdef _OPENMP
    wstart = omp_get_wtime ( );
#endif
#endif //SHOW_TIME
    //scalar factor
	for (int num=0;num<iter;++num)
	{
		for (int i=0;i<size;++i)
		{
			s[i] = t[i]+alpha*u[i];
		}

	}
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    wend = omp_get_wtime ( );
    	std::cout<<"[NuTo::VectorOperations: linear combination scalar factor with Array] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<") - "<<(wend-wstart)/diff<<"\n";
#else
    std::cout<<"[NuTo::VectorOperations: linear combination scalar factor with Array] " << difftime(end,start)/CLOCKS_PER_SEC << "sec  " << std::endl;
#endif
#endif
*/
    delete[] p;
    delete[] q;
    delete[] r;

}
