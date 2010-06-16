#include <iostream>
#include <ctime>
#include <cstdio>

#include "VirtualFunctionsVsDynamicCastClasses.h"

int main()
{
    unsigned int maxLoop(100000000);

	AB* OrigPtr = new AB;
	Base* basePtr = OrigPtr;
    double r[3];

    std::clock_t start,end;

    // direct access via public variables
    start=clock();
    for (unsigned int count=0; count<maxLoop ; count++)
    {
    	r[0] = OrigPtr->mB[0];
    	r[1] = OrigPtr->mB[1];
       	r[2] = OrigPtr->mB[2];
    }
    end=clock();
    std::cout<<"time required for "<< maxLoop << " direct access as public variable " << difftime(end,start)/CLOCKS_PER_SEC << " ("<< 1.0 << ")" << std::endl;
    //store as reference
    double ref(difftime(end,start));

    // access via base class and virtual function call
    start=clock();
    for (unsigned int count=0; count<maxLoop ; count++)
    	basePtr->FunctionB(r);
    end=clock();
    std::cout<<"time required for "<< maxLoop << " virtual function calls           " << difftime(end,start)/CLOCKS_PER_SEC << " ("<< difftime(end,start)/ref << ")" << std::endl;

    // access via original class and virtual function call
    start=clock();
    for (unsigned int count=0; count<maxLoop ; count++)
    	OrigPtr->FunctionB(r);
    end=clock();
    std::cout<<"time required for "<< maxLoop << " direct function calls            " << difftime(end,start)/CLOCKS_PER_SEC << " ("<< difftime(end,start)/ref << ")" << std::endl;

    // dynamic cast to derived ptr
    start=clock();
    for (unsigned int count=0; count<maxLoop ; count++)
    	dynamic_cast<AB*>(basePtr)->FunctionB(r);
    end=clock();
    std::cout<<"time required for "<< maxLoop << " dynamic casts                    " << difftime(end,start)/CLOCKS_PER_SEC << " ("<< difftime(end,start)/ref << ")" << std::endl;

    // explicit cast to derived ptr
    start=clock();
    for (unsigned int count=0; count<maxLoop ; count++)
    	reinterpret_cast<AB*>(basePtr)->FunctionB(r);
    end=clock();
    std::cout<<"time required for "<< maxLoop << " reinterpret_casts                " << difftime(end,start)/CLOCKS_PER_SEC << " ("<< difftime(end,start)/ref << ")" << std::endl;

    // call via asAB ptr
    start=clock();
    for (unsigned int count=0; count<maxLoop ; count++)
    	basePtr->asAB()->FunctionB(r);
    end=clock();
    std::cout<<"time required for "<< maxLoop << " via as DerivedPtr                " << difftime(end,start)/CLOCKS_PER_SEC << " ("<< difftime(end,start)/ref << ")" << std::endl;

    // just for comparison a simple multipication
    start=clock();
    double t(0.999);
    for (unsigned int count=0; count<maxLoop ; count++)
    	t = 0.999*t;
    end=clock();
    std::cout<<"time required for "<< maxLoop << " multiplication (double)          " << difftime(end,start)/CLOCKS_PER_SEC << " ("<< difftime(end,start)/ref << ")" << std::endl;
    if (t>1)
         printf("t %g",t);

    // just for comparison a simple add
    start=clock();
    t = 2;
    for (unsigned int count=0; count<maxLoop ; count++)
    	t = t+t;
    if (t>100)
    	t=0;
    end=clock();
    std::cout<<"time required for "<< maxLoop << " add (double)                     " << difftime(end,start)/CLOCKS_PER_SEC << " ("<< difftime(end,start)/ref << ")" << std::endl;
    if (t>100)
          printf("d %g",t);

    // just for comparison a simple multipication
    start=clock();
    long i(2);
    for (unsigned int count=0; count<maxLoop ; count++)
    {
    	i = 2*i;
    	if (i>100)
    		i=2;
    }
    end=clock();
    std::cout<<"time required for "<< maxLoop << " multiplication (long)            " << difftime(end,start)/CLOCKS_PER_SEC << " ("<< difftime(end,start)/ref << ")" << std::endl;
    if (i>100)
          printf("d %d",(int)i);

    // just for comparison a simple add
    start=clock();
    i = 2;
    for (unsigned int count=0; count<maxLoop ; count++)
    	i = i+i;
    if (i>100)
    	i=0;
    end=clock();
    std::cout<<"time required for "<< maxLoop << " add (long)                       " << difftime(end,start)/CLOCKS_PER_SEC << " ("<< difftime(end,start)/ref << ")" << std::endl;
    if (t>100)
          printf("d %d",(int)i);
}

