// $Id$

#include <iostream>
#include "nuto/base/NuToObject.h"
#include "nuto/math/MortonOrder.h"

int main()
{
    // test morton
	uint32_t x=2;
	uint32_t y=1;
	uint32_t z=1;
	std::cout<<"Encode 3D coordinates (integer) as Morton key.\nEnter x y z value:\n";
	std::cin>>x>>y>>z;;
	uint32_t key=NuTo::MortonOrder::EncodeMorton3D(x,y,z);
	std::cout<<" result "<<key<<"\n";
	return 0;
}
