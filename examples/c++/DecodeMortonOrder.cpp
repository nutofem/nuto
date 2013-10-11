// $Id$

#include <iostream>
#include "nuto/base/NuToObject.h"
#include "nuto/math/MortonOrder.h"

int main()
{
    // test decode morton key
	uint32_t x=0,
			y=0,
			z=0,
			key=0;
	std::cout<<"Decode 3D coordinates (integer) from Morton key.\nEnter key:\n";
	std::cin>>key;;
	x=NuTo::MortonOrder::DecodeMorton3X(key);
	y=NuTo::MortonOrder::DecodeMorton3Y(key);
	z=NuTo::MortonOrder::DecodeMorton3Z(key);
	std::cout<<" result: x="<<x<<", y="<<y<<", z="<<z<<" \n";
	return 0;
}
