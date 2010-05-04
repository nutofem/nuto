// $ld: $ 
// VisualizeComponentNonlocalWeight.cpp
// created Apr 27, 2010 by Joerg F. Unger

#include "nuto/visualize/VisualizeComponentNonlocalWeight.h"
#include "nuto/visualize/VisualizeException.h"
#include <sstream>

NuTo::VisualizeComponentNonlocalWeight::VisualizeComponentNonlocalWeight() : VisualizeComponentBase::VisualizeComponentBase()
{}

void NuTo::VisualizeComponentNonlocalWeight::SetElementIp(int rElementId, int rIp)
{
    mElementId = rElementId;
    mIp = rIp;
}


int NuTo::VisualizeComponentNonlocalWeight::GetElementId()const
{
	return mElementId;
}

int NuTo::VisualizeComponentNonlocalWeight::GetIp()const
{
	return mIp;
}

std::string NuTo::VisualizeComponentNonlocalWeight::GetComponentName()const
{
	std::stringstream out;
	out << "NONLOCAL_WEIGHT_ELEM_" << mElementId << "_IP_" << mIp;
	return out.str();
}
