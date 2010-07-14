// $ld: $ 
// VisualizeComponentNonlocalWeight.cpp
// created Apr 27, 2010 by Joerg F. Unger

#include "nuto/visualize/VisualizeComponentNonlocalWeight.h"
#include "nuto/visualize/VisualizeException.h"
#include <sstream>

NuTo::VisualizeComponentNonlocalWeight::VisualizeComponentNonlocalWeight(const ElementBase* rElement, int rElementId, int rIp) : VisualizeComponentBase::VisualizeComponentBase()
{
    mElement = rElement;
    mElementId = rElementId;
    mIp = rIp;
}

int NuTo::VisualizeComponentNonlocalWeight::GetElementId()const
{
	return mElementId;
}

const NuTo::ElementBase* NuTo::VisualizeComponentNonlocalWeight::GetElement()const
{
	return mElement;
}

int NuTo::VisualizeComponentNonlocalWeight::GetIp()const
{
	return mIp;
}

std::string NuTo::VisualizeComponentNonlocalWeight::GetComponentName()const
{
	std::stringstream out;
	out << "NonlocalWeight_Element" << mElementId << "_Ip_" << mIp;
	return out.str();
}
