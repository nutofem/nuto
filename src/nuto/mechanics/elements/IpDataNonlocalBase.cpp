// $ld: $ 
// IpDataNonlocalBase.cpp
// created Apr 28, 2010 by Joerg F. Unger

#include <assert.h>
#include <iostream>
#include "nuto/mechanics/elements/IpDataNonlocalBase.h"

NuTo::IpDataNonlocalBase::~IpDataNonlocalBase()
{
	//std::cout << std::endl << "call Desctructor [NuTo::IpDataNonlocalBase::~IpDataNonlocalBase]." << std::endl;
}


//! @brief adds the weight to an integration point, eventually reallocates the data
//! @param rNonlocalElement the Element (local number from the nonlocal elements)
//! @param rNonlocalIp integration point of the nonlocal element
//! @param rNumIps number of integration points of the nonlocal element (for allocation purpose of not existing)
//! @param rWeight nonlocal weight
void NuTo::IpDataNonlocalBase::SetNonlocalWeight(int rNonlocalElement,int rNonlocalIp,int rNumIps, double rWeight)
{
    if (rNonlocalElement<(int)mNonlocalWeights.size())
    {
    	assert((int)mNonlocalWeights[rNonlocalElement].size()==rNumIps);
    	assert(rNonlocalIp<rNumIps);
    	mNonlocalWeights[rNonlocalElement][rNonlocalIp]=rWeight;
    }
    else
    {
    	// the nonlocal element has to be added
    	assert(rNonlocalElement==(int)mNonlocalWeights.size());
    	mNonlocalWeights.push_back(std::vector<double>(rNumIps,0));
    	std::cout << "[NuTo::IpDataNonlocalBase::SetNonlocalWeight] " << rNonlocalIp << " " << rNumIps << std::endl;
    	assert(rNonlocalIp<rNumIps);
    	mNonlocalWeights[rNonlocalElement][rNonlocalIp]=rWeight;
    }
}
//! @brief return the nonlocal weights
//! @param rNonlocalElement nonlocal element (between 0 and nonlocal elements.size stored in nonlocal element data)
//! @return nonlocal weights
const std::vector<double>& NuTo::IpDataNonlocalBase::GetNonlocalWeights(int rNonlocalElement)const
{
    assert(rNonlocalElement>=0 && rNonlocalElement<(int)mNonlocalWeights.size());
	return mNonlocalWeights[rNonlocalElement];
}
