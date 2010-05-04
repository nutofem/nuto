// $ld: $ 
// IpDataNonlocalBase.cpp
// created Apr 28, 2010 by Joerg F. Unger

#include "nuto/mechanics/elements/IpDataNonlocalBase.h"
#include <assert.h>

void NuTo::IpDataNonlocalBase::SetNonlocalWeights(int rNonlocalElement, const std::vector<double>& rWeights)
{
    if (rNonlocalElement<(int)mNonlocalWeights.size())
    	mNonlocalWeights[rNonlocalElement]=rWeights;
    else
    	mNonlocalWeights.push_back(rWeights);
}

const std::vector<double>& NuTo::IpDataNonlocalBase::GetNonlocalWeights(int rNonlocalElement)const
{
    assert(rNonlocalElement>=0 && rNonlocalElement<(int)mNonlocalWeights.size());
	return mNonlocalWeights[rNonlocalElement];
}
