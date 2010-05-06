// $ld: $ 
// IpDataBase.cpp
// created Apr 29, 2010 by Joerg F. Unger
#include <iostream>
#include <string>
#include <sstream>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/IpDataBase.h"

NuTo::IpDataBase::~IpDataBase()
{
	std::cout << "call Desctructor [NuTo::IpDataBase::~IpDataBase]." << std::endl;
}

//! @brief adds the weight to an integration point, eventually reallocates the data
//! @param rNonlocalElement the Element (local number from the nonlocal elements)
//! @param rNonlocalIp integration point of the nonlocal element
//! @param rNumIps number of integration points of the nonlocal element (for allocation purpose of not existing)
//! @param rWeight nonlocal weight
void NuTo::IpDataBase::SetNonlocalWeight(int rElement,int rNonlocalIp,int rNumIps, double rWeight)
{
	throw NuTo::MechanicsException("[IpDataBase::SetNonlocalWeight] This Ip data type cannot store nonlocal weights - check the ip data type.");
}

//! @brief return the nonlocal weights
//! @param rNonlocalElement nonlocal element (between 0 and nonlocal elements.size stored in nonlocal element data)
//! @return nonlocal weights
const std::vector<double>& NuTo::IpDataBase::GetNonlocalWeights(int rNonlocalElement)const
{
	throw NuTo::MechanicsException("[IpDataBase::GetNonlocalWeights] This Ip data type cannot store nonlocal weights - check the ip data type.");
}
