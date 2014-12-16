// $Id: ElementDataNonlocalBase.cpp 551 2011-06-24 18:24:23Z unger3 $
// ElementDataNonlocalBase.cpp
// created Apr 22, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"


NuTo::ElementOutputBase::ElementOutputBase()
{}


NuTo::ElementOutputBase::~ElementOutputBase()
{}

NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& NuTo::ElementOutputBase::GetFullMatrixDouble()
{
	throw MechanicsException("[ElementOutputBase::GetFullMatrixDouble] element output matrix is not of type FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>");
}

NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& NuTo::ElementOutputBase::GetFullMatrixInt()
{
	throw MechanicsException("[ElementOutputBase::GetFullMatrixDouble] element output matrix is not of type FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>");
}

NuTo::FullVector<double,Eigen::Dynamic>& NuTo::ElementOutputBase::GetFullVectorDouble()
{
	throw MechanicsException("[ElementOutputBase::GetFullVectorDouble] element output matrix is not of type FullVector<double,Eigen::Dynamic>");
}

NuTo::FullVector<int,Eigen::Dynamic>& NuTo::ElementOutputBase::GetFullVectorInt()
{
	throw MechanicsException("[ElementOutputBase::GetFullVectorInt] element output matrix is not of type FullVector<double,Eigen::Dynamic>");
}

std::vector<int>& NuTo::ElementOutputBase::GetVectorInt()
{
	throw MechanicsException("[ElementOutputBase::GetFullMatrixDouble] element output matrix is not of type std::vector<int>");
}

NuTo::IpData::eIpStaticDataType NuTo::ElementOutputBase::GetIpDataType()
{
	throw MechanicsException("[ElementOutputBase::GetIpDataType] ipdata is not stored.");
}

NuTo::FatigueData::eFatigueDataType NuTo::ElementOutputBase::GetFatigueDataType()
{
	throw MechanicsException("[ElementOutputBase::GetFatigueDataType] fatigue data is not stored.");
}

void NuTo::ElementOutputBase::SetSymmetry(bool rSymmetric)
{
	throw MechanicsException("[ElementOutputBase::SetSymmetry] symmetry is not stored.");
}

bool NuTo::ElementOutputBase::GetSymmetry()const
{
	throw MechanicsException("[ElementOutputBase::SetSymmetry] symmetry is not stored.");
}

void NuTo::ElementOutputBase::SetConstant(bool rConstant)
{
	throw MechanicsException("[ElementOutputBase::SetConstant] constness is not stored.");
}

bool NuTo::ElementOutputBase::GetConstant()const
{
	throw MechanicsException("[ElementOutputBase::GetConstant] constness is not stored.");
}

NuTo::ElementOutputBase* NuTo::ElementOutputBase::Clone() const
{
	return 0;
}


NuTo::ElementOutputBase* NuTo::new_clone( const NuTo::ElementOutputBase& o)
{
	return o.Clone();
}

