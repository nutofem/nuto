// $Id$
#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux3D.h"

// constructor
NuTo::ConstitutiveOutputBase::ConstitutiveOutputBase()
{
}

NuTo::EngineeringStrain3D& NuTo::ConstitutiveOutputBase::GetEngineeringStrain3D()
{
	throw MechanicsException("[NuTo::ConstitutiveOutputBase::GetEngineeringStrain3D] not implemented for this output object.");
}

NuTo::EngineeringStress1D& NuTo::ConstitutiveOutputBase::GetEngineeringStress1D()
{
	throw MechanicsException("[NuTo::ConstitutiveOutputBase::GetEngineeringStress1D] not implemented for this output object.");
}

NuTo::EngineeringStress2D& NuTo::ConstitutiveOutputBase::GetEngineeringStress2D()
{
	throw MechanicsException("[NuTo::ConstitutiveOutputBase::GetEngineeringStress2D] not implemented for this output object.");
}

NuTo::EngineeringStress3D& NuTo::ConstitutiveOutputBase::GetEngineeringStress3D()
{
	throw MechanicsException("[NuTo::ConstitutiveOutputBase::GetEngineeringStress3D] not implemented for this output object.");
}

NuTo::ConstitutiveTangentLocal6x6& NuTo::ConstitutiveOutputBase::GetConstitutiveTangentLocal6x6()
{
	throw MechanicsException("[NuTo::ConstitutiveOutputBase::GetConstitutiveTangentLocal6x6] not implemented for this output object.");
}

NuTo::ConstitutiveTangentLocal3x3& NuTo::ConstitutiveOutputBase::GetConstitutiveTangentLocal3x3()
{
	throw MechanicsException("[NuTo::ConstitutiveOutputBase::GetConstitutiveTangentLocal3x3] not implemented for this output object.");
}

NuTo::Damage& NuTo::ConstitutiveOutputBase::GetDamage()
{
	throw MechanicsException("[NuTo::ConstitutiveOutputBase::GetDamage] not implemented for this output object.");
}

NuTo::HeatFlux3D& NuTo::ConstitutiveOutputBase::GetHeatFlux3D()
{
	throw MechanicsException("[NuTo::ConstitutiveOutputBase::GetHeatFlux3D] not implemented for this output object.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
NuTo::ConstitutiveTangentLocal<1,1>& NuTo::ConstitutiveOutputBase::AsConstitutiveTangentLocal_1x1()
{
    throw MechanicsException("[ConstitutiveTangentBase::ConstitutiveTangentBase] matrix is not of type tangent local 1x1.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
NuTo::ConstitutiveTangentLocal<2,2>& NuTo::ConstitutiveOutputBase::AsConstitutiveTangentLocal_2x2()
{
    throw MechanicsException("[ConstitutiveTangentBase::ConstitutiveTangentBase] matrix is not of type tangent local 2x2.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
NuTo::ConstitutiveTangentLocal<3,3>& NuTo::ConstitutiveOutputBase::AsConstitutiveTangentLocal_3x3()
{
    throw MechanicsException("[ConstitutiveTangentBase::ConstitutiveTangentBase] matrix is not of type tangent local 3x3.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
NuTo::ConstitutiveTangentLocal<6,6>& NuTo::ConstitutiveOutputBase::AsConstitutiveTangentLocal_6x6()
{
    throw MechanicsException("[ConstitutiveTangentBase::ConstitutiveTangentBase] matrix is not of type tangent local 6x6.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
NuTo::ConstitutiveTangentLocal<6,1>& NuTo::ConstitutiveOutputBase::AsConstitutiveTangentLocal_6x1()
{
    throw MechanicsException("[ConstitutiveTangentBase::ConstitutiveTangentBase] matrix is not of type tangent local 6x1.");
}

//! @brief return part of the nonlocal matrix
NuTo::ConstitutiveTangentLocal<1,1>& NuTo::ConstitutiveOutputBase::GetSubMatrix_1x1(int rSubMatrix)
{
    throw MechanicsException("[ConstitutiveTangentBase::GetSubMatrix_1x1] matrix is not of type tangent nonlocal  1x1.");
}

//! @brief return part of the nonlocal matrix
NuTo::ConstitutiveTangentLocal<2,2>& NuTo::ConstitutiveOutputBase::GetSubMatrix_2x2(int rSubMatrix)
{
    throw MechanicsException("[ConstitutiveTangentBase::GetSubMatrix_2x2] matrix is not of type tangent nonlocal  2x2.");
}

//! @brief return part of the nonlocal matrix
NuTo::ConstitutiveTangentLocal<3,3>& NuTo::ConstitutiveOutputBase::GetSubMatrix_3x3(int rSubMatrix)
{
    throw MechanicsException("[ConstitutiveTangentBase::GetSubMatrix_3x3] matrix is not of type tangent nonlocal  3x3.");
}

//! @brief return part of the nonlocal matrix
NuTo::ConstitutiveTangentLocal<6,6>& NuTo::ConstitutiveOutputBase::GetSubMatrix_6x1(int rSubMatrix)
{
    throw MechanicsException("[ConstitutiveTangentBase::GetSubMatrix_6x1] matrix is not of type tangent nonlocal  6x1.");
}

//! @brief return part of the nonlocal matrix
NuTo::ConstitutiveTangentLocal<6,1>& NuTo::ConstitutiveOutputBase::GetSubMatrix_6x6(int rSubMatrix)
{
    throw MechanicsException("[ConstitutiveTangentBase::GetSubMatrix_6x6] matrix is not of type tangent nonlocal  6x6.");
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstitutiveOutputBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveOutputBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveOutputBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveOutputBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveOutputBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveOutputBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveOutputBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveOutputBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveOutputBase" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
