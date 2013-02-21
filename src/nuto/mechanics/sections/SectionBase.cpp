// $Id$
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION
#include <boost/noncopyable.hpp>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/sections/SectionBase.h"


//! @brief ... constructor
NuTo::SectionBase::SectionBase()
{
	mInputConstitutiveIsTemperature = false;
	mInputConstitutiveIsTemperatureGradient = false;
	mInputConstitutiveIsDamage = false;
	mInputConstitutiveIsDeformationGradient = true;
    mDisplacementDof = true;
    mTemperatureDof = false;
    mDamageDof = false;
}

int NuTo::SectionBase::GetNumInputConstitutive()const
{
	return mInputConstitutiveIsTemperature+mInputConstitutiveIsTemperatureGradient+mInputConstitutiveIsDeformationGradient;
}

int NuTo::SectionBase::GetNumDofs()const
{
	return mDisplacementDof+mTemperatureDof;
}

//! @brief... get if temperatures are dof
//! @return ... true, if temperatures are dofs
bool NuTo::SectionBase::GetIsTemperatureDof()const
{
	return mTemperatureDof;
}

//! @brief... set if damage are dofs
//! @param rFlag ... true, if damage are dofs
void NuTo::SectionBase::SetIsDamageDof(bool rFlag)
{
	mDamageDof = rFlag;
}

//! @brief... get if damage are dof
//! @return ... true, if damage are dofs
bool NuTo::SectionBase::GetIsDamageDof()const
{
	return mDamageDof;
}

//! @brief... set if temperatures are dofs
//! @param rFlag ... true, if temperatures are dofs
void NuTo::SectionBase::SetIsTemperatureDof(bool rFlag)
{
	mTemperatureDof = rFlag;
}

//! @brief... get if temperatures are to be used as input to the constitutive model
//! @return ... true, if temperatures are to be used as input to the constitutive model
bool NuTo::SectionBase::GetInputConstitutiveIsTemperature()const
{
	return mInputConstitutiveIsTemperature;
}
//! @brief... set if temperatures are to be used as input to the constitutive model
//! @param rFlag ... true, if temperatures are to be used as input to the constitutive model
void NuTo::SectionBase::SetInputConstitutiveIsTemperature(bool rFlag)
{
	mInputConstitutiveIsTemperature = rFlag;
}

//! @brief... get if damage is to be used as input to the constitutive model
//! @return ... true, if damage is to be used as input to the constitutive model
bool NuTo::SectionBase::GetInputConstitutiveIsDamage()const
{
	return mInputConstitutiveIsDamage;
}
//! @brief... set if damage is to be used as input to the constitutive model
//! @param rFlag ... true, if damage is to be used as input to the constitutive model
void NuTo::SectionBase::SetInputConstitutiveIsDamage(bool rFlag)
{
	mInputConstitutiveIsDamage = rFlag;
}

//! @brief... get if temperature gradients are to be used as input to the constitutive model
//! @return ... true, if temperature gradients are to be used as input to the constitutive model
bool NuTo::SectionBase::GetInputConstitutiveIsTemperatureGradient()const
{
	return mInputConstitutiveIsTemperatureGradient;
}
//! @brief... set if temperature gradients are to be used as input to the constitutive model
//! @param rFlag ... true, if temperature gradients are to be used as input to the constitutive model
void NuTo::SectionBase::SetInputConstitutiveIsTemperatureGradient(bool rFlag)
{
	mInputConstitutiveIsTemperatureGradient = rFlag;
}

//! @brief... get if displacements are dof
//! @return ... true, if displacements are dofs
bool NuTo::SectionBase::GetIsDisplacementDof()const
{
	return mDisplacementDof;
}

//! @brief... set if displacements are dofs
//! @param rFlag ... true, if displacements are dofs
void NuTo::SectionBase::SetIsDisplacementDof(bool rFlag)
{
	mDisplacementDof = rFlag;
}

//! @brief... get if rotations are dof
//! @return ... true, if rotations are dofs
bool NuTo::SectionBase::GetIsRotationDof()const
{
	return mDisplacementDof;
}

//! @brief... set if rotations are dofs
//! @param rFlag ... true, if rotations are dofs
void NuTo::SectionBase::SetIsRotationDof(bool rFlag)
{
	mDisplacementDof = rFlag;
}

//! @brief... get if deformation Gradient is to be used as input to the constitutive model
//! @return ... true, if deformation Gradient are to be used as input to the constitutive model
bool NuTo::SectionBase::GetInputConstitutiveIsDeformationGradient()const
{
	return mInputConstitutiveIsDeformationGradient;
}

//! @brief... set if deformation Gradient are to be used as input to the constitutive model
//! @param rFlag ... true, if deformation Gradient are to be used as input to the constitutive model
void NuTo::SectionBase::SetInputConstitutiveIsDeformationGradient(bool rFlag)
{
	mInputConstitutiveIsDeformationGradient = rFlag;
}

double NuTo::SectionBase::GetArea() const
{
    throw NuTo::MechanicsException("[NuTo::SectionBase::GetArea] section type has no cross-section area.") ;
}

void NuTo::SectionBase::SetArea(double rArea)
{
    throw NuTo::MechanicsException("[NuTo::SectionBase::SetArea] section type has no cross-section area.") ;
}

double NuTo::SectionBase::GetThickness() const
{
    throw NuTo::MechanicsException("[NuTo::SectionBase::GetThickness] section type has no thickness.");
}

void NuTo::SectionBase::SetThickness(double rThickness)
{
    throw NuTo::MechanicsException("[NuTo::SectionBase::SetThickness] section type has no thickness.");
}

void NuTo::SectionBase::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    section pointer: " << this << std::endl;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::SectionBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::SectionBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::SectionBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::SectionBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::SectionBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::SectionBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::SectionBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SectionBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_NVP(mInputConstitutiveIsDamage)
       & BOOST_SERIALIZATION_NVP(mInputConstitutiveIsTemperature)
       & BOOST_SERIALIZATION_NVP(mInputConstitutiveIsTemperatureGradient)
       & BOOST_SERIALIZATION_NVP(mInputConstitutiveIsDeformationGradient)
       & BOOST_SERIALIZATION_NVP(mDamageDof)
       & BOOST_SERIALIZATION_NVP(mDisplacementDof)
       & BOOST_SERIALIZATION_NVP(mTemperatureDof);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SectionBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SectionBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::SectionBase)
#endif // ENABLE_SERIALIZATION
