// $Id: ConstitutiveLatticeStressStrain.cpp 112 2009-11-17 16:43:15Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveLatticeStressStrain.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStress3D.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/elements/ElementBase.h"

NuTo::ConstitutiveLatticeStressStrain::ConstitutiveLatticeStressStrain() : ConstitutiveBase()
{
	mEnergyFlag = true;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveLatticeStressStrain::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveLatticeStressStrain::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveLatticeStressStrain::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveLatticeStressStrain::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveLatticeStressStrain::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveLatticeStressStrain::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveLatticeStressStrain::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveLatticeStressStrain" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
       & BOOST_SERIALIZATION_NVP(mEnergyFlag);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveLatticeStressStrain" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveLatticeStressStrain)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstitutiveLatticeStressStrain)
#endif // ENABLE_SERIALIZATION


//  Lattice strain /////////////////////////////////////
//! @brief ... convert an Lattice strain from 2D into 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain2D ... Lattice strain 2D (input)
//! @param rLatticeStrain3D ... Lattice strain 3D (output)
NuTo::Error::eError NuTo::ConstitutiveLatticeStressStrain::GetLatticeStrain(const ElementBase* rElement, int rIp,
                                  const LatticeStrain2D& rLatticeStrain2D, LatticeStrain3D& rLatticeStrain3D) const
{
	rLatticeStrain3D.mLatticeStrain[0] = rLatticeStrain2D.mLatticeStrain[0];
	rLatticeStrain3D.mLatticeStrain[1] = rLatticeStrain2D.mLatticeStrain[1];
	rLatticeStrain3D.mLatticeStrain[2] = rLatticeStrain2D.mLatticeStrain[2];
	rLatticeStrain3D.mLatticeStrain[3] = 0;

	return Error::SUCCESSFUL;
}

//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeStressStrain::GetInternalEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
		const LatticeStrain2D& rLatticeStrain, double& rEnergy) const
{
    throw MechanicsException("[ConstitutiveLatticeStressStrain::GetInternalEnergy_LatticeStress_LatticeStrain] to be implemented.");
}

//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeStressStrain::GetInternalEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
		const LatticeStrain3D& rLatticeStrain, double& rEnergy) const
{
    throw MechanicsException("[ConstitutiveLatticeStressStrain::GetInternalEnergy_LatticeStress_LatticeStrain] to be implemented.");
}

//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeStressStrain::GetElasticEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
        const LatticeStrain2D& rLatticeStrain, double& rEnergy) const
{
    throw MechanicsException("[ConstitutiveLatticeStressStrain::GetElasticEnergy_LatticeStress_LatticeStrain] to be implemented.");
}

//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeStressStrain::GetElasticEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
        const LatticeStrain3D& rLatticeStrain, double& rEnergy) const
{
    throw MechanicsException("[ConstitutiveLatticeStressStrain::GetElasticEnergy_LatticeStress_LatticeStrain] to be implemented.");
}
//! @brief ... avoid dynamic cast
//! @return ... see brief explanation
NuTo::ConstitutiveLatticeStressStrain* NuTo::ConstitutiveLatticeStressStrain::AsConstitutiveLatticeStressStrain()
{
	return this;
}

//! @brief ... avoid dynamic cast
//! @return ... see brief explanation
const NuTo::ConstitutiveLatticeStressStrain* NuTo::ConstitutiveLatticeStressStrain::AsConstitutiveLatticeStressStrain()const
{
	return this;
}

//! @brief ... checks, if a model has to be switched from linear to nonlinear, and then performs the adaption
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rLatticeStrain ... deformation gradient
NuTo::Error::eError NuTo::ConstitutiveLatticeStressStrain::MultiscaleSwitchToNonlinear(ElementBase* rElement, int rIp,
        const LatticeStrain2D& rLatticeStrain)const
{
	throw MechanicsException("[NuTo::ConstitutiveLatticeStressStrain::MultiscaleSwitchToNonlinear] not implemented for this material routine.");
}
