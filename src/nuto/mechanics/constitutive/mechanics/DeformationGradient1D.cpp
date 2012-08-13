// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"

#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"

#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"

#include "nuto/mechanics/MechanicsException.h"

// constructor
NuTo::DeformationGradient1D::DeformationGradient1D(): ConstitutiveInputBase::ConstitutiveInputBase()
{
    this->mDeformationGradient = 0.0;
}

//! @brief ... copy constructor
NuTo::DeformationGradient1D::DeformationGradient1D(const DeformationGradient1D& rOther)
{
	mDeformationGradient = rOther.mDeformationGradient;
}

// get number of components
unsigned int NuTo::DeformationGradient1D::GetNumberOfComponents() const
{
    return 1;
}

// get deformation gradient
const NuTo::DeformationGradient1D& NuTo::DeformationGradient1D::GetDeformationGradient1D() const
{
    return *this;
}

/*// get deformation gradient
void NuTo::DeformationGradient1D::GetDeformationGradient(NuTo::DeformationGradient1D& rDeformationGradient) const
{
    rDeformationGradient = NuTo::DeformationGradient1D(*this);
}

// get deformation gradient
void NuTo::DeformationGradient1D::GetDeformationGradient(NuTo::DeformationGradient2D& rDeformationGradient) const
{
    rDeformationGradient = NuTo::DeformationGradient2D(*this);
}


void NuTo::DeformationGradient1D::GetDeformationGradient(NuTo::DeformationGradient3D& rDeformationGradient) const
{
	rDeformationGradient = NuTo::DeformationGradient3D(*this);
}
*/

// set deformation gradient
void NuTo::DeformationGradient1D::SetDeformationGradient1D(const double* rDeformationGradient)
{
    this->mDeformationGradient = rDeformationGradient[0];
}

// calculate engineering strain
void NuTo::DeformationGradient1D::GetEngineeringStrain(NuTo::EngineeringStrain1D& rEngineeringStrain) const
{
    rEngineeringStrain.mEngineeringStrain = mDeformationGradient -1.;
}

void NuTo::DeformationGradient1D::GetEngineeringStrain(NuTo::EngineeringStrain2D& rEngineeringStrain) const
{
	rEngineeringStrain.mEngineeringStrain[0] = mDeformationGradient -1.;
	rEngineeringStrain.mEngineeringStrain[1] = 0.;
	rEngineeringStrain.mEngineeringStrain[2] = 0.;
}

void NuTo::DeformationGradient1D::GetEngineeringStrain(NuTo::EngineeringStrain3D& rEngineeringStrain) const
{
	rEngineeringStrain.mEngineeringStrain[0] = mDeformationGradient -1.;
	rEngineeringStrain.mEngineeringStrain[1] = 0.;
	rEngineeringStrain.mEngineeringStrain[2] = 0.;
	rEngineeringStrain.mEngineeringStrain[3] = 0.;
	rEngineeringStrain.mEngineeringStrain[4] = 0.;
	rEngineeringStrain.mEngineeringStrain[5] = 0.;
}

// calculate Green strain
void NuTo::DeformationGradient1D::GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain1D& rGreenLagrangeStrain) const
{
	throw MechanicsException("[NuTo::DeformationGradient1D::GetGreenLagrangeStrain] to be implemented.");
}

void NuTo::DeformationGradient1D::GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain2D& rGreenLagrangeStrain) const
{
	throw MechanicsException("[NuTo::DeformationGradient1D::GetGreenLagrangeStrain] to be implemented.");
}

void NuTo::DeformationGradient1D::GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain3D& rGreenLagrangeStrain) const
{
	throw MechanicsException("[NuTo::DeformationGradient1D::GetGreenLagrangeStrain] to be implemented.");
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::DeformationGradient1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::DeformationGradient1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize DeformationGradient1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveInputBase)
       & BOOST_SERIALIZATION_NVP(mDeformationGradient);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize DeformationGradient1D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

