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
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"

// constructor
NuTo::DeformationGradient2D::DeformationGradient2D()
{
    this->mDeformationGradient[0] = 0.0;
    this->mDeformationGradient[1] = 0.0;
    this->mDeformationGradient[2] = 0.0;
    this->mDeformationGradient[3] = 0.0;
}

//! @brief ... copy constructor
NuTo::DeformationGradient2D::DeformationGradient2D(const DeformationGradient1D& rOther)
{
	mDeformationGradient[0] = rOther.GetDeformationGradient1D()[0];
    mDeformationGradient[1] = 0.0;
    mDeformationGradient[2] = 0.0;
    mDeformationGradient[3] = 0.0;
}

//! @brief ... copy constructor
NuTo::DeformationGradient2D::DeformationGradient2D(const DeformationGradient2D& rOther)
{
	mDeformationGradient[0] = rOther.mDeformationGradient[0];
    mDeformationGradient[1] = rOther.mDeformationGradient[1];
    mDeformationGradient[2] = rOther.mDeformationGradient[2];
    mDeformationGradient[3] = rOther.mDeformationGradient[3];
}


// get number of components
unsigned int NuTo::DeformationGradient2D::GetNumberOfComponents() const
{
    return 4;
}

// get deformation gradient
const double* NuTo::DeformationGradient2D::GetDeformationGradient2D() const
{
    return this->mDeformationGradient;
}

void NuTo::DeformationGradient2D::GetDeformationGradient(NuTo::DeformationGradient2D& rDeformationGradient) const
{
    rDeformationGradient = NuTo::DeformationGradient2D(*this);
}

void NuTo::DeformationGradient2D::GetDeformationGradient(NuTo::DeformationGradient3D& rDeformationGradient) const
{
	rDeformationGradient = NuTo::DeformationGradient3D(*this);
}

// set deformation gradient
void NuTo::DeformationGradient2D::SetDeformationGradient2D(const double* rDeformationGradient)
{
    this->mDeformationGradient[0] = rDeformationGradient[0];
    this->mDeformationGradient[1] = rDeformationGradient[1];
    this->mDeformationGradient[2] = rDeformationGradient[2];
    this->mDeformationGradient[3] = rDeformationGradient[3];
}

// calculate engineering strain
void NuTo::DeformationGradient2D::GetEngineeringStrain(NuTo::EngineeringStrain2D& rEngineeringStrain) const
{
	rEngineeringStrain =  NuTo::EngineeringStrain2D(*this);
}

void NuTo::DeformationGradient2D::GetEngineeringStrain(NuTo::EngineeringStrain3D& rEngineeringStrain) const
{
	rEngineeringStrain =  NuTo::EngineeringStrain3D(*this);}

// calculate Green strain
void NuTo::DeformationGradient2D::GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain2D& rGreenLagrangeStrain) const
{
	rGreenLagrangeStrain = NuTo::GreenLagrangeStrain2D(*this);
}

void NuTo::DeformationGradient2D::GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain3D& rGreenLagrangeStrain) const
{
	rGreenLagrangeStrain = NuTo::GreenLagrangeStrain3D(*this);
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::DeformationGradient2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::DeformationGradient2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize DeformationGradient2D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_NVP(mDeformationGradient);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize DeformationGradient2D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

