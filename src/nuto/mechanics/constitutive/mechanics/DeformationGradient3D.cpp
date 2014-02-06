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
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"

#include "nuto/mechanics/MechanicsException.h"

// constructor
NuTo::DeformationGradient3D::DeformationGradient3D(): ConstitutiveInputBase::ConstitutiveInputBase()
{
    mDeformationGradient[0] = 0.0;
    mDeformationGradient[1] = 0.0;
    mDeformationGradient[2] = 0.0;
    mDeformationGradient[3] = 0.0;
    mDeformationGradient[4] = 0.0;
    mDeformationGradient[5] = 0.0;
    mDeformationGradient[6] = 0.0;
    mDeformationGradient[7] = 0.0;
    mDeformationGradient[8] = 0.0;
}

//! @brief ... copy constructor
NuTo::DeformationGradient3D::DeformationGradient3D(const DeformationGradient3D& rOther)
{
	mDeformationGradient[0] = rOther.mDeformationGradient[0];
    mDeformationGradient[1] = rOther.mDeformationGradient[1];
    mDeformationGradient[2] = rOther.mDeformationGradient[2];
    mDeformationGradient[3] = rOther.mDeformationGradient[3];
    mDeformationGradient[4] = rOther.mDeformationGradient[4];
    mDeformationGradient[5] = rOther.mDeformationGradient[5];
    mDeformationGradient[6] = rOther.mDeformationGradient[6];
    mDeformationGradient[7] = rOther.mDeformationGradient[7];
    mDeformationGradient[8] = rOther.mDeformationGradient[8];
}

// get number of components
unsigned int NuTo::DeformationGradient3D::GetNumberOfComponents() const
{
    return 9;
}

// get deformation gradient
const NuTo::DeformationGradient3D& NuTo::DeformationGradient3D::GetDeformationGradient3D() const
{
    return *this;
}

// set deformation gradient
void NuTo::DeformationGradient3D::SetDeformationGradient3D(const double* rDeformationGradient)
{
    for (unsigned int count = 0; count < 9; count++)
    {
        this->mDeformationGradient[count] = rDeformationGradient[count];
    }
}

// info routine
void NuTo::DeformationGradient3D::Info() const
{
    std::cout << "    components of deformation gradient (vector notation, columns): " << std::endl
              << mDeformationGradient[0] << ", " << mDeformationGradient[3] << ", " << mDeformationGradient[6] << ", " << std::endl
              << mDeformationGradient[1] << ", " << mDeformationGradient[4] << ", " << mDeformationGradient[7] << ", " << std::endl
              << mDeformationGradient[2] << ", " << mDeformationGradient[5] << ", " << mDeformationGradient[8] << std::endl;
}

// calculate engineering strain
void NuTo::DeformationGradient3D::GetEngineeringStrain(NuTo::EngineeringStrain3D& rEngineeringStrain) const
{
	rEngineeringStrain = NuTo::EngineeringStrain3D(*this);
}

// calculate Green strain
void NuTo::DeformationGradient3D::GetGreenLagrangeStrain(NuTo::GreenLagrangeStrain3D& rGreenLagrangeStrain) const
{
	rGreenLagrangeStrain = NuTo::GreenLagrangeStrain3D(*this);
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::DeformationGradient3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::DeformationGradient3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::DeformationGradient3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize DeformationGradient3D" << std::endl;
#endif
   ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveInputBase)
      & BOOST_SERIALIZATION_NVP(mDeformationGradient);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize DeformationGradient3D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

