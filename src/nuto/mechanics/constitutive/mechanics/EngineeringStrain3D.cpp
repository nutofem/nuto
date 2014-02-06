// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"

NuTo::EngineeringStrain3D::EngineeringStrain3D() : ConstitutiveOutputBase::ConstitutiveOutputBase(),ConstitutiveInputBase::ConstitutiveInputBase(), FullVector<double,6>()
{
	(*this)[0] = 0.0;
	(*this)[1] = 0.0;
	(*this)[2] = 0.0;
	(*this)[3] = 0.0;
	(*this)[4] = 0.0;
	(*this)[5] = 0.0;
}

NuTo::EngineeringStrain3D::EngineeringStrain3D(const DeformationGradient3D& rDeformationGradient)
{
	(*this)[0] = rDeformationGradient.mDeformationGradient[0] -1.;
	(*this)[1] = rDeformationGradient.mDeformationGradient[4] -1.;
	(*this)[2] = rDeformationGradient.mDeformationGradient[8] -1.;
	(*this)[3] = rDeformationGradient.mDeformationGradient[1]+rDeformationGradient.mDeformationGradient[3];
	(*this)[4] = rDeformationGradient.mDeformationGradient[5]+rDeformationGradient.mDeformationGradient[7];
	(*this)[5] = rDeformationGradient.mDeformationGradient[2]+rDeformationGradient.mDeformationGradient[6];
}

//! @brief ... get number of strain components
//! @return ... number of strain components
unsigned int NuTo::EngineeringStrain3D::GetNumberOfComponents() const
{
    return 6;
}

//! @brief ... get Engineering Strain
//! @return ... Engineering Strain (exx,eyy,ezz,gxy,gyz,gzx)
//! @sa mDeformationGradient
const double* NuTo::EngineeringStrain3D::GetData() const
{
    return data();
}

//! @brief ... get Engineering Strain
//! @return ... Engineering Strain (exx,eyy,ezz,gxy,gyz,gzx)
//! @sa mDeformationGradient
void NuTo::EngineeringStrain3D::SetData(const double rData[6])
{
	(*this)[0] = rData[0];
	(*this)[1] = rData[1];
	(*this)[2] = rData[2];
	(*this)[3] = rData[3];
	(*this)[4] = rData[4];
	(*this)[5] = rData[5];
}

//! @brief ... calculates the norm of the stress tensor in 3D case
#define sqrt_2div3 0.81649658
double NuTo::EngineeringStrain3D::Norm() const
{
	NuTo::FullVector<double,6> Strain;
	Strain(0) = (*this)[0];  // Strain(0) = eps_11
	Strain(1) = (*this)[1];  // Strain(1) = eps_22
	Strain(2) = (*this)[2];  // Strain(2) = eps_33
	Strain(3) = (*this)[3];  // Strain(3) = 2*eps_23
	Strain(4) = (*this)[4];  // Strain(4) = 2*eps_13
	Strain(5) = (*this)[5];  // Strain(5) = 2*eps_12

    //*******************************************************************
    //*    NorM = sqrt ( (2/3) * Strain : Strain)                                                *
    //*******************************************************************
    double invariante_2 = Strain(0)*Strain(0) + Strain(1)*Strain(1) + Strain(2)*Strain(2) +
    		0.5*(Strain(3)*Strain(3) + Strain(4)*Strain(4) + Strain(5)*Strain(5));

    return sqrt_2div3*std::sqrt(invariante_2);
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::EngineeringStrain3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStrain3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::EngineeringStrain3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize EngineeringStrain3D" << std::endl;
#endif
    ar &  BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveOutputBase);
    ar & boost::serialization::make_nvp ("EngineeringStrain3DEigen",boost::serialization::base_object< FullVector<double,6> > ( *this ));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize EngineeringStrain3D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
