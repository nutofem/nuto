// $ld: $ 
// NonlocalDamagePlasticity.cpp
// created Apr 26, 2010 by Joerg F. Unger
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <eigen2/Eigen/LU>
#include <eigen2/Eigen/Array>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal3x3.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal6x6.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataNonlocalDamagePlasticity3D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalDamagePlasticity.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"
#include "nuto/mechanics/elements/ElementWithDataBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

#define sqrt3 1.732050808
//#define ENABLE_DEBUG

NuTo::NonlocalDamagePlasticity::NonlocalDamagePlasticity() : ConstitutiveEngineeringStressStrain()
{
    std::cout << "[NuTo::NonlocalDamagePlasticity::NonlocalDamagePlasticity]" << std::endl;
    mE = 0.;
    mNu = 0.;
    mNonlocalRadius = 1.;
    mTensileStrength = 0.;
    mCompressiveStrength = 0.;
    mBiaxialCompressiveStrength = 0.;
    mFractureEnergy = 0.;
    mYieldSurface = Constitutive::eNonlocalDamageYieldSurface::COMBINED_ROUNDED;
    mDamage = false;
    SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::NonlocalDamagePlasticity::serialize(Archive & ar, const unsigned int version)
    {
        std::cout << "start serialization of linear elastic" << std::endl;
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveEngineeringStressStrain)
           & BOOST_SERIALIZATION_NVP(mE)
           & BOOST_SERIALIZATION_NVP(mNu)
           & BOOST_SERIALIZATION_NVP(mNonlocalRadius)
           & BOOST_SERIALIZATION_NVP(mTensileStrength)
           & BOOST_SERIALIZATION_NVP(mCompressiveStrength)
           & BOOST_SERIALIZATION_NVP(mBiaxialCompressiveStrength)
           & BOOST_SERIALIZATION_NVP(mFractureEnergy)
           & BOOST_SERIALIZATION_NVP(mYieldSurface)
           & BOOST_SERIALIZATION_NVP(mDamage);
        std::cout << "finish serialization of linear elastic" << std::endl;
    }
#endif // ENABLE_SERIALIZATION

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::NonlocalDamagePlasticity::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    //TODO check this routine if it is necessary
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::NonlocalDamagePlasticity::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    //TODO check this routine if it is necessary
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::NonlocalDamagePlasticity::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    //TODO check this routine if it is necessary
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
              const DeformationGradient1D& rDeformationGradient, EngineeringStress1D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Material model not implemented for 1D.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... engineering stress
void NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
          const DeformationGradient1D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Material model not implemented for 1D.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
              const DeformationGradient2D& rDeformationGradient, EngineeringStress2D& rEngineeringStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    assert(rElement->GetSection()!=0);
    if (rElement->GetSection()->GetType()==Section::PLANE_STRAIN)
    {
        ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain *oldStaticData = (rElement->GetStaticData(rIp))->AsNonlocalDamagePlasticity2DPlaneStrain();
        double rNewEqPlasticStrain;
        Eigen::Matrix<double,4,4> rdEpsilonPdEpsilon;
        Eigen::Matrix<double,4,1> rNewEpsilonP;

        // no new history variables, since no update has to be performed
        this->ReturnMapping2D(
                engineeringStrain,
                oldStaticData->mEpsilonP,
                oldStaticData->mPrevStrain,
                oldStaticData->mKappa,
                rNewEpsilonP,
                rNewEqPlasticStrain,
                rdEpsilonPdEpsilon);
/*
        const EngineeringStrain2D& rStrain,
        double l_eq_plane,
        double l_eq_circ,
        const double rPrevPlasticStrain[4],
        const double rPrevTotalStrain[3],
        double rPrevEqPlasticStrain,
        Eigen::Matrix<double,4,1>& rEpsilonP,
        double& rNewEqPlasticStrain,
        Eigen::Matrix<double,4,4>& rdEpsilonPdEpsilon
*/
        // calculate Engineering stress
        //rEngineeringStress.mEngineeringStress(0) = C11 * engineeringStrain.mEngineeringStrain(0) + C12 * engineeringStrain.mEngineeringStrain(1);
        //rEngineeringStress.mEngineeringStress(1) = C11 * engineeringStrain.mEngineeringStrain(1) + C12 * engineeringStrain.mEngineeringStrain(0);
        //rEngineeringStress.mEngineeringStress(2) = C33 * engineeringStrain.mEngineeringStrain(2) ;

    }
    else
    {
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Plane stress is not implemented.");
    }
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
              const DeformationGradient2D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
        //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    assert(rElement->GetSection()!=0);
    if (rElement->GetSection()->GetType()==Section::PLANE_STRAIN)
    {
        // calculate coefficients of the material matrix
        double C11, C12, C33;
        this->CalculateCoefficients3D(C11, C12, C33);

        // calculate Engineering stress
        rEngineeringStress.mEngineeringStress[0] = C11 * engineeringStrain.mEngineeringStrain[0] + C12 * engineeringStrain.mEngineeringStrain[1];
        rEngineeringStress.mEngineeringStress[1] = C11 * engineeringStrain.mEngineeringStrain[1] + C12 * engineeringStrain.mEngineeringStrain[0];
        rEngineeringStress.mEngineeringStress[2] = C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[1]);
        rEngineeringStress.mEngineeringStress[3] = C33 * engineeringStrain.mEngineeringStrain[2] ;
        rEngineeringStress.mEngineeringStress[4] = 0.;
        rEngineeringStress.mEngineeringStress[5] = 0.;
    }
    else
    {
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Plane stress is to be implemented.");
    }
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
              const DeformationGradient3D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain3D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    // calculate coefficients of the material matrix
    double C11, C12, C44;
    this->CalculateCoefficients3D(C11, C12, C44);

    // calculate Engineering stress
    rEngineeringStress.mEngineeringStress[0] = C11 * engineeringStrain.mEngineeringStrain[0] + C12 * (engineeringStrain.mEngineeringStrain[1]+engineeringStrain.mEngineeringStrain[2]);
    rEngineeringStress.mEngineeringStress[1] = C11 * engineeringStrain.mEngineeringStrain[1] + C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[2]);
    rEngineeringStress.mEngineeringStress[2] = C11 * engineeringStrain.mEngineeringStrain[2] + C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[1]);
    rEngineeringStress.mEngineeringStress[3] = C44 * engineeringStrain.mEngineeringStrain[3] ;
    rEngineeringStress.mEngineeringStress[4] = C44 * engineeringStrain.mEngineeringStrain[4] ;
    rEngineeringStress.mEngineeringStress[5] = C44 * engineeringStrain.mEngineeringStrain[5] ;
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient,
        ConstitutiveTangentLocal1x1& rTangent) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // store tangent at the output object
    rTangent.mTangent = mE;
    rTangent.SetSymmetry(true);
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient,
        ConstitutiveTangentLocal3x3& rTangent) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    const SectionBase* theSection(rElement->GetSection());
    if (theSection==0)
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain] No section defined for element.");
    if (theSection->GetType()==Section::PLANE_STRAIN)
    {
        // calculate coefficients of the material matrix
        double C11, C12, C33;
        this->CalculateCoefficients3D(C11, C12, C33);

        // store tangent at the output object
         rTangent.mTangent[ 0] = C11;
         rTangent.mTangent[ 1] = C12;
         rTangent.mTangent[ 2] = 0;

         rTangent.mTangent[ 3] = C12;
         rTangent.mTangent[ 4] = C11;
         rTangent.mTangent[ 5] = 0;

         rTangent.mTangent[ 6] = 0.;
         rTangent.mTangent[ 7] = 0.;
         rTangent.mTangent[ 8] = C33;
    }
    else
    {
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetEngineeringStressFromEngineeringStrain] Plane stress is to be implemented.");
    }

    rTangent.SetSymmetry(true);
}


//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient,
        ConstitutiveTangentLocal6x6& rTangent) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // calculate coefficients of the material matrix
    double C11, C12, C44;
    this->CalculateCoefficients3D(C11, C12, C44);

    // store tangent at the output object
    rTangent.mTangent[ 0] = C11;
    rTangent.mTangent[ 1] = C12;
    rTangent.mTangent[ 2] = C12;
    rTangent.mTangent[ 3] = 0;
    rTangent.mTangent[ 4] = 0;
    rTangent.mTangent[ 5] = 0;

    rTangent.mTangent[ 6] = C12;
    rTangent.mTangent[ 7] = C11;
    rTangent.mTangent[ 8] = C12;
    rTangent.mTangent[ 9] = 0;
    rTangent.mTangent[10] = 0;
    rTangent.mTangent[11] = 0;

    rTangent.mTangent[12] = C12;
    rTangent.mTangent[13] = C12;
    rTangent.mTangent[14] = C11;
    rTangent.mTangent[15] = 0;
    rTangent.mTangent[16] = 0;
    rTangent.mTangent[17] = 0;

    rTangent.mTangent[18] = 0;
    rTangent.mTangent[19] = 0;
    rTangent.mTangent[20] = 0;
    rTangent.mTangent[21] = C44;
    rTangent.mTangent[22] = 0;
    rTangent.mTangent[23] = 0;

    rTangent.mTangent[24] = 0;
    rTangent.mTangent[25] = 0;
    rTangent.mTangent[26] = 0;
    rTangent.mTangent[27] = 0;
    rTangent.mTangent[28] = C44;
    rTangent.mTangent[29] = 0;

    rTangent.mTangent[30] = 0;
    rTangent.mTangent[31] = 0;
    rTangent.mTangent[32] = 0;
    rTangent.mTangent[33] = 0;
    rTangent.mTangent[34] = 0;
    rTangent.mTangent[35] = C44;

    rTangent.SetSymmetry(true);
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::NonlocalDamagePlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    //no static data required -> empty routine
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::NonlocalDamagePlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    //no static data required -> empty routine
}


//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::NonlocalDamagePlasticity::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
       //no static data required -> empty routine
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::NonlocalDamagePlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain1D(
        const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::NonlocalDamagePlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Nonlocal damage plasticity model not implemented for 1D.");
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::NonlocalDamagePlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain2D(
        const ElementBase* rElement) const
{
    if (rElement->GetSection()==0)
    	throw MechanicsException("[NuTo::NonlocalDamagePlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Section required to distinguish between plane stress and plane strain and thickness information.");
	if (rElement->GetSection()->GetType()==NuTo::Section::PLANE_STRESS)
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Nonlocal damage plasticity model not implemented for plane stress.");
    else
        return new ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain();
}


//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::NonlocalDamagePlasticity::AllocateStaticDataEngineeringStress_EngineeringStrain3D(
        const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataNonlocalDamagePlasticity3D();
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // calculate engineering strain
    EngineeringStrain1D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);
    return 0.5 * engineeringStrain.mEngineeringStrain * this->mE * engineeringStrain.mEngineeringStrain;
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    EngineeringStress2D engineeringStress;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    const SectionBase* theSection(rElement->GetSection());
    if (theSection==0)
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTangent_EngineeringStress_EngineeringStrain] No section defined for element.");
    if (theSection->GetType()==Section::PLANE_STRAIN)
    {
        // calculate coefficients of the material matrix
        double C11, C12, C33;
        this->CalculateCoefficients3D(C11, C12, C33);

        // calculate Engineering stress
        engineeringStress.mEngineeringStress[0] = C11 * engineeringStrain.mEngineeringStrain[0] + C12 * engineeringStrain.mEngineeringStrain[1];
        engineeringStress.mEngineeringStress[1] = C11 * engineeringStrain.mEngineeringStrain[1] + C12 * engineeringStrain.mEngineeringStrain[0];
        engineeringStress.mEngineeringStress[2] = C33 * engineeringStrain.mEngineeringStrain[2] ;
    }
    else
    {
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain] Plane stress is to be implemented.");
    }
    return 0.5*(
            engineeringStrain.mEngineeringStrain[0]*engineeringStress.mEngineeringStress[0]
           +engineeringStrain.mEngineeringStrain[1]*engineeringStress.mEngineeringStress[1]
           +engineeringStrain.mEngineeringStrain[2]*engineeringStress.mEngineeringStress[2]);
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain3D engineeringStrain;
    EngineeringStress3D engineeringStress;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    // calculate coefficients of the material matrix
    double C11, C12, C44;
    this->CalculateCoefficients3D(C11, C12, C44);

    // calculate Engineering stress
    engineeringStress.mEngineeringStress[0] = C11 * engineeringStrain.mEngineeringStrain[0] + C12 * (engineeringStrain.mEngineeringStrain[1]+engineeringStrain.mEngineeringStrain[2]);
    engineeringStress.mEngineeringStress[1] = C11 * engineeringStrain.mEngineeringStrain[1] + C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[2]);
    engineeringStress.mEngineeringStress[2] = C11 * engineeringStrain.mEngineeringStrain[2] + C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[1]);
    engineeringStress.mEngineeringStress[3] = C44 * engineeringStrain.mEngineeringStrain[3] ;
    engineeringStress.mEngineeringStress[4] = C44 * engineeringStrain.mEngineeringStrain[4] ;
    engineeringStress.mEngineeringStress[5] = C44 * engineeringStrain.mEngineeringStrain[5] ;

    return 0.5*(
            engineeringStrain.mEngineeringStrain[0]*engineeringStress.mEngineeringStress[0]
           +engineeringStrain.mEngineeringStrain[1]*engineeringStress.mEngineeringStress[1]
           +engineeringStrain.mEngineeringStrain[2]*engineeringStress.mEngineeringStress[2]
           +engineeringStrain.mEngineeringStrain[3]*engineeringStress.mEngineeringStress[3]
           +engineeringStrain.mEngineeringStrain[4]*engineeringStress.mEngineeringStress[4]
           +engineeringStrain.mEngineeringStrain[5]*engineeringStress.mEngineeringStress[5]);
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    return GetTotalEnergy_EngineeringStress_EngineeringStrain(rElement, rIp, rDeformationGradient);
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    return GetTotalEnergy_EngineeringStress_EngineeringStrain(rElement, rIp, rDeformationGradient);
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::NonlocalDamagePlasticity::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
    return GetTotalEnergy_EngineeringStress_EngineeringStrain(rElement, rIp, rDeformationGradient);
}


//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStrain1D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain] this method is only required for \
NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain, which is reimplemented for Linear elastic -- consequently, this metods is not required.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStrain2D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain] this method is only required for \
NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain, which is reimplemented for Linear elastic -- consequently, this metods is not required.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::NonlocalDamagePlasticity::GetDeltaElasticEngineeringStrain] this method is only required for \
NuTo::ConstitutiveEngineeringStressStrain::GetTotalEnergy_EngineeringStress_EngineeringStrain, which is reimplemented for Linear elastic -- consequently, this metods is not required.");
}

// calculate coefficients of the material matrix
void NuTo::NonlocalDamagePlasticity::CalculateCoefficients3D(double& C11, double& C12, double& C44) const
{
    double factor = this->mE/((1.0 + this->mNu) * (1.0 - 2.0 * this->mNu));
    C11 = factor * (1.0 - this->mNu);
    C12 = factor * this->mNu;
    C44 = this->mE/(2*(1.0 + this->mNu));
}

// parameters /////////////////////////////////////////////////////////////
//! @brief ... get Young's modulus
//! @return ... Young's modulus
double NuTo::NonlocalDamagePlasticity::GetYoungsModulus() const
{
    return mE;
}


//! @brief ... set Young's modulus
//! @param rE ... Young's modulus
void NuTo::NonlocalDamagePlasticity::SetYoungsModulus(double rE)
{
    this->CheckYoungsModulus(rE);
    this->mE = rE;
    this->SetParametersValid();
}


//! @brief ... get Poisson's ratio
//! @return ... Poisson's ratio
double NuTo::NonlocalDamagePlasticity::GetPoissonsRatio() const
{
    return mNu;
}

//! @brief ... set Poisson's ratio
//! @param rNu ... Poisson's ratio
void NuTo::NonlocalDamagePlasticity::SetPoissonsRatio(double rNu)
{
    this->CheckPoissonsRatio(rNu);
    this->mNu = rNu;
    this->SetParametersValid();
}

//! @brief ... get nonlocal radius
//! @return ... nonlocal radius
double NuTo::NonlocalDamagePlasticity::GetNonlocalRadius() const
{
    return mNonlocalRadius;
}

//! @brief ... set nonlocal radius
//! @param rRadius...  nonlocal radius
void NuTo::NonlocalDamagePlasticity::SetNonlocalRadius(double rNonlocalRadius)
{
    this->CheckNonlocalRadius(rNonlocalRadius);
    this->mNonlocalRadius = rNonlocalRadius;
    this->SetParametersValid();
}
//! @brief ... get tensile strength
//! @return ... tensile strength
double NuTo::NonlocalDamagePlasticity::GetTensileStrength() const
{
    return mTensileStrength;
}

//! @brief ... set tensile strength
//! @param rTensileStrength...  tensile strength
void NuTo::NonlocalDamagePlasticity::SetTensileStrength(double rTensileStrength)
{
    this->CheckNonlocalRadius(rTensileStrength);
    this->mTensileStrength = rTensileStrength;
    this->SetParametersValid();
}

//! @brief ... get compressive strength
//! @return ... compressive strength
double NuTo::NonlocalDamagePlasticity::GetCompressiveStrength() const
{
    return mCompressiveStrength;
}

//! @brief ... set compressive strength
//! @param rCompressiveStrength...  compressive strength
void NuTo::NonlocalDamagePlasticity::SetCompressiveStrength(double rCompressiveStrength)
{
    this->CheckNonlocalRadius(rCompressiveStrength);
    this->mCompressiveStrength = rCompressiveStrength;
    this->SetParametersValid();
}

//! @brief ... get biaxial compressive strength
//! @return ... biaxial compressive strength
double NuTo::NonlocalDamagePlasticity::GetBiaxialCompressiveStrength() const
{
    return mBiaxialCompressiveStrength;
}

//! @brief ... set biaxial compressive strength
//! @param rBiaxialCompressiveStrength...  biaxial compressive strength
void NuTo::NonlocalDamagePlasticity::SetBiaxialCompressiveStrength(double rBiaxialCompressiveStrength)
{
    this->CheckNonlocalRadius(rBiaxialCompressiveStrength);
    this->mBiaxialCompressiveStrength = rBiaxialCompressiveStrength;
    this->SetParametersValid();
}

//! @brief ... get fracture energy
//! @return ... fracture energy
double NuTo::NonlocalDamagePlasticity::GetFractureEnergy() const
{
    return mFractureEnergy;
}

//! @brief ... set fracture energy
//! @param rFractureEnergy... fracture energy
void NuTo::NonlocalDamagePlasticity::SetFractureEnergy(double rFractureEnergy)
{
    this->CheckNonlocalRadius(rFractureEnergy);
    this->mFractureEnergy = rFractureEnergy;
    this->SetParametersValid();
}
///////////////////////////////////////////////////////////////////////////


//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::NonlocalDamagePlasticity::GetType() const
{
    return NuTo::Constitutive::NONLOCAL_DAMAGE_PLASTICITY;
}


//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::NonlocalDamagePlasticity::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::BRICK8N:
        return true;
    case NuTo::Element::PLANE2D3N:
        return true;
    case NuTo::Element::PLANE2D4N:
        return true;
    case NuTo::Element::PLANE2D6N:
        return true;
    case NuTo::Element::TETRAHEDRON4N:
        return true;
    case NuTo::Element::TETRAHEDRON10N:
        return true;
    case NuTo::Element::TRUSS1D2N:
        return true;
    case NuTo::Element::TRUSS1D3N:
        return true;
    default:
        return false;
    }
}

//! @brief ... check if Young's modulus is positive
//! @param rE ... Young's modulus
void NuTo::NonlocalDamagePlasticity::CheckYoungsModulus(double rE) const
{
    if (rE <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckYoungsModulus] The Young's modulus must be a positive value.");
    }
}

//! @brief ... check if Poisson's ratio is valid \f$ (-1.0 < \nu < 0.5) \f$
//! @param rE ... Poisson's ratio
void NuTo::NonlocalDamagePlasticity::CheckPoissonsRatio(double rNu) const
{
    if (rNu <= -1.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckPoissonsRatio] Poisson's ratio must be greater or equal to -1.0.");
    }
    if (rNu >= 0.5)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckPoissonsRatio] Poisson's ratio must be smaller or equal to 0.5.");
    }
}

//! @brief ... check if the nonlocal radius is positive
//! @param rRadius ... nonlocal radius
void NuTo::NonlocalDamagePlasticity::CheckNonlocalRadius(double rRadius) const
{
    if (rRadius <= 0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckNonlocalRadius] Nonlocal radius must be positive.");
    }
}
//! @brief ... check if tensile strength is positive
//! @param rTensileStrength ... nonlocal radius
void NuTo::NonlocalDamagePlasticity::CheckTensileStrength(double rTensileStrength) const
{
    if (rTensileStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckTensileStrength] The tensile strength must be a positive value.");
    }
}

//! @brief ... check if compressive strength is positive
//! @param rRadius ... compressive strength
void NuTo::NonlocalDamagePlasticity::CheckCompressiveStrength(double rCompressiveStrength) const
{
    if (rCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckCompressiveStrength] The compressive strength must be a positive value.");
    }
}

//! @brief ... check if biaxial compressive strength is positive
//! @param rBiaxialCompressiveStrength ... biaxial compressive strength
void NuTo::NonlocalDamagePlasticity::CheckBiaxialCompressiveStrength(double rBiaxialCompressiveStrength) const
{
    if (rBiaxialCompressiveStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckBiaxialCompressiveStrength] The biaxial compressive strength must be a positive value.");
    }
    if (rBiaxialCompressiveStrength <= mCompressiveStrength)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckBiaxialCompressiveStrength] The biaxial compressive strength must be higher than the uniaxial compressive strength.");
    }
}

//! @brief ... check if fracture energy is positive
//! @param rFractureEnergy ... fracture energy
void NuTo::NonlocalDamagePlasticity::CheckFractureEnergy(double rFractureEnergy) const
{
    if (rFractureEnergy <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::NonlocalDamagePlasticity::CheckFractureEnergy] The fracture energy must be a positive value.");
    }
}


//! @brief ... calculate the length of the element in plane coordinates (square root of area)
void NuTo::NonlocalDamagePlasticity::CalculateEquivalentLength(const ElementBase* rElement, double& l_eq_plane, double& l_eq_circ) const
{
    double l_element(sqrt(rElement->CalculateArea()*0.5));
    l_eq_plane = mNonlocalRadius;
    l_eq_circ  = mNonlocalRadius*rElement->GetSection()->GetThickness()/l_element;
}

#define ACTIVE true
#define INACTIVE false
#define toleranceResidual 1e-8      //tolerance to decide whether the Newton iteration has converged
#define toleranceYieldSurface 1e-10 //tolerance whether a point is on the yield surface or not (multiplied by the YoungsModulus)
#define toleranceDeterminant 1e-12  //tolerance to decide if a matrix is not invertible (only required in the debugging version)
#define maxSteps 25                 //maximum number of Newton iterations, until it is decided that there is no convergence and a cutback is performed
#define minCutbackFactor 1e-3       //minimum cutback factor for the application of the total strain in steps
#define minCutbackFactorLS 1e-3     //minimum cutback factor used for the linesearch in the Newton iteration
//! @brief ... performs the return mapping procedure for the plasticity model
//! @param rStrain              ... current total strain
//! @param rPrevPlasticStrain   ... previous plastic strain (history variable)
//! @param rPrevTotalStrain     ... previous total strain (history variable)
//! @param rPrevEqPlasticStrain ... previous equiavalente plastic strain (history variable)
//! @param rEpsilonP            ... new plastic strain after return mapping
//! @param rEqPlasticStrain     ... new equivalente olastic strain after return mapping
//! @param rdEpsilonPdEpsilon   ... new derivative of current plastic strain with respect to the total strain
void NuTo::NonlocalDamagePlasticity::ReturnMapping2D(
        const EngineeringStrain2D& rStrain,
        const double rPrevPlasticStrain[4],
        const double rPrevTotalStrain[3],
        double rPrevEqPlasticStrain,
        Eigen::Matrix<double,4,1>& rEpsilonP,
        double& rEqPlasticStrain,
        Eigen::Matrix<double,4,4>& rdEpsilonPdEpsilon)const
{
    double e_mod = mE; //modify that one in the case of random fields
    double   nu  = mNu;
    double f_ct  = mTensileStrength;
    double f_c1  = mCompressiveStrength;
    double f_c2  = mBiaxialCompressiveStrength;

    assert(f_c2>f_c1);

    //******************************************************************
    //*    F_BETA:    required by DRUCKER-PRAGER yield surface               *
    //******************************************************************
    double BETA = sqrt3*(f_c2-f_c1) / (2*f_c2-f_c1);
    double H_P  = f_c2*f_c1 / (sqrt3*(2*f_c2-f_c1));

    //! @brief strain currently solved for the plastic strains, in general equal to  rStrain, but if applied in steps, it's smaller
    Eigen::Vector4d curTotalStrain;
    //! @brief previous plastic strain, either from previous equilibrium (static data) or if applied in steps, previous converged state
    Eigen::Vector4d lastPlastStrain;
    //! @brief plastic strain in the line search
    Eigen::Vector4d epsilonPLS;
    //! @brief previous eq plastic strain, either from previous equilibrium (static data) or if applied in steps, previous converged state
    double lastEqPlasticStrain;
    //! @brief residual in the return mapping procedure
    Eigen::Vector4d residual;
    //! @brief residual in the return mapping procedure within linesearch
    Eigen::Vector4d residualLS;
    //! @brief full stress increment within one iteration of the return mapping, might be applied in steps in the succeeding line search
    Eigen::Vector4d deltaStress;
    //! @brief total strain increment between strain from previous static data and new total strain
    Eigen::Vector4d deltaStrain;
    //! @brief plastic strain increment within the linesearch
    Eigen::Vector4d deltaPlasticStrainLS;
    //! @brief trial stress of the first iteration
    Eigen::Vector4d initTrialStress;
    //! @brief trial stress in the line search
    Eigen::Vector4d stressLS;
    //! @brief trial stress
    Eigen::Vector4d trialStress;
    //! @brief elastic strain
    Eigen::Vector4d elasticStrain;
    //! @brief elastic strain in line search
    Eigen::Vector4d elasticStrainLS;
    //! @brief plastic multiplier
    Eigen::Matrix<double,2,1> deltaGamma;
    //! @brief plastic multiplier
    Eigen::Matrix<double,2,1> deltaGammaLS;
    //! @brief increment of plastic multiplier in return mapping procedure
    Eigen::Matrix<double,2,1> delta2Gamma;
    //! @brief yield condition
    Eigen::Matrix<double,2,1> yieldCondition;
    //! @brief yield condition in line search
    Eigen::Matrix<double,2,1> yieldConditionLS;
    //! @brief yield condition at the first iteration
    Eigen::Matrix<double,2,1> initYieldCondition;
    //! @brief flag that indicates if a yield function is active or not
    Eigen::Matrix<bool,2,1> yieldConditionFlag;
    //! @brief (dF/dsigma)T * Hessian * dF/dsigma
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matG;
    //! @brief ((dF/dsigma)T * Hessian * dF/dsigma )^-1
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matGInv;
    //! @brief algorithmic modulus DElastInv+delta_gamm*d2F/d2Sigma
    Eigen::Matrix<double,4,4> hessian;
    //! @brief temporary matrix
    Eigen::Matrix<double,4,4> tmpStiffness;
    //! @brief first derivatives of the yield functions
    std::vector<Eigen::Matrix<double,4,1> > dF_dsigma;
    //! @brief second derivatives of the yield functions
    std::vector<Eigen::Matrix<double,4,4> > d2F_d2sigma;
    //! @brief algorithmic modulus * dF_dsigma
    std::vector<Eigen::Matrix<double,4,1> > vectorN;
    //! @brief number of active yield functions
    int numActiveYieldFunctions;

    // for the application of strains in steps, calculate the total strain increment to be applied
    deltaStrain(0) = rStrain.mEngineeringStrain[0]-rPrevTotalStrain[0];
    deltaStrain(1) = rStrain.mEngineeringStrain[1]-rPrevTotalStrain[1];
    deltaStrain(2) = rStrain.mEngineeringStrain[2]-rPrevTotalStrain[2];
    deltaStrain(3) = 0;

    // initialize last plastic strain and last converged stress
    lastPlastStrain << rPrevPlasticStrain[0] , rPrevPlasticStrain[1] ,rPrevPlasticStrain[2] ,rPrevPlasticStrain[3];
	lastEqPlasticStrain = rPrevEqPlasticStrain;

    //*****************************************************************
    //                  elastic matrix generation                       *
    //*****************************************************************
    //! @brief elastic stiffness
    Eigen::Matrix<double,4,4> dElast;
    //! @brief inverse elastic stiffness
    Eigen::Matrix<double,4,4> dElastInv;

    {
		double factor = e_mod/((1.+nu)*(1.-2.*nu));
		double oneminusnufactor = (1-nu)*factor;
		double nufactor = nu*factor;

		dElast <<  oneminusnufactor , nufactor         , 0.              , nufactor ,
				   nufactor         , oneminusnufactor , 0.              , nufactor ,
				   0.               , 0.               , (0.5-nu)*factor , 0.,
				   nufactor         , nufactor         , 0.              ,oneminusnufactor;

		factor = 1./e_mod;
		double minusnufactor = -nu*factor;
		dElastInv   << factor        , minusnufactor , 0.               , minusnufactor,
					   minusnufactor , factor        , 0.               , minusnufactor,
					   0.            , 0.            , 2.*factor*(1.+nu), 0.,
					   minusnufactor , minusnufactor , 0                , factor;

    }
    //! @brief delta load factor for the previous iteration
    double deltaCutbackFactorExternal(1.);

    //! @brief current load factor (between 0 and 1) to apply the total strain increment in steps
    double cutbackFactorExternal(deltaCutbackFactorExternal);

    //! @brief flag to determine if the iteration is finished (converged at  cutbackFactorExternal=1)
    bool convergedExternal(false);

    int numberOfExternalCutbacks(0);
    int numberOfInternalIterations(0);
    int prevNumberOfInternalIterations(0);
    while (cutbackFactorExternal>minCutbackFactor && !convergedExternal)
    {
    	numberOfExternalCutbacks++;

        curTotalStrain(0) = rPrevTotalStrain[0]+cutbackFactorExternal*deltaStrain(0);
        curTotalStrain(1) = rPrevTotalStrain[1]+cutbackFactorExternal*deltaStrain(1);
        curTotalStrain(2) = rPrevTotalStrain[2]+cutbackFactorExternal*deltaStrain(2);
        curTotalStrain(3) = 0.;

#ifdef ENABLE_DEBUG
        std::cout << std::endl << "curTotalStrain" << curTotalStrain.transpose() << std::endl << std::endl;
#endif

        // checks the convergence of the Newton iteration for a prescribed current strain
        bool convergedInternal(false);
        try
        {
            //resize yield condition vector
        	int numYieldSurfaces=2;
			yieldCondition.setZero(numYieldSurfaces);

			//elastic strain, stress and d_matrix
			elasticStrain = curTotalStrain - lastPlastStrain;

			//TODO just use the upper part for new EigenVersion 3.0
			//trialStress = dElast.selfadjointView<Eigen::Upper>()*elasticStrain;
			initTrialStress = dElast*elasticStrain;
#ifdef ENABLE_DEBUG
			std::cout << "initTrialStress " << std::endl << initTrialStress.transpose() << std::endl << std::endl;
#endif

			//calculate yield condition
			//Drucker Prager
			initYieldCondition(0) = YieldSurfaceDruckerPrager2D(initTrialStress, BETA, H_P);

			//rounded Rankine
			initYieldCondition(1) = YieldSurfaceRankine2DRounded(initTrialStress, f_ct);

			if (initYieldCondition(0)<-toleranceYieldSurface*e_mod && initYieldCondition(1)<-toleranceYieldSurface*e_mod)
			{
				//*************************************************
				//*  thus we have elastic -------------> elastic  *
				//*************************************************
			    convergedInternal = true;
				rEpsilonP =  lastPlastStrain;
				rEqPlasticStrain = lastEqPlasticStrain;
				trialStress = initTrialStress;
				if (cutbackFactorExternal==1)
				{
					rdEpsilonPdEpsilon.setZero(4,4);
				}
#ifdef ENABLE_DEBUG
				std::cout << "linear elastic step" << std::endl << std::endl;
#endif
			}

			// perform return mapping
            dF_dsigma.resize(numYieldSurfaces);
            d2F_d2sigma.resize(numYieldSurfaces);

            // initialize plastic multiplier
			deltaGamma.setZero(numYieldSurfaces);
			delta2Gamma.setZero(numYieldSurfaces);
			yieldConditionFlag.setZero(numYieldSurfaces);

            for (int fixedYieldConditions=0; fixedYieldConditions<3 && convergedInternal==false; fixedYieldConditions++)
            {
    			switch(fixedYieldConditions)
				{
				case 0:
#ifdef ENABLE_DEBUG
					std::cout<< "pure Rankine" << std::endl;
#endif
					deltaGamma(1) = 0;
					if (initYieldCondition(1)<-toleranceYieldSurface*e_mod)
						continue;
					yieldConditionFlag(0) = INACTIVE;
					yieldConditionFlag(1) = ACTIVE;
					numActiveYieldFunctions = 1;
					break;
				case 1:
#ifdef ENABLE_DEBUG
					std::cout<< "combined" << std::endl;
#endif
					if (initYieldCondition(0)<-toleranceYieldSurface*e_mod || initYieldCondition(1)<-toleranceYieldSurface*e_mod)
						continue;
					yieldConditionFlag(0) = ACTIVE;
					yieldConditionFlag(1) = ACTIVE;
					deltaGamma(0) = 0;
					deltaGamma(1) = 0;
					numActiveYieldFunctions = 2;
					break;
				case 2:
#ifdef ENABLE_DEBUG
					std::cout<< "pure DP" << std::endl;
#endif
					if (initYieldCondition(0)<-toleranceYieldSurface*e_mod)
						continue;
					yieldConditionFlag(0) = ACTIVE;
					yieldConditionFlag(1) = INACTIVE;
					deltaGamma(0) = 0;
					numActiveYieldFunctions = 1;
					break;
				default:
					throw MechanicsException("[NuTo::NonlocalDamagePlasticity::ReturnMapping2D] programming error - should not happen.");
				}

            	for (int iteration = 0; iteration < maxSteps && convergedInternal==false; iteration++)
				{
            		numberOfInternalIterations++;
            		if (iteration==0)
					{
						rEpsilonP = lastPlastStrain;
						rEqPlasticStrain = lastEqPlasticStrain;
						yieldCondition= initYieldCondition;
						trialStress = initTrialStress;
					}
					else
					{
						//trial stress is the last stress state from the previous line search
						if (yieldConditionFlag(0)==ACTIVE)
							yieldCondition(0) = yieldConditionLS(0);
						else
							yieldCondition(0) = YieldSurfaceDruckerPrager2D(trialStress, BETA, H_P);

						if (yieldConditionFlag(1)==ACTIVE)
							yieldCondition(1) = yieldConditionLS(1);
						else
							yieldCondition(1) = YieldSurfaceRankine2DRounded(trialStress, f_ct);

					}
#ifdef ENABLE_DEBUG
                    std::cout << "trialStress " <<  std::endl << trialStress.transpose() << std::endl << std::endl;
					std::cout << "yieldCondition " <<  std::endl << yieldCondition.transpose() << std::endl << std::endl;
#endif

					// DP
					if (yieldConditionFlag(0)==ACTIVE)
					{
						if (!YieldSurfaceDruckerPrager2DDerivatives(dF_dsigma[0],&(d2F_d2sigma[0]),trialStress,BETA))
						{
							//no convergence, decrease line search step
							iteration = maxSteps;
							continue;
						}
#ifdef ENABLE_DEBUG
						std::cout << "dF_dsigma[0] (DP) " <<  std::endl << dF_dsigma[0].transpose() << std::endl << std::endl;
#endif
					}

					// Rounded Rankine
					if (yieldConditionFlag(1)==ACTIVE)
					{
						YieldSurfaceRankine2DRoundedDerivatives(dF_dsigma[1],&(d2F_d2sigma[1]),trialStress);
#ifdef ENABLE_DEBUG
						std::cout << "dF_dsigma[1] (Rankine)" <<  std::endl << dF_dsigma[1].transpose() << std::endl << std::endl;
#endif
					}


					//************************************************************************
					// residual
					//************************************************************************
					residual = lastPlastStrain-rEpsilonP;

					for (int count=0; count<numYieldSurfaces; count++)
					{
						if (yieldConditionFlag[count] == ACTIVE)
						{
							residual += deltaGamma(count)*dF_dsigma[count];
						}
					}

#ifdef ENABLE_DEBUG
                    std::cout << "residual " <<  std::endl << residual.transpose() << std::endl << std::endl;
#endif

					//this is just for scaling with a relative norm
					double epsilonEq = curTotalStrain.norm();
					double absResidual = residual.norm()/epsilonEq;

#ifdef ENABLE_DEBUG
                    std::cout << iteration <<" residual " << abs_residual << " yield conditionc " << yieldCondition.transpose() << std::endl << std::endl;
#endif

					// in case of PERFECT PLASTICITY [A] = hessian
					hessian = dElastInv;

					for (int count=0; count<numYieldSurfaces; count++)
					{
						if (yieldConditionFlag(count)==ACTIVE)
						{
							hessian+=deltaGamma(count)*d2F_d2sigma[count];
						}
					}

					assert(fabs(hessian.determinant())>toleranceDeterminant);
					hessian = hessian.inverse();

					// check abs_residual and yieldCondition for active yield surfaces
					bool convergenceFlagYieldCondition(true);
					if (absResidual > toleranceResidual)
						convergenceFlagYieldCondition = false;
					else
					{
						for (int count=0; count<numYieldSurfaces; count++)
						{
							if (yieldConditionFlag(count) == INACTIVE)
								continue;
							if (fabs(yieldCondition(count)) > toleranceYieldSurface*e_mod)
							{
								convergenceFlagYieldCondition = false;
								break;
							}
						}
					}

					if (convergenceFlagYieldCondition==true)
					{
						// convergence is achieved - now check of the deltaGamma is nonnegative and all other yield surfaces are valid
						for (int count=0; count<numYieldSurfaces; count++)
						{
							if (yieldConditionFlag(count) == INACTIVE)
							{
								if (yieldCondition(count) > toleranceYieldSurface*e_mod)
								{
									convergenceFlagYieldCondition = false;
									iteration=maxSteps;
								    break;
								}
							}
							else
							{
								if (deltaGamma(count)<0)
							    {
									convergenceFlagYieldCondition = false;
									iteration=maxSteps;
								    break;
							    }
							}
						}

						if (convergenceFlagYieldCondition)
						{
							convergedInternal = true;
#ifdef ENABLE_DEBUG
							std::cout << "convergence after " << iteration << " iterations" << std::endl << std::endl;
#endif
						}
					}
					if (convergedInternal && cutbackFactorExternal==1)
					{
						// compute elasto plastic d_matrix
						// compute G_matrix
						int curYieldFunction = 0;
						matG.setZero(numActiveYieldFunctions,numActiveYieldFunctions);
						vectorN.resize(numYieldSurfaces);
						for (int count=0; count<numYieldSurfaces; count++)
						{
							if (yieldConditionFlag(count)==INACTIVE)
								continue;
							int curYieldFunction2 = 0;
							for (int count2=0; count2<=count; count2++)
							{
								if (yieldConditionFlag(count2)==INACTIVE)
									continue;

								matG(curYieldFunction,curYieldFunction2) = (dF_dsigma[count].transpose() * hessian * dF_dsigma[count2]).lazy()(0);
								// copy symmetric part
								if (count!=count2)
									matG(curYieldFunction2,curYieldFunction) = matG(curYieldFunction,curYieldFunction2);

								curYieldFunction2++;
							}

							// N
							vectorN[count] = hessian * dF_dsigma[count];
							curYieldFunction++;
						}

						// solve linearized system of equations for G_inv
						assert(fabs(matG.determinant())>toleranceDeterminant);
						matG.computeInverse(&matGInv);

						// compute elasto_plastic matrix
						tmpStiffness = hessian;
						curYieldFunction = 0;
						for (int count=0; count<numYieldSurfaces; count++)
						{
							if (yieldConditionFlag(count)==INACTIVE)
								continue;
							int curYieldFunction2 = 0;
							for (int count2=0; count2<numYieldSurfaces; count2++)
							{
								if (yieldConditionFlag(count2)==INACTIVE)
									continue;
								tmpStiffness-=matGInv(curYieldFunction,curYieldFunction2)*vectorN[count]*vectorN[count2].transpose();
								curYieldFunction2++;
							}
							curYieldFunction++;
						}

						//update new history variables
						rdEpsilonPdEpsilon = dElastInv * (dElast - tmpStiffness);

#ifdef ENABLE_DEBUG
						std::cout << "numberOfExternalSteps (totalCurstrainincreases) " << numberOfExternalCutbacks <<std::endl;
						std::cout << "numberOfInternalIterations (Newton iterations to find delta2Gamma and deltaStress) " << numberOfInternalIterations <<std::endl;
#endif
						return;

	/*
						// add damage part

						if (mDamage)
						{
							double factor  = (curPlastStrain(0)-curPlastStrain(1))/2.;
							double help_scalar = sqrt(factor*factor+curPlastStrain(2)*curPlastStrain(2)*0.25);
							double epsilon_p_max_plane;
							bool mixed_p_max;
							if ((curPlastStrain(0)+curPlastStrain(1))*0.5-help_scalar>0)
							{
								// mixed interpolation
								epsilon_p_max_plane = 0.5 * sqrt(4*(curPlastStrain(0)*curPlastStrain(0)+curPlastStrain(1)*curPlastStrain(1))+2*curPlastStrain(2)*curPlastStrain(2));
								mixed_p_max = true;
							}
							else
							{
								epsilon_p_max_plane = (curPlastStrain(0)+curPlastStrain(1))/2.+sqrt(factor*factor+curPlastStrain(2)*curPlastStrain(2)*0.25);
								mixed_p_max = false;
							}

							double l_eq,x_s,y_s,tan_beta;
							if (fabs(epsilon_p_max_plane)>1e-10)
							{
								tan_beta = curPlastStrain(3)/epsilon_p_max_plane;
								x_s = sqrt(1./(1./(l_eq_plane*l_eq_plane)+tan_beta*tan_beta/(l_eq_circ*l_eq_circ)));
								y_s = x_s * tan_beta;
								l_eq = sqrt(x_s*x_s+y_s*y_s);
							}
							else
							{
								l_eq = l_eq_circ;
							}

							double kappa_d = kappa_d_unscaled/l_eq;
							rNewStaticData->mOmega = 1.-exp(-rNewStaticData->mKappa/kappa_d);

							//limit damage variable to a maximal value in order to prevent zero eigenvalues
							if (rNewStaticData->mOmega>=MAX_OMEGA)
							{
								rNewStaticData->mOmega = MAX_OMEGA;
								if (rStiffness!=0)
								{
									// secante
									*rStiffness *= 1-MAX_OMEGA;
								}
								if (rStress!=0)
									*rStress *= 1-MAX_OMEGA;
								return;
							}


							if (rStiffness!=0)
							{

								if (rNewStaticData->mOmega>=(rOldStaticData->mOmega))
								{
									Eigen::Matrix<double,4,4> d_epsilon_pl_d_epsilon;
									d_epsilon_pl_d_epsilon = dElastInv * (dElast - tmpStiffness);

									double d_omega_d_kappa_eq = 1./kappa_d*exp(-rNewStaticData->mKappa/kappa_d);

									Eigen::Matrix<double,4,1> d_kappa_eq_depsilon_p;
									factor = 1./rNewStaticData->mKappa;
									d_kappa_eq_depsilon_p(0) = factor * curPlastStrain(0);
									d_kappa_eq_depsilon_p(1) = factor * curPlastStrain(1);
									d_kappa_eq_depsilon_p(2) = factor * 0.5*curPlastStrain(2); // w.r.t. gamma
									d_kappa_eq_depsilon_p[3] = factor * curPlastStrain(3);

									//memcpy(help_matrix,result,16*sizeof(double));

									// derivative w.r.t. kappa_d

									Eigen::Matrix<double,4,1> d_omega_d_epsilon;
									if (fabs(epsilon_p_max_plane)>1e-10) // same tolerance as above
									{
										double d_kappad_d_leq = -kappa_d_unscaled/(l_eq*l_eq);
										double d_leq_d_xs = x_s/l_eq;
										double d_leq_d_ys = y_s/l_eq;
										double d_xs_d_tanbeta = -tan_beta*x_s*x_s*x_s/(l_eq_circ*l_eq_circ);
										double d_tan_beta_d_epsilon_p_zz = 1./epsilon_p_max_plane;
										double d_tan_beta_d_epsilon_p_max_plane = -curPlastStrain[3]/(epsilon_p_max_plane*epsilon_p_max_plane);
										double d_epsilon_p_max_plane_d_epsilon_p_xx, d_epsilon_p_max_plane_d_epsilon_p_yy, d_epsilon_p_max_plane_d_epsilon_p_xy;
										if (!mixed_p_max)
										{
											factor=1./(2.*sqrt((curPlastStrain(0)-curPlastStrain(1))*(curPlastStrain(0)-curPlastStrain(1))+curPlastStrain(2)*curPlastStrain(2)));
											d_epsilon_p_max_plane_d_epsilon_p_xx = 0.5+(curPlastStrain(0)-curPlastStrain(1))*factor;
											d_epsilon_p_max_plane_d_epsilon_p_yy = 0.5+(curPlastStrain(1)-curPlastStrain(0))*factor;
											d_epsilon_p_max_plane_d_epsilon_p_xy = curPlastStrain(2)*factor;
										}
										else
										{
											d_epsilon_p_max_plane_d_epsilon_p_xx = curPlastStrain(0)/epsilon_p_max_plane;
											d_epsilon_p_max_plane_d_epsilon_p_yy = curPlastStrain(1)/epsilon_p_max_plane;
											d_epsilon_p_max_plane_d_epsilon_p_xy = 0.5*curPlastStrain(2)/epsilon_p_max_plane;
										}

										Eigen::Matrix<double,4,1> d_epsilon_p_max_plane_d_epsilon;
										d_epsilon_p_max_plane_d_epsilon (0) =
											d_epsilon_p_max_plane_d_epsilon_p_xx *d_epsilon_pl_d_epsilon(0)+
											d_epsilon_p_max_plane_d_epsilon_p_yy *d_epsilon_pl_d_epsilon(1)+
											d_epsilon_p_max_plane_d_epsilon_p_xy *d_epsilon_pl_d_epsilon(2);

										d_epsilon_p_max_plane_d_epsilon (1) =
											d_epsilon_p_max_plane_d_epsilon_p_xx *d_epsilon_pl_d_epsilon(4)+
											d_epsilon_p_max_plane_d_epsilon_p_yy *d_epsilon_pl_d_epsilon(5)+
											d_epsilon_p_max_plane_d_epsilon_p_xy *d_epsilon_pl_d_epsilon(6);

										d_epsilon_p_max_plane_d_epsilon (2) =
											d_epsilon_p_max_plane_d_epsilon_p_xx *d_epsilon_pl_d_epsilon(8)+
											d_epsilon_p_max_plane_d_epsilon_p_yy *d_epsilon_pl_d_epsilon(9)+
											d_epsilon_p_max_plane_d_epsilon_p_xy *d_epsilon_pl_d_epsilon(10);

										d_epsilon_p_max_plane_d_epsilon (3) =
											d_epsilon_p_max_plane_d_epsilon_p_xx *d_epsilon_pl_d_epsilon(12)+
											d_epsilon_p_max_plane_d_epsilon_p_yy *d_epsilon_pl_d_epsilon(13)+
											d_epsilon_p_max_plane_d_epsilon_p_xy *d_epsilon_pl_d_epsilon(14);

										Eigen::Matrix<double,4,1> d_tanbeta_d_epsilon;
										d_tanbeta_d_epsilon (0) = d_tan_beta_d_epsilon_p_zz*d_epsilon_pl_d_epsilon(3)+
																  d_tan_beta_d_epsilon_p_max_plane *d_epsilon_p_max_plane_d_epsilon (0);
										d_tanbeta_d_epsilon (1) = d_tan_beta_d_epsilon_p_zz*d_epsilon_pl_d_epsilon(7)+
																  d_tan_beta_d_epsilon_p_max_plane *d_epsilon_p_max_plane_d_epsilon (1);
										d_tanbeta_d_epsilon (2) = d_tan_beta_d_epsilon_p_zz*d_epsilon_pl_d_epsilon(11)+
																  d_tan_beta_d_epsilon_p_max_plane *d_epsilon_p_max_plane_d_epsilon (2);
										d_tanbeta_d_epsilon (3) = d_tan_beta_d_epsilon_p_zz*d_epsilon_pl_d_epsilon(15)+
																  d_tan_beta_d_epsilon_p_max_plane *d_epsilon_p_max_plane_d_epsilon [3];


										Eigen::Matrix<double,4,1> d_xs_d_epsilon;
										d_xs_d_epsilon (0) = d_xs_d_tanbeta * d_tanbeta_d_epsilon(0);
										d_xs_d_epsilon (1) = d_xs_d_tanbeta * d_tanbeta_d_epsilon(1);
										d_xs_d_epsilon (2) = d_xs_d_tanbeta * d_tanbeta_d_epsilon(2);
										d_xs_d_epsilon (3) = d_xs_d_tanbeta * d_tanbeta_d_epsilon(3);

										Eigen::Matrix<double,4,1> d_ys_d_epsilon;
										d_ys_d_epsilon (0) = d_xs_d_epsilon (0)*tan_beta + x_s * d_tanbeta_d_epsilon(0);
										d_ys_d_epsilon (1) = d_xs_d_epsilon (1)*tan_beta + x_s * d_tanbeta_d_epsilon(1);
										d_ys_d_epsilon (2) = d_xs_d_epsilon (2)*tan_beta + x_s * d_tanbeta_d_epsilon(2);
										d_ys_d_epsilon (3) = d_xs_d_epsilon (3)*tan_beta + x_s * d_tanbeta_d_epsilon(3);

										Eigen::Matrix<double,4,1> d_leq_d_epsilon;
										d_leq_d_epsilon (0) = d_leq_d_xs*d_xs_d_epsilon (0) + d_leq_d_ys*d_ys_d_epsilon (0);
										d_leq_d_epsilon (1) = d_leq_d_xs*d_xs_d_epsilon (1) + d_leq_d_ys*d_ys_d_epsilon (1);
										d_leq_d_epsilon (2) = d_leq_d_xs*d_xs_d_epsilon (2) + d_leq_d_ys*d_ys_d_epsilon (2);
										d_leq_d_epsilon (3) = d_leq_d_xs*d_xs_d_epsilon (3) + d_leq_d_ys*d_ys_d_epsilon (3);


										Eigen::Matrix<double,4,1> d_kappad_d_epsilon;
										d_kappad_d_epsilon (0) = d_kappad_d_leq * d_leq_d_epsilon (0);
										d_kappad_d_epsilon (1) = d_kappad_d_leq * d_leq_d_epsilon (1);
										d_kappad_d_epsilon (2) = d_kappad_d_leq * d_leq_d_epsilon (2);
										d_kappad_d_epsilon (3) = d_kappad_d_leq * d_leq_d_epsilon (3);

										factor = rNewStaticData->mKappa/(kappa_d*kappa_d)*exp(-rNewStaticData->mKappa/kappa_d);
										d_omega_d_epsilon (0)=- factor * d_kappad_d_epsilon (0);
										d_omega_d_epsilon (1)=- factor * d_kappad_d_epsilon (1);
										d_omega_d_epsilon (2)=- factor * d_kappad_d_epsilon (2);  // w.r.t. gamma
										d_omega_d_epsilon (3)=- factor * d_kappad_d_epsilon (3);
									}
									else
									{
										d_omega_d_epsilon.setZero(4);
									}

									// derivative w.r.t. epsilon_xx
									d_omega_d_epsilon (0) += d_omega_d_kappa_eq*(
																 d_kappa_eq_depsilon_p(0)*d_epsilon_pl_d_epsilon(0)+
																 d_kappa_eq_depsilon_p(1)*d_epsilon_pl_d_epsilon(1)+
																 d_kappa_eq_depsilon_p(2)*d_epsilon_pl_d_epsilon(2)+
																 d_kappa_eq_depsilon_p(3)*d_epsilon_pl_d_epsilon(3));

									// derivative w.r.t. epsilon_yy
									d_omega_d_epsilon (1) += d_omega_d_kappa_eq*(
																 d_kappa_eq_depsilon_p(0)*d_epsilon_pl_d_epsilon[4]+
																 d_kappa_eq_depsilon_p(1)*d_epsilon_pl_d_epsilon[5]+
																 d_kappa_eq_depsilon_p(2)*d_epsilon_pl_d_epsilon[6]+
																 d_kappa_eq_depsilon_p(3)*d_epsilon_pl_d_epsilon[7]);

									// derivative w.r.t. gamma_xy
									d_omega_d_epsilon (2) += d_omega_d_kappa_eq*(
																 d_kappa_eq_depsilon_p(0)*d_epsilon_pl_d_epsilon(8)+
																 d_kappa_eq_depsilon_p(1)*d_epsilon_pl_d_epsilon(9)+
																 d_kappa_eq_depsilon_p(2)*d_epsilon_pl_d_epsilon(10)+
																 d_kappa_eq_depsilon_p[3]*d_epsilon_pl_d_epsilon(11));

									// derivative w.r.t. sigma_zz
									d_omega_d_epsilon [3] += d_omega_d_kappa_eq*(
																 d_kappa_eq_depsilon_p(0)*d_epsilon_pl_d_epsilon(12)+
																 d_kappa_eq_depsilon_p(1)*d_epsilon_pl_d_epsilon(13)+
																 d_kappa_eq_depsilon_p(2)*d_epsilon_pl_d_epsilon(14)+
																 d_kappa_eq_depsilon_p(3)*d_epsilon_pl_d_epsilon(15));

									*rStiffness *= 1.-rNewStaticData->mOmega;
									*rStiffness -= trialStress*d_omega_d_epsilon.transpose().lazy();

								}
								else
								{
									if (rStiffness!=0)
										*rStiffness *= 1.-rNewStaticData->mOmega;
								}
							}
							if (rStress!=0)
							{
								*rStress *= 1.-rNewStaticData->mOmega;
							}
						}
						break; // end of damage_flag
						*/
					}
					else
					{
						if (iteration<maxSteps)
						{
							//compute delta gamma
							int curYieldFunction = 0;
							matG.setZero(numActiveYieldFunctions,numActiveYieldFunctions);
							for (int count=0; count<numYieldSurfaces; count++)
							{
								if (yieldConditionFlag(count)==INACTIVE)
									continue;
								int curYieldFunction2 = 0;
								for (int count2=0; count2<=count; count2++)
								{
									if (yieldConditionFlag(count2)==INACTIVE)
										continue;

									matG(curYieldFunction,curYieldFunction2) = (dF_dsigma[count].transpose()*hessian*dF_dsigma[count2])(0);
									// copy symmetric part
									if (count!=count2)
									{
										matG(curYieldFunction2,curYieldFunction) = matG(curYieldFunction,curYieldFunction2);
									}

									curYieldFunction2++;
								}

								curYieldFunction++;
							}
							assert(curYieldFunction==numActiveYieldFunctions);

							// solve linearized system of equations for G_inv
							assert(fabs(matG.determinant())>toleranceDeterminant);

							matG.computeInverse(&matGInv);

							// compute deltaGamma
							Eigen::Matrix<double,4,1> helpVector = hessian * residual;
							Eigen::Matrix<double,4,1> helpVector2;
							for (int count=0; count<numYieldSurfaces; count++)
							{
								if (yieldConditionFlag(count)==INACTIVE)
									continue;

								helpVector2(count) = yieldCondition(count) - dF_dsigma[count].dot(helpVector);
							}

							curYieldFunction = 0;
							for (int count=0; count<numYieldSurfaces; count++)
							{
								if (yieldConditionFlag(count)==INACTIVE)
									continue;
								delta2Gamma(count) = 0.;
								int curYieldFunction2 = 0;
								for (int count2=0; count2<numYieldSurfaces; count2++)
								{
									if (yieldConditionFlag(count2)==INACTIVE)
										continue;

									delta2Gamma(count) += matGInv(curYieldFunction,curYieldFunction2)*helpVector2(count2);

									curYieldFunction2++;
								}
								curYieldFunction++;
							}

#ifdef ENABLE_DEBUG
							std::cout << "delta2Gamma " << delta2Gamma.transpose() << std::endl<< std::endl;
#endif
							//******************************************************************
							// compute increments for stress
							//******************************************************************
							helpVector = residual;

							for (int count=0; count<numYieldSurfaces; count++)
							{
								if (yieldConditionFlag(count)==INACTIVE)
									continue;
#ifdef ENABLE_DEBUG
                                std::cout<<"dfdsigma " << dF_dsigma[count].transpose()<<std::endl<< std::endl;
#endif
								helpVector += delta2Gamma(count)*dF_dsigma[count];
							}

							deltaStress =  hessian * helpVector;
							//internal line search convergedExternal
							double deltaCutbackFactor(1.);
							double cutbackFactorLS(deltaCutbackFactor);
							bool convergedLS(false);

							//norm of the residual in the local Newton iteration, used as convergence indicator
							// in the local iteration, the unknowns are deltaSigma and delta2Gamma
							// whereas the residuals are the difference between (total -elastic) and (plastic strain)
							// and the yield conditions
		                    double normInit = residual.squaredNorm();
		                    for (int count=0; count<numYieldSurfaces; count++)
		                    {
		                        if (yieldConditionFlag(count)==INACTIVE)
		                            continue;
		                        normInit +=yieldCondition(count)*yieldCondition(count);
	                        }
		                    int numberOfLinesearchSteps(0);
		                    while (!convergedLS)
							{
		                    	numberOfLinesearchSteps++;
		                    	convergedLS = true;

								deltaGammaLS= deltaGamma + cutbackFactorLS * delta2Gamma;
		                        deltaPlasticStrainLS = cutbackFactorLS * (dElastInv * deltaStress);
								epsilonPLS = rEpsilonP + deltaPlasticStrainLS;
								stressLS = trialStress - cutbackFactorLS*deltaStress;

#ifdef ENABLE_DEBUG
		                    	std::cout << "delta2Gamma " << delta2Gamma.transpose() << std::endl<< std::endl;
		                        std::cout << "deltaPlasticStrainLS " << deltaPlasticStrainLS.transpose() << std::endl<< std::endl;
		                        std::cout << "epsilonPLS " << epsilonPLS.transpose() << std::endl<< std::endl;
		                		std::cout << "stressLS " << stressLS.transpose() << std::endl<< std::endl;
#endif


								//* calculate yield condition and yield condition flag (active or not)
								//* Drucker Prager
								if (yieldConditionFlag(0)==ACTIVE)
								{
								    yieldConditionLS(0) = YieldSurfaceDruckerPrager2D(stressLS, BETA, H_P);
									if (!YieldSurfaceDruckerPrager2DDerivatives(dF_dsigma[0],0,stressLS,BETA))
									{
										//no convergence, decrease line search step
										convergedLS =false;
									}
#ifdef ENABLE_DEBUG
									std::cout << "dF_dsigma[0] " <<  std::endl << dF_dsigma[0].transpose() << std::endl << std::endl;
#endif
								}

								// Rounded Rankine
								if (yieldConditionFlag(1)==ACTIVE)
								{
								    yieldConditionLS(1) = YieldSurfaceRankine2DRounded(stressLS, f_ct);
									YieldSurfaceRankine2DRoundedDerivatives(dF_dsigma[1],0,stressLS);
#ifdef ENABLE_DEBUG
									std::cout << "dF_dsigma[1] " <<  std::endl << dF_dsigma[1].transpose() << std::endl << std::endl;
#endif
								}

								// residual in line search
								residualLS = lastPlastStrain-epsilonPLS;

								for (int count=0; count<numYieldSurfaces; count++)
								{
									if (yieldConditionFlag[count] == ACTIVE)
									{
										residualLS += deltaGammaLS(count)*dF_dsigma[count];
									}
								}

#ifdef ENABLE_DEBUG
			                    std::cout << "residual linesearch" <<  std::endl << residualLS.transpose() << std::endl << std::endl;
#endif
			                    double normCurr = residualLS.squaredNorm();
			                    for (int count=0; count<numYieldSurfaces; count++)
			                    {
			                        if (yieldConditionFlag(count)==INACTIVE)
			                            continue;
			                        normCurr +=yieldConditionLS(count)*yieldConditionLS(count);
		                        }
#ifdef ENABLE_DEBUG
			                    std::cout << "normInit " << normInit << "normCurr " << normCurr<<std::endl<< std::endl;
#endif

			                    // check of relative norm of residual is sufficiently decreasing or if it is almost zero
			                    if ((normCurr/normInit)*(normCurr/normInit)>1-0.5*cutbackFactorLS && normCurr>toleranceResidual*toleranceResidual)
    	                        	convergedLS = false;
		                        if (cutbackFactorLS<=minCutbackFactorLS)
		                        	convergedLS = true;

		                        if (!convergedLS)
		                        {
		                        	deltaCutbackFactor*=0.5;
		                        	cutbackFactorLS-=deltaCutbackFactor;
		                        }
							}
		                    trialStress = stressLS;
		                    deltaGamma = deltaGammaLS;
							rEpsilonP = epsilonPLS;

							//update equivalente plastic strain
							rEqPlasticStrain += sqrt(deltaPlasticStrainLS[0]*deltaPlasticStrainLS(0)+deltaPlasticStrainLS(1)*deltaPlasticStrainLS(1)+
											0.5 *(deltaPlasticStrainLS(2)*deltaPlasticStrainLS(2))+deltaPlasticStrainLS(3)*deltaPlasticStrainLS(3));
#ifdef ENABLE_DEBUG
							double tmp(sqrt(rEpsilonP(0)*rEpsilonP(0) + rEpsilonP(1)*rEpsilonP(1)+ 0.5*rEpsilonP(2)*rEpsilonP(2)+ rEpsilonP(3)*rEpsilonP(3)));
	                        std::cout << std::endl << "rEqPlasticStrain " << rEqPlasticStrain << "Norm of epsilonp " << tmp << "delta between " << rEqPlasticStrain - tmp << std::endl<< std::endl;
	                        if (yieldConditionFlag(0)==ACTIVE)
	                            std::cout << "dF_dsigma[0] at end of line search " <<  std::endl << dF_dsigma[0].transpose() << std::endl << std::endl;
	                        if (yieldConditionFlag(1)==ACTIVE)
	                            std::cout << "dF_dsigma[1] at end of line search " <<  std::endl << dF_dsigma[1].transpose() << std::endl << std::endl;
	                        std::cout << "numberOfLinesearchSteps " << numberOfLinesearchSteps << std::endl;
#endif
						}
					}
				} // end of loop
            }
        }
        catch (NuTo::MechanicsException& e)
        {
            e.AddMessage("[NuTo::NonlocalDamagePlasticity::ReturnMapping2D] Error performing return mapping procedure.");
            throw e;
        }
        catch (...)
        {
               throw MechanicsException("[NuTo::NonlocalDamagePlasticity::ReturnMapping2D] Error performing return mapping procedure.");
        }

        if (convergedInternal)
        {
#ifdef ENABLE_DEBUG
			std::cout << "numberOfInternalIterations " << numberOfInternalIterations -  prevNumberOfInternalIterations<< "(" << numberOfInternalIterations << ")" << std::endl;
            std::cout << "convergence for external cutback factor" << std::endl;
#endif
			prevNumberOfInternalIterations = numberOfInternalIterations;

        	if (cutbackFactorExternal==1.)
                convergedExternal=true;
            else
            {
				lastPlastStrain = rEpsilonP;
				lastEqPlasticStrain = rEqPlasticStrain;
                cutbackFactorExternal +=deltaCutbackFactorExternal;
                if (cutbackFactorExternal>1)
                {
                    deltaCutbackFactorExternal -= cutbackFactorExternal -1.;
                    cutbackFactorExternal = 1.;
                }
            }
        }
        else
        {
#ifdef ENABLE_DEBUG
            std::cout << "decrease external cutback factor" << std::endl;
#endif
        	// no convergence in return mapping
        	deltaCutbackFactorExternal*=0.5;
            cutbackFactorExternal-=deltaCutbackFactorExternal;
        }

    }
    if (cutbackFactorExternal<=minCutbackFactor)
    {
        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::ReturnMapping2D] No convergence can be obtained in the return mapping procedure.");
    }
}

//! @brief calculates the first and second derivative of the second Rankine yield surface with respect to the stress
//! @param rStress current stress
//! @param rSigma_1 first principal stress
//! @param rSigma_2 second principal stress
//! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm 0.5*value_sqrt$
//! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm 0.5*value_sqrt$
//! @return yield condition
double NuTo::NonlocalDamagePlasticity::YieldSurfaceRankine2DRounded(Eigen::Matrix<double,4,1>& rStress, double rFct)const
{
	double value_sum  = 0.5*(rStress(0)+rStress(1));
	double help_scalar = (rStress(0)-rStress(1));
	double value_sqrt = sqrt(help_scalar*help_scalar+4.*rStress(2)*rStress(2));
	double sigma_1(value_sum+0.5*value_sqrt);
	double sigma_2(value_sum-0.5*value_sqrt);
	//* (rounded) Rankine
#ifdef ENABLE_DEBUG
	std::cout << "p1 " << sigma_1 << " p2 " << sigma_2 << " p3 " << rStress[3] << std::endl << std::endl;
#endif
	if (rStress[3]<0)
	{
		//sigma_3 is negative
		if (sigma_1<0)
		{
			// sigma_1 is negative and as a consequence sigma_2 is also negative
			if (rStress[3]>sigma_1)
			{
#ifdef ENABLE_DEBUG
				std::cout << std::endl << " all negative f1" << std::endl;
#endif
				return rStress[3]-rFct;
			}
			else
			{
#ifdef ENABLE_DEBUG
				std::cout << std::endl << " all negative f3" << std::endl;
#endif
				return sigma_1-rFct;
			}
		}
		else
		{
			//sigma_1 is positive
			if (sigma_2<0)
			{
				//sigma_2 is negative
#ifdef ENABLE_DEBUG
				std::cout << std::endl << " f1" << std::endl;
#endif
				return sigma_1 - rFct;
			}
			else
			{
				//sigma_2 is positive
#ifdef ENABLE_DEBUG
				std::cout << std::endl << " f1,f2" << std::endl;
#endif
				return sqrt(sigma_1*sigma_1 + sigma_2*sigma_2)-rFct;
			}
		}
	}
	else
	{
		//sigma_3 is positive
		if (sigma_1<0)
		{
			// sigma_1 is negative and as a consequence sigma_2 is also negative
#ifdef ENABLE_DEBUG
			std::cout << std::endl << " f3" << std::endl;
#endif
			return rStress[3]-rFct;
		}
		else
		{
		    //sigma_1 is positive
			if (sigma_2<0)
			{
				//sigma_2 is negative
#ifdef ENABLE_DEBUG
				std::cout << std::endl << " f1,f3" << std::endl;
#endif
				return sqrt(sigma_1*sigma_1 + rStress[3] * rStress[3]) - rFct;
			}
			else
			{
				//sigma_2 is positive
#ifdef ENABLE_DEBUG
				std::cout << std::endl << " f1,f2,f3" << std::endl;
#endif
				return sqrt(sigma_1*sigma_1 + sigma_2*sigma_2 + rStress[3] * rStress[3])-rFct;
			}
		}
	}
}

//! @brief calculates the first and second derivative of the Rankine yield surface with respect to the stress
//! @param rdF_dSigma return value (first derivative)
//! @param rd2F_d2Sigma return value (second derivative)
//! @param rStress current stress
//! @param value_sqrt second term for the calculation of the principal stresses in 2D $sigma_{1,2}= \dfrac{s1+s2}{2} \pm 0.5*value_sqrt$
void NuTo::NonlocalDamagePlasticity::YieldSurfaceRankine2DRoundedDerivatives(Eigen::Matrix<double,4,1>& rdF_dSigma,Eigen::Matrix<double,4,4>* rd2F_d2Sigma,
		Eigen::Matrix<double,4,1>& rStress)const
{
	double value_sum  = 0.5*(rStress(0)+rStress(1));
	double help_scalar = (rStress(0)-rStress(1));
	double value_sqrt = sqrt(help_scalar*help_scalar+4.*rStress(2)*rStress(2));
	double sigma_1(value_sum+0.5*value_sqrt);
	double sigma_2(value_sum-0.5*value_sqrt);

	double factor;
	if (rStress[3]<0)
	{
		//sigma_3 is negative
		if (sigma_1<0)
		{
			// sigma_1 is negative and as a consequence sigma_2 is also negative
			if (rStress[3]>sigma_1)
			{
		        rdF_dSigma(0) = 0.;
		        rdF_dSigma(1) = 0.;
		        rdF_dSigma(2) = 0.;
		        rdF_dSigma(3) = 1.;
		        //second derivatives are not required
			}
			else
			{
			    // f = f(sigma_1)
				if (value_sqrt<1e-12)
			        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::YieldSurfaceRankine2DRounded] value_sqrt<1e-12 should not happen, since sigma_1>0 and sigma_2<0");

			    rdF_dSigma(0) = (value_sqrt+rStress(0)-rStress(1))/(2.*value_sqrt);
			    rdF_dSigma(1) = (value_sqrt-rStress(0)+rStress(1))/(2.*value_sqrt);
			    rdF_dSigma(2) = 2*rStress(2)/value_sqrt;
			    rdF_dSigma(3) = 0.;
			    //second derivatives are not required
			}
		}
		else
		{
			//sigma_1 is positive
			if (sigma_2<0)
			{
			    // f = f(sigma_1)
				if (value_sqrt<1e-12)
			        throw MechanicsException("[NuTo::NonlocalDamagePlasticity::YieldSurfaceRankine2DRounded] value_sqrt<1e-12 should not happen, since sigma_1>0 and sigma_2<0");

			    rdF_dSigma(0) = (value_sqrt+rStress(0)-rStress(1))/(2.*value_sqrt);
			    rdF_dSigma(1) = (value_sqrt-rStress(0)+rStress(1))/(2.*value_sqrt);
			    rdF_dSigma(2) = 2*rStress(2)/value_sqrt;
			    rdF_dSigma(3) = 0.;

			    // store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
			    if (rd2F_d2Sigma!=0)
			    {
					double factor = 1./(value_sqrt*value_sqrt*value_sqrt);
					(*rd2F_d2Sigma)(0,0) = factor*2.*rStress[2]*rStress[2];
					(*rd2F_d2Sigma)(0,1) = factor*(-2.)*rStress[2]*rStress[2];
					(*rd2F_d2Sigma)(0,2) = factor*(2.)*(rStress[1]-rStress[0])*rStress[2];
					(*rd2F_d2Sigma)(0,3) = 0;
					(*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
					(*rd2F_d2Sigma)(1,1) = factor*2.*rStress[2]*rStress[2];
					(*rd2F_d2Sigma)(1,2) = factor*(2.)*(rStress[0]-rStress[1])*rStress[2];
					(*rd2F_d2Sigma)(1,3) = 0;
					(*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
					(*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
					(*rd2F_d2Sigma)(2,2) = factor*(2.)*(rStress[0]-rStress[1])*(rStress[0]-rStress[1]);
					(*rd2F_d2Sigma)(2,3) = 0.;
					(*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
					(*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
					(*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
					(*rd2F_d2Sigma)(3,3) = 0.;
			    }
			}
			else
			{
				// f = f(sigma_1, sigma_2)
				//sigma_2 is positive
		        factor = 1./sqrt(rStress[0]*rStress[0]+rStress[1]*rStress[1]+
		                         2*rStress[2]*rStress[2]);
		        rdF_dSigma(0) = rStress[0]*factor;
		        rdF_dSigma(1) = rStress[1]*factor;
		        rdF_dSigma(2) = 2.*rStress[2]*factor;
		        rdF_dSigma(3) = 0.;

			    if (rd2F_d2Sigma!=0)
			    {
					// store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
					factor = sqrt(rStress[0]*rStress[0]+rStress[1]*rStress[1]+2*rStress[2]*rStress[2]);
					factor = 1./(factor*factor*factor);

					(*rd2F_d2Sigma)(0,0) = factor*(rStress[1]*rStress[1]+2.*rStress[2]*rStress[2]);
					(*rd2F_d2Sigma)(0,1) = -factor*rStress[0]*rStress[1];
					(*rd2F_d2Sigma)(0,2) = (-2.)*factor*(rStress[0]*rStress[2]);
					(*rd2F_d2Sigma)(0,3) = 0;
					(*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
					(*rd2F_d2Sigma)(1,1) = factor*(rStress[0]*rStress[0]+2.*rStress[2]*rStress[2]);
					(*rd2F_d2Sigma)(1,2) = (-2.)*factor*(rStress[1]*rStress[2]);
					(*rd2F_d2Sigma)(1,3) = 0;
					(*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
					(*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
					(*rd2F_d2Sigma)(2,2) = (2.)*factor*(rStress[0]*rStress[0]+rStress[1]*rStress[1]);
					(*rd2F_d2Sigma)(2,3) = 0.;
					(*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
					(*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
					(*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
					(*rd2F_d2Sigma)(3,3) = 0.;
			    }
			}
		}
	}
	else
	{
		//sigma_3 is positive
		if (sigma_1<0)
		{
			// sigma_1 is negative and as a consequence sigma_2 is also negative
			// f = f(sigma_3)
	        rdF_dSigma(0) = 0.;
	        rdF_dSigma(1) = 0.;
	        rdF_dSigma(2) = 0.;
	        rdF_dSigma(3) = 1.;

		    if (rd2F_d2Sigma!=0)
		    {
		    	(*rd2F_d2Sigma).setZero(4,4);
		    }
		}
		else
		{
			//sigma_1 is positive
			if (sigma_2<0)
			{
				//sigma_2 is negative
				//f = f( sigma_1,sigma_3)

				double yield = sqrt(sigma_1*sigma_1 + rStress(3)*rStress(3));
				factor = 1./yield;
				rdF_dSigma(0) = factor*sigma_1*(0.5+0.5*(rStress(0)-rStress(1))/value_sqrt);
		        rdF_dSigma(1) = factor*sigma_1*(0.5+0.5*(rStress(1)-rStress(0))/value_sqrt);
		        rdF_dSigma(2) = factor*sigma_1*(2.*rStress(2)/value_sqrt);
		        rdF_dSigma(3) = factor*rStress(3);

			    if (rd2F_d2Sigma!=0)
			    {
					factor *= factor;
					double dsigma1_dx = (value_sqrt+rStress(0)-rStress(1))/(2.*value_sqrt);
					double dsigma1_dy = (value_sqrt-rStress(0)+rStress(1))/(2.*value_sqrt);
					double dsigma1_dxy = 2.*rStress(2)/value_sqrt;

					double factor2 = 1./(value_sqrt*value_sqrt*value_sqrt);
					double dsigma1_dx2    = 2.*rStress(2)*rStress(2)*factor2;
					double dsigma1_dxdy   = -dsigma1_dx2;
					double dsigma1_dxdxy  = 2.*(rStress(1)-rStress(0))*rStress(2)*factor2;
					double dsigma1_dy2    = dsigma1_dx2;
					double dsigma1_dydxy  = 2.*(rStress(0)-rStress(1))*rStress(2)*factor2;
					double dsigma1_dxy2   = 2.*(rStress(0)-rStress(1))*(rStress(0)-rStress(1))*factor2;

					// store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
					(*rd2F_d2Sigma)(0,0) = factor*((dsigma1_dx*dsigma1_dx +sigma_1*dsigma1_dx2  )*yield-sigma_1*dsigma1_dx*rdF_dSigma(0));
					(*rd2F_d2Sigma)(0,1) = factor*((dsigma1_dx*dsigma1_dy +sigma_1*dsigma1_dxdy )*yield-sigma_1*dsigma1_dx*rdF_dSigma(1));
					(*rd2F_d2Sigma)(0,2) = factor*((dsigma1_dx*dsigma1_dxy+sigma_1*dsigma1_dxdxy)*yield-sigma_1*dsigma1_dx*rdF_dSigma(2));
					(*rd2F_d2Sigma)(0,3) = factor*rStress(3)*rdF_dSigma(0);
					(*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
					(*rd2F_d2Sigma)(1,1) = factor*((dsigma1_dy*dsigma1_dy +sigma_1*dsigma1_dy2  )*yield-sigma_1*dsigma1_dy*rdF_dSigma(1));
					(*rd2F_d2Sigma)(1,2) = factor*((dsigma1_dy*dsigma1_dxy+sigma_1*dsigma1_dydxy)*yield-sigma_1*dsigma1_dy*rdF_dSigma(2));
					(*rd2F_d2Sigma)(1,3) = factor*rStress(3)*rdF_dSigma(1);
					(*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
					(*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
					(*rd2F_d2Sigma)(2,2) = factor*((dsigma1_dxy*dsigma1_dxy +sigma_1*dsigma1_dxy2  )*yield-sigma_1*dsigma1_dxy*rdF_dSigma(2));
					(*rd2F_d2Sigma)(2,3) = factor*rStress(3)*rdF_dSigma(2);
					(*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
					(*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
					(*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
					(*rd2F_d2Sigma)(3,3) = factor*(rStress(3)*rdF_dSigma(3)-yield);
			    }
			}
			else
			{
				//sigma_2 is positive
				//f = f( sigma_1,sigma_2, sigma_3)
				factor = 1./sqrt(rStress(0)*rStress(0) + rStress(1)*rStress(1) + 2.*rStress(2)*rStress(2) + rStress(3)*rStress(3));
				rdF_dSigma(0) = factor*rStress(0);
		        rdF_dSigma(1) = factor*rStress(1);
		        rdF_dSigma(2) = factor*2.*rStress(2);
		        rdF_dSigma(3) = factor*rStress(3);

			    if (rd2F_d2Sigma!=0)
			    {

			    	factor*=factor*factor;
			    	(*rd2F_d2Sigma)(0,0) = factor * (rStress(1)*rStress(1) + 2.*rStress(2)*rStress(2) + rStress(3)*rStress(3));
			    	(*rd2F_d2Sigma)(0,1) = -factor * rStress(1) * rStress(0);
			    	(*rd2F_d2Sigma)(0,2) = -2.*factor * rStress(2) * rStress(0);
			    	(*rd2F_d2Sigma)(0,3) = -factor * rStress(3) * rStress(0);
			    	(*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
			    	(*rd2F_d2Sigma)(1,1) = factor * (rStress(0)*rStress(0) + 2.*rStress(2)*rStress(2) + rStress(3)*rStress(3));
			    	(*rd2F_d2Sigma)(1,2) = -2.*factor * rStress(2) * rStress(1);
			    	(*rd2F_d2Sigma)(1,3) = -factor * rStress(3) * rStress(1);
			    	(*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
			    	(*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
			    	(*rd2F_d2Sigma)(2,2) = 2.*factor*(rStress(0)*rStress(0) + rStress(1)*rStress(1) + rStress(3)*rStress(3));
			    	(*rd2F_d2Sigma)(2,3) = -2.*factor * rStress(3) * rStress(2);
			    	(*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
			    	(*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
			    	(*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
			    	(*rd2F_d2Sigma)(3,3) = factor * (rStress(0)*rStress(0) + rStress(1)*rStress(1) + 2.*rStress(2)*rStress(2));
			    }
			}
		}
	}
}

//! @brief calculates the first and second derivative of the second Rankine yield surface with respect to the stress
//! @param rStress current stress
//! @param rBeta parameter of the Drucker-Prager yield surface
//! @param rHP parameter of the Drucker-Prager yield surface
//! @return yield condition
double NuTo::NonlocalDamagePlasticity::YieldSurfaceDruckerPrager2D(Eigen::Matrix<double,4,1>& rStress, double rBeta, double rHP)const
{
	double invariante_1 = rStress(0)+rStress(1)+rStress(3);

	//*******************************************************************
	//*    second invariante                                                 *
	//*******************************************************************
	double invariante_2 = ((rStress(0)-rStress(1))*(rStress(0)-rStress(1))+
							(rStress(1)-rStress(3))*(rStress(1)-rStress(3))+
							(rStress(0)-rStress(3))*(rStress(0)-rStress(3)))/6.+
							 rStress(2)*rStress(2);

	//*******************************************************************
	//*    F_DEV:    term of stress-deviator with/without kinematic hard.    *
	//*******************************************************************
	double F_DEV = sqrt( invariante_2 ) ;

	//*******************************************************************
	//*    check yield first time                                            *
	//*******************************************************************
	double F_BETA = invariante_1 * rBeta/3. ;
	double F_FLOW = -rHP;

	//* calculate yield condition and yield condition flag (active or not)
	//* Drucker Prager
	return F_BETA + F_DEV + F_FLOW;
}

//! @brief calculates the first and second derivative of the Rankine yield surface with respect to the stress
//! @param rdF_dSigma return value (first derivative)
//! @param rd2F_d2Sigma return value (second derivative)
//! @param rStress current stress
//! @param rBETA parameter of the Drucker Prager yield surface
//! @return false, if the stress is on the hydrostatic axis otherwise true
bool NuTo::NonlocalDamagePlasticity::YieldSurfaceDruckerPrager2DDerivatives(Eigen::Matrix<double,4,1>& rdF_dSigma,Eigen::Matrix<double,4,4>* rd2F_d2Sigma,
        Eigen::Matrix<double,4,1>& rStress, double rBETA)const
{
	//*******************************************************************
	//*    second invariante                                                 *
	//*******************************************************************
	double invariante_2 = ((rStress(0)-rStress(1))*(rStress(0)-rStress(1))+
							(rStress(1)-rStress(3))*(rStress(1)-rStress(3))+
							(rStress(0)-rStress(3))*(rStress(0)-rStress(3)))/6.+
							 rStress(2)*rStress(2);

	if (fabs(invariante_2)<1e-12)
        return false;
    double factor = 1./(sqrt(invariante_2)*6.);

    rdF_dSigma(0) = factor * (2.*rStress[0]-rStress[1]-rStress[3])+rBETA/3.;
    rdF_dSigma(1) = factor * (2.*rStress[1]-rStress[0]-rStress[3])+rBETA/3.;
    rdF_dSigma(2) = factor * (6.*rStress[2]); /* vector notation from second order tensor */
    rdF_dSigma(3) = factor * (2.*rStress[3]-rStress[0]-rStress[1])+rBETA/3.;

    if (rd2F_d2Sigma!=0)
    {

		//TODO store upper part for new eigen version 3.0 rd2F_d2Sigma.selfadjointView<Upper>()
		factor = 1./(invariante_2*sqrt(invariante_2)*12.);
		(*rd2F_d2Sigma)(0,0) = factor * (rStress[1]*rStress[1]-2.*rStress[1]*rStress[3]+rStress[3]*rStress[3]+4.*rStress[2]*rStress[2]);
		(*rd2F_d2Sigma)(0,1) = factor * (-rStress[0]*rStress[1]+rStress[0]*rStress[3]+rStress[1]*rStress[3]-rStress[3]*rStress[3]-2.*rStress[2]*rStress[2]);
		(*rd2F_d2Sigma)(0,2) = factor * (-2.*rStress[2]*(2.*rStress[0]-rStress[1]-rStress[3]));
		(*rd2F_d2Sigma)(0,3) = factor * (rStress[0]*rStress[1]-rStress[0]*rStress[3]-rStress[1]*rStress[1]+rStress[1]*rStress[3]-2.*rStress[2]*rStress[2]);
		(*rd2F_d2Sigma)(1,0) = (*rd2F_d2Sigma)(0,1);
		(*rd2F_d2Sigma)(1,1) = factor * (rStress[0]*rStress[0]-2.*rStress[0]*rStress[3]+rStress[3]*rStress[3]+4.*rStress[2]*rStress[2]);
		(*rd2F_d2Sigma)(1,2) = factor * (2.*rStress[2]*(rStress[0]-2.*rStress[1]+rStress[3]));
		(*rd2F_d2Sigma)(1,3) = factor * (-rStress[0]*rStress[0]+rStress[0]*rStress[1]+rStress[0]*rStress[3]-rStress[1]*rStress[3]-2.*rStress[2]*rStress[2]);
		(*rd2F_d2Sigma)(2,0) = (*rd2F_d2Sigma)(0,2);
		(*rd2F_d2Sigma)(2,1) = (*rd2F_d2Sigma)(1,2);
		(*rd2F_d2Sigma)(2,2) = factor * 4.*(rStress[0]*rStress[0]-rStress[0]*rStress[1]-rStress[0]*rStress[3]+
										rStress[1]*rStress[1]-rStress[1]*rStress[3]+rStress[3]*rStress[3]);
		(*rd2F_d2Sigma)(2,3) = factor * 2. *(rStress[0]+rStress[1]-2.*rStress[3])*rStress[2];
		(*rd2F_d2Sigma)(3,0) = (*rd2F_d2Sigma)(0,3);
		(*rd2F_d2Sigma)(3,1) = (*rd2F_d2Sigma)(1,3);
		(*rd2F_d2Sigma)(3,2) = (*rd2F_d2Sigma)(2,3);
		(*rd2F_d2Sigma)(3,3) = factor * (rStress[0]*rStress[0]-2*rStress[0]*rStress[1]+rStress[1]*rStress[1]+4.*rStress[2]*rStress[2]);
    }
    return true;
}


//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::NonlocalDamagePlasticity::Info(unsigned short rVerboseLevel) const
{
    this->ConstitutiveBase::Info(rVerboseLevel);
    std::cout << "    Young's modulus      : " << this->mE << std::endl;
    std::cout << "    Poisson's ratio      : " << this->mNu << std::endl;
    std::cout << "    nonlocal radius      : " << this->mNonlocalRadius << std::endl;
    std::cout << "    tensile strength     : " << this->mTensileStrength << std::endl;
    std::cout << "    compressive strength : " << this->mCompressiveStrength << std::endl;
    std::cout << "    biaxial compressive strength : " << this->mBiaxialCompressiveStrength << std::endl;
    std::cout << "    fracture energy      : " << this->mFractureEnergy << std::endl;
}

// check parameters
void NuTo::NonlocalDamagePlasticity::CheckParameters()const
{
    this->CheckBiaxialCompressiveStrength(this->mNu);
    this->CheckYoungsModulus(this->mE);
    this->CheckPoissonsRatio(this->mNu);
    this->CheckNonlocalRadius(this->mNu);
    this->CheckTensileStrength(this->mNu);
    this->CheckCompressiveStrength(this->mNu);
    this->CheckFractureEnergy(this->mNu);
}


