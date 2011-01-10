// $Id: Multiscale.cpp 342 2010-10-18 12:39:08Z arnold2 $
// Multiscale.cpp
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

#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal3x3.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataMultiscale2DPlaneStrain.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/Multiscale.h"
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
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/IpDataBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/structures/unstructured/StructureIp.h"

#define sqrt3 1.732050808
#define MAX_OMEGA 0.999
//#define ENABLE_DEBUG

NuTo::Multiscale::Multiscale() : ConstitutiveEngineeringStressStrain()
{
    SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::Multiscale::serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
       std::cout << "start serialize Multiscale" << std::endl;
#endif
/*       ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveEngineeringStressStrain)
           & BOOST_SERIALIZATION_NVP(mE)
           & BOOST_SERIALIZATION_NVP(mNu)
           & BOOST_SERIALIZATION_NVP(mNonlocalRadius)
           & BOOST_SERIALIZATION_NVP(mTensileStrength)
           & BOOST_SERIALIZATION_NVP(mCompressiveStrength)
           & BOOST_SERIALIZATION_NVP(mBiaxialCompressiveStrength)
           & BOOST_SERIALIZATION_NVP(mFractureEnergy)
           & BOOST_SERIALIZATION_NVP(mYieldSurface)
           & BOOST_SERIALIZATION_NVP(mDamage);
*/
#ifdef DEBUG_SERIALIZATION
       std::cout << "finish serialize Multiscale" << std::endl;
#endif
    }
    BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Multiscale)
#endif // ENABLE_SERIALIZATION

//************ constitutive routines    ***********
//**  defined in structures/StructureIpConstitutive.cpp *********
//*************************************************
//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::Multiscale::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringPlasticStrain] not implemented.");
}
//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::Multiscale::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringPlasticStrain] not implemented.");
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::Multiscale::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringPlasticStrain] not implemented.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStress1D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] not implemented for 1D.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] not implemented for 1D.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStress2D& rEngineeringStress) const
{

    //solve
    double tolerance(1e-6);
    // this is a somehow weird situation, since for a normal material law nothing should be changed
    // since the material law is a full structure whose bc change, this can either be implemented with a cast (as I did)
    // or by removing the const flag from all material routines (which I do not consider is good)
    std::stringstream previousState;
    bool updatePerformed;
    const ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();

    Solve(rElement, rIp, rDeformationGradient,tolerance,previousState,updatePerformed);
    const StructureIp *fineScaleStructure = staticData->GetFineScaleStructure();


    //calculate average stress
    NuTo::FullMatrix<double> averageStress;

    double mlX, mlY;
    mlX = fineScaleStructure->GetDimensionX();
    mlY = fineScaleStructure->GetDimensionY();
    fineScaleStructure->ElementTotalGetAverageStress(mlX*mlY,averageStress);
    rEngineeringStress.mEngineeringStress[0] = averageStress(0,0);
    rEngineeringStress.mEngineeringStress[1] = averageStress(1,0);
    rEngineeringStress.mEngineeringStress[2] = averageStress(3,0);

    if (updatePerformed)
    {
        //restore previous state
        //don't be astonished, a const pointer does only mean, you are not allowed to change the pointer
        // but the object it points to can be changed - weird, but that's the way it is
//TODO test, if I need the const_cast
        const_cast<StructureIp*>(fineScaleStructure)->Restore(previousState);
   }
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] not implemented.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] not implemented.");
}

//  Damage /////////////////////////////////////
//! @brief ... calculate isotropic damage from deformation gradient in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDamage ... damage variable
void NuTo::Multiscale::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient1D& rDeformationGradient, double& rDamage) const
{
    throw MechanicsException("[NuTo::Multiscale::GetDamage] not implemented.");
}

//  Damage /////////////////////////////////////
//! @brief ... calculate isotropic damage from deformation gradient in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDamage ... damage variable
void NuTo::Multiscale::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient2D& rDeformationGradient, double& rDamage) const
{
    throw MechanicsException("[NuTo::Multiscale::GetDamage] not implemented.");
}

//  Damage /////////////////////////////////////
//! @brief ... calculate isotropic damage from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDamage ... damage variable
void NuTo::Multiscale::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient3D& rDeformationGradient, double& rDamage) const
{
    throw MechanicsException("[NuTo::Multiscale::GetDamage] not implemented.");
}

//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient,
        ConstitutiveTangentBase* rTangent) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient,
        ConstitutiveTangentBase* rTangent) const
{
    //solve
    double tolerance(1e-6);
    // this is a somehow weird situation, since for a normal material law nothing should be changed
    // since the material law is a full structure whose bc change, this can either be implemented with a cast (as I did)
    // or by removing the const flag from all material routines (which I do not consider is good)
    std::stringstream previousState;
    bool updatePerformed;
    const ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();

    Solve(rElement, rIp, rDeformationGradient,tolerance,previousState,updatePerformed);

    // use Schur complement to calculate the stiffness
    /*
    const StructureIp *fineScaleStructure = staticData->GetFineScaleStructure();

    //calculate average stress
    NuTo::FullMatrix<double> averageStress;

    double mlX, mlY;
    mlX = fineScaleStructure->GetDimensionX();
    mlY = fineScaleStructure->GetDimensionY();
    fineScaleStructure->ElementTotalGetAverageStress(mlX*mlY,averageStress);
    rEngineeringStress.mEngineeringStress[0] = averageStress(0,0);
    rEngineeringStress.mEngineeringStress[1] = averageStress(1,0);
    rEngineeringStress.mEngineeringStress[2] = averageStress(3,0);

    if (updatePerformed)
    {
        //restore previous state
        //don't be astonished, a const pointer does only mean, you are not allowed to change the pointer
        // but the object it points to can be changed - weird, but that's the way it is
//TODO test, if I need the const_cast
        const_cast<StructureIp*>(fineScaleStructure)->Restore(previousState);
   }

*/
    throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] not implemented add Schur complement.");
}

//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient,
        ConstitutiveTangentBase* rTangent) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain1D] not implemented.");
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const
{
    if (rElement->GetSection()==0)
        throw MechanicsException("[NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Section required to distinguish between plane stress and plane strain and thickness information.");
    if (rElement->GetSection()->GetType()==NuTo::Section::PLANE_STRESS)
        throw MechanicsException("[NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Nonlocal damage plasticity model not implemented for plane stress.");
    else
        return new ConstitutiveStaticDataMultiscale2DPlaneStrain();
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain3D] not implemented.");
}

//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
    /*
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // calculate engineering strain
    EngineeringStrain1D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);
    return 0.5 * engineeringStrain.mEngineeringStrain * this->mE * engineeringStrain.mEngineeringStrain;
*/
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
/*
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    EngineeringStress2D engineeringStress;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    const SectionBase* theSection(rElement->GetSection());
    if (theSection==0)
        throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] No section defined for element.");
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
        throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Plane stress is to be implemented.");
    }
    return 0.5*(
            engineeringStrain.mEngineeringStrain[0]*engineeringStress.mEngineeringStress[0]
           +engineeringStrain.mEngineeringStrain[1]*engineeringStress.mEngineeringStress[1]
           +engineeringStrain.mEngineeringStrain[2]*engineeringStress.mEngineeringStress[2]);
*/
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
/*
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
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
*/
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
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
void NuTo::Multiscale::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStrain1D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::Multiscale::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStrain2D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetDeltaElasticEngineeringStrain] not implemented.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::Multiscale::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetDeltaElasticEngineeringStrain] not implemented.");
}

#define MAXNUMNEWTONITERATIONS 20
#define PRINTRESULT true
#define MIN_DELTA_STRAIN_FACTOR 1e-7
void NuTo::Multiscale::Solve(const ElementBase* rElement, int rIp, const NuTo::DeformationGradient2D &rDeformationGradient, double rTolerance,
        std::stringstream& rStringStreamBeforeSolve, bool& rStringStreamBeforeSolveWritten)const
{
    const ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();
    StructureIp *fineScaleStructure = const_cast<StructureIp*>(staticData->GetFineScaleStructure());

    // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    // write data to stringstream
    rStringStreamBeforeSolveWritten = false;
    //if there is a problem with memory this could be replaced by writing it to a file
    //std::cout << "sizeof string " << stringStreamBeforeSolve.str().length() << std::endl;

    // start analysis
    double deltaStrainFactor(1.0);
    double curStrainFactor(1.0);

    const ConstitutiveStaticDataMultiscale2DPlaneStrain *oldStaticData = rElement->GetStaticData(rIp)->AsMultiscale2DPlaneStrain();
    EngineeringStrain2D prevStrain(oldStaticData->GetPrevStrain());

    EngineeringStrain2D totalEngineeringStrain(engineeringStrain);

    EngineeringStrain2D deltaEngineeringStrain;
    deltaEngineeringStrain = engineeringStrain-prevStrain;

    //update conre mat
    fineScaleStructure->NodeBuildGlobalDofs();

    //update tmpstatic data with zero displacements
    fineScaleStructure->ElementTotalUpdateTmpStaticData();

    //init some auxiliary variables
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
    NuTo::FullMatrix<double> dispForceVector;
    NuTo::FullMatrix<double> intForceVector;
    NuTo::FullMatrix<double> extForceVector;
    NuTo::FullMatrix<double> rhsVector;

    //allocate solver
    NuTo::SparseDirectSolverMUMPS mySolver;
    mySolver.SetShowTime(false);

    //calculate stiffness
    fineScaleStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    {
        std::cout << "test of stiffness still included " << std::endl;
        NuTo::FullMatrix<double> stiffnessMatrixCSRVector2Full(stiffnessMatrixCSRVector2);
        //std::cout<<"stiffness matrix" << std::endl;
        //stiffnessMatrixCSRVector2Full.Info(10,3);
        double interval(1e-7);
        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
        NuTo::FullMatrix<double> stiffnessMatrixCSRVector2_CDF(stiffnessMatrixCSRVector2.GetNumRows(), stiffnessMatrixCSRVector2.GetNumColumns());
        NuTo::FullMatrix<double> intForceVector1, intForceVector2, intForceVectorCDF(stiffnessMatrixCSRVector2.GetNumRows(),1);
        double energy1,energy2;
        fineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        fineScaleStructure->ElementTotalUpdateTmpStaticData();
        fineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector1);
        energy1 = fineScaleStructure->ElementTotalGetTotalEnergy();
        for (int count=0; count<displacementsActiveDOFsCheck.GetNumRows(); count++)
        {
            displacementsActiveDOFsCheck(count,0)+=interval;
            fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
            fineScaleStructure->ElementTotalUpdateTmpStaticData();
            fineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector2);
            energy2 = fineScaleStructure->ElementTotalGetTotalEnergy();
            stiffnessMatrixCSRVector2_CDF.SetColumn(count,(intForceVector2-intForceVector1)*(1./interval));
            intForceVectorCDF(count,0) = (energy2-energy1)/interval;
            displacementsActiveDOFsCheck(count,0)-=interval;
        }
        fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        if ((stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).Abs().Max()>1e-2)
        {
            std::cout << "globalStiffnessMatrix algo" << std::endl;
            stiffnessMatrixCSRVector2Full.Info(10,3);
            std::cout<< std::endl << "globalStiffnessMatrix cdf" << std::endl;
            stiffnessMatrixCSRVector2_CDF.Info(10,3);
            std::cout<< std::endl << "error" << std::endl;
            (stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).Info(10);
            std::cout << "maximum error is " << (stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).Abs().Max() << std::endl;
            std::cout<< std::endl << "intForceVector algo" << std::endl;
            intForceVector1.Trans().Info(10);
            std::cout<< std::endl << "intForceVector cdf" << std::endl;
            intForceVectorCDF.Trans().Info(10);
            throw MechanicsException("[NuTo::Multiscale::Solve] Stiffness matrix is not correct.");
        }
        else
            std::cout << "stiffness is OK "<< std::endl;
    }
    //set the total strain and calculate from the existing crack opening the homogeneous strain
    EngineeringStrain2D curEngineeringStrain(prevStrain+deltaEngineeringStrain*curStrainFactor);
    fineScaleStructure->SetTotalEngineeringStrain(curEngineeringStrain);

    //update conre mat
    fineScaleStructure->NodeBuildGlobalDofs();

    //update displacements of all nodes according to the new conre mat
    {
        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
        fineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        fineScaleStructure->ElementTotalUpdateTmpStaticData();
    }

    // build global external load vector and RHS vector
    fineScaleStructure->BuildGlobalExternalLoadVector(extForceVector);
    fineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector);
    rhsVector = extForceVector + dispForceVector - intForceVector;

    //calculate absolute tolerance for matrix entries to be not considered as zero
    double maxValue, minValue, ToleranceZeroStiffness;
    stiffnessMatrixCSRVector2.Max(maxValue);
    stiffnessMatrixCSRVector2.Min(minValue);
    //std::cout << "min and max " << minValue << " , " << maxValue << std::endl;

    ToleranceZeroStiffness = (1e-14) * (fabs(maxValue)>fabs(minValue) ?  fabs(maxValue) : fabs(minValue));
    fineScaleStructure->SetToleranceStiffnessEntries(ToleranceZeroStiffness);
    int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
    int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
    //std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;


    //repeat until max displacement is reached
    bool convergenceStatusLoadSteps(false);
    int loadstep(1);
    NuTo::FullMatrix<double> displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged;
    fineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
    while (!convergenceStatusLoadSteps)
    {

        double normResidual(1);
        double maxResidual(1);
        int numNewtonIterations(0);
        double normRHS(1.);
        double alpha(1.);
        int convergenceStatus(0);
        //0 - not converged, continue Newton iteration
        //1 - converged
        //2 - stop iteration, decrease load step
        while(convergenceStatus==0)
        {
            numNewtonIterations++;

            if (numNewtonIterations>MAXNUMNEWTONITERATIONS)
            {
                if (PRINTRESULT)
                {
                    std::cout << "numNewtonIterations (" << numNewtonIterations << ") > MAXNUMNEWTONITERATIONS (" << MAXNUMNEWTONITERATIONS << ")" << std::endl;
                }
                convergenceStatus = 2; //decrease load step
                break;
            }

            normRHS = rhsVector.Norm();

            // solve
            NuTo::FullMatrix<double> deltaDisplacementsActiveDOFs;
            NuTo::FullMatrix<double> oldDisplacementsActiveDOFs;
            NuTo::FullMatrix<double> displacementsActiveDOFs;
            NuTo::FullMatrix<double> displacementsDependentDOFs;
            NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrixCSRVector2);
            stiffnessMatrixCSR.SetOneBasedIndexing();
            mySolver.Solve(stiffnessMatrixCSR, rhsVector, deltaDisplacementsActiveDOFs);

            // write displacements to node
            fineScaleStructure->NodeExtractDofValues(oldDisplacementsActiveDOFs, displacementsDependentDOFs);

            //perform a linesearch
            alpha = 1.;
            do
            {
                //add new displacement state
                displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha;
                deltaDisplacementsActiveDOFs.Trans().Info(10,3);
                fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFs);
                fineScaleStructure->ElementTotalUpdateTmpStaticData();

                // calculate residual
                fineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector);
                rhsVector = extForceVector - intForceVector;
                normResidual = rhsVector.Norm();
                std::cout << "alpha " << alpha << ", normResidual " << normResidual << ", normResidualInit "<< normRHS << ", normRHS*(1-0.5*alpha) " << normRHS*(1-0.5*alpha) << std::endl;
                alpha*=0.5;
            }
            while(alpha>1e-3 && normResidual>normRHS*(1-0.5*alpha) && normResidual>rTolerance);
            if (normResidual>normRHS*(1-0.5*alpha) && normResidual>rTolerance)
            {
                convergenceStatus=2;
                break;
            }

            maxResidual = rhsVector.Max();

            //std::cout << std::endl << "Newton iteration " << numNewtonIterations << ", final alpha " << 2*alpha << ", normResidual " << normResidual<< ", maxResidual " << maxResidual<<std::endl;

            //check convergence
            if (normResidual<rTolerance || maxResidual<rTolerance)
            {
                if (PRINTRESULT)
                {
                    std::cout << "Convergence after " << numNewtonIterations << " Newton iterations, curStrainFactor " << curStrainFactor << ", deltaStrainFactor "<< deltaStrainFactor << std::endl<< std::endl;
                }
                convergenceStatus=1;
                break;
            }

            //convergence status == 0 (continue Newton iteration)
            //build new stiffness matrix
            fineScaleStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
            int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;
        }

        if (deltaStrainFactor<1e-7)
            throw NuTo::MechanicsException("[NuTo::Multiscale::Solve] No convergence, delta strain factor < 1e-7");

        if (convergenceStatus==1)
        {
            // visualize results

//#ifdef ENABLE_VISUALIZE
//            std::cout << " store element id and ip in output file" << std::endl;
//            this->ExportVtkDataFile(mIPName+std::string(".vtk"));
//#endif
//#ifdef ENABLE_SERIALIZATION
//            this->Save(mIPName+std::string("bin"),"BINARY");
//#endif // ENABLE_SERIALIZATION

            fineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
            if (curStrainFactor==1)
            {
                convergenceStatusLoadSteps=true;
            }
            else
            {
                if (rStringStreamBeforeSolveWritten==false)
                {
                    fineScaleStructure->Save(rStringStreamBeforeSolve);
                    rStringStreamBeforeSolveWritten=true;
                }

                // the update is only required to allow for a stepwise solution procedure in the fine scale model
                // a final update is only required for an update on the macroscale, otherwise,the original state has
                // to be reconstructed.
                fineScaleStructure->ElementTotalUpdateStaticData();

                //eventually increase load step
                if (numNewtonIterations<MAXNUMNEWTONITERATIONS/3)
                {
                    deltaStrainFactor*=1.5;
                }

                //increase displacement
                curStrainFactor+=deltaStrainFactor;
                if (curStrainFactor>1)
                {
                    deltaStrainFactor -= curStrainFactor -1.;
                    curStrainFactor=1;
                }

                curEngineeringStrain = prevStrain + deltaEngineeringStrain * curStrainFactor;

                //old stiffness matrix is used in first step of next load increment in order to prevent spurious problems at the boundary
                //std::cout << "press enter to next load increment, delta strain factor " << deltaStrainFactor << " max delta strain factor " <<  maxDeltaStrainFactor << std::endl << std::endl;
                //char cDummy[100]="";
                //std::cin.getline(cDummy, 100);
            }
            loadstep++;
        }
        else
        {
            assert(convergenceStatus==2);
            //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
            //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
            //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
            curStrainFactor-=deltaStrainFactor;
            curEngineeringStrain = prevStrain + deltaEngineeringStrain * curStrainFactor;

            //set the total strain and calculate from the existing crack opening and the total strain the homogeneous strain
            fineScaleStructure->SetTotalEngineeringStrain(curEngineeringStrain);

            // build global dof numbering
            fineScaleStructure->NodeBuildGlobalDofs();

            //set previous converged displacements
            fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
            fineScaleStructure->ElementTotalUpdateTmpStaticData();

            //decrease load step
            deltaStrainFactor*=0.5;
            curStrainFactor+=deltaStrainFactor;
            curEngineeringStrain = prevStrain + deltaEngineeringStrain * curStrainFactor;

            //check for minimum delta (this mostly indicates an error in the software
            if (deltaStrainFactor<MIN_DELTA_STRAIN_FACTOR)
            {
                deltaStrainFactor = 0;
                //throw NuTo::MechanicsException("Example ConcurrentMultiscale : No convergence, delta strain factor < 1e-7");
            }

            //std::cout << "press enter to reduce load increment" << std::endl;
            //char cDummy[100]="";
            //std::cin.getline(cDummy, 100);;
        }

        if (!convergenceStatusLoadSteps)
        {
            //update new displacement of RHS
            //set the total strain and calculate from the existing crack opening and the total strain the homogeneous strain
            fineScaleStructure->SetTotalEngineeringStrain(curEngineeringStrain);

            // build global dof numbering
            fineScaleStructure->NodeBuildGlobalDofs();

            //update stiffness in order to calculate new dispForceVector
            fineScaleStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
            int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            //std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;

            //update rhs vector for next Newton iteration
            rhsVector = dispForceVector + extForceVector - intForceVector;

            //update displacements of all nodes according to the new conre mat
            NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
            NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
            fineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
            fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
            fineScaleStructure->ElementTotalUpdateTmpStaticData();

            // calculate initial residual for next load step
            fineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector);

            //update rhs vector for next Newton iteration
            rhsVector = dispForceVector + extForceVector - intForceVector;
        }
    }
}
///////////////////////////////////////////////////////////////////////////

//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::Multiscale::GetType() const
{
    return NuTo::Constitutive::MULTISCALE;
}

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::Multiscale::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::PLANE2D3N:
        return true;
    case NuTo::Element::PLANE2D4N:
        return true;
    case NuTo::Element::PLANE2D6N:
        return true;
    default:
        return false;
    }
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::Multiscale::Info(unsigned short rVerboseLevel) const
{
}

// check parameters
void NuTo::Multiscale::CheckParameters()const
{
}

//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool NuTo::Multiscale::HaveTmpStaticData() const
{
    return false;
}

//! @brief ... returns true, if a material model has is nonlocal (stiffness is of dynamic size, nonlocal averaging)
//! @return ... see brief explanation
bool NuTo::Multiscale::IsNonlocalModel()const
{
    return false;
}
