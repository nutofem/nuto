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
#include "nuto/math/SparseMatrixCSRVector2General.h"
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
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"

#define sqrt3 1.732050808
#define MAX_OMEGA 0.999
//#define ENABLE_DEBUG

NuTo::Multiscale::Multiscale() : ConstitutiveEngineeringStressStrain()
{
    SetParametersValid();
    mToleranceResidualForce=1e-5;
    mMaxDeltaLoadFactor = 1;
    mMaxNumNewtonIterations=20;
    mDecreaseFactor=0.5;
    mMinNumNewtonIterations=7;
    mIncreaseFactor=1.5;
    mMinLoadFactor=1e-3;
    mMinLineSearchFactor=1e-3;
    mAugmentedLagrangeStiffnessCrackOpening = 1e3;
    mCrackTransitionRadius = 0;

    mTensileStrength = 0;

    mScalingFactorCrackAngle = 0;
    mScalingFactorCrackOpening = 0;
    mScalingFactorEpsilon = 0;

    mPenaltyStiffnessCrackAngle = 0;
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
       ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveEngineeringStressStrain)
          & BOOST_SERIALIZATION_NVP(mElasticStiffness)
          & BOOST_SERIALIZATION_NVP(mFileName)
          & BOOST_SERIALIZATION_NVP(mAugmentedLagrangeStiffnessCrackOpening)
          & BOOST_SERIALIZATION_NVP(mTensileStrength)
          & BOOST_SERIALIZATION_NVP(mScalingFactorCrackAngle)
          & BOOST_SERIALIZATION_NVP(mScalingFactorCrackOpening)
          & BOOST_SERIALIZATION_NVP(mScalingFactorEpsilon)
          & BOOST_SERIALIZATION_NVP(mPenaltyStiffnessCrackAngle)
          & BOOST_SERIALIZATION_NVP(mCrackTransitionRadius)
          & BOOST_SERIALIZATION_NVP(mToleranceResidualForce)
          & BOOST_SERIALIZATION_NVP(mMaxDeltaLoadFactor)
          & BOOST_SERIALIZATION_NVP(mMaxNumNewtonIterations)
          & BOOST_SERIALIZATION_NVP(mDecreaseFactor)
          & BOOST_SERIALIZATION_NVP(mMinNumNewtonIterations)
          & BOOST_SERIALIZATION_NVP(mIncreaseFactor)
          & BOOST_SERIALIZATION_NVP(mMinLoadFactor)
          & BOOST_SERIALIZATION_NVP(mMinLineSearchFactor)
          & BOOST_SERIALIZATION_NVP(mResultDirectory)
          & BOOST_SERIALIZATION_NVP(mLoadStepMacro);
#ifdef DEBUG_SERIALIZATION
       std::cout << "finish serialize Multiscale" << std::endl;
#endif
    }
    BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Multiscale)
#endif // ENABLE_SERIALIZATION

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
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
    }
    // this is a somehow weird situation, since for a normal material law nothing should be changed
    // since the material law is a full structure whose bc change, this can either be implemented with a cast (as I did)
    // or by removing the const flag from all material routines (which I do not consider is good)
    const ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();
    if (staticData->NonlinearSolutionOn()==false)
    {
        // linear solution

        // calculate engineering strain
        EngineeringStrain2D engineeringStrain;
        rDeformationGradient.GetEngineeringStrain(engineeringStrain);

        // calculate Engineering stress
        rEngineeringStress.mEngineeringStress[0] = mElasticStiffness(0,0) * engineeringStrain.mEngineeringStrain[0] + mElasticStiffness(0,1) * engineeringStrain.mEngineeringStrain[1] + mElasticStiffness(0,2) * engineeringStrain.mEngineeringStrain[2];
        rEngineeringStress.mEngineeringStress[1] = mElasticStiffness(1,0) * engineeringStrain.mEngineeringStrain[0] + mElasticStiffness(1,1) * engineeringStrain.mEngineeringStrain[1] + mElasticStiffness(1,2) * engineeringStrain.mEngineeringStrain[2];
        rEngineeringStress.mEngineeringStress[2] = mElasticStiffness(2,0) * engineeringStrain.mEngineeringStrain[0] + mElasticStiffness(2,1) * engineeringStrain.mEngineeringStrain[1] + mElasticStiffness(2,2) * engineeringStrain.mEngineeringStrain[2];
    }
    else
    {
        //nonlinear solution
        StructureMultiscale *fineScaleStructure = const_cast<StructureMultiscale*>(staticData->GetFineScaleStructure());
        NuTo::FullMatrix<double> activeDOF, dependentDOF;

        //change the previous angle of the alpha-constraint
        fineScaleStructure->SetPrevCrackAngleElastic(staticData->GetPrevCrackAngleElastic());
        fineScaleStructure->SetPrevCrackAngle(staticData->GetPrevCrackAngle());

        //Get and set previous total strain
        EngineeringStrain2D prevStrain(staticData->GetPrevStrain());
        fineScaleStructure->SetPrevTotalEngineeringStrain(prevStrain);

        fineScaleStructure->SetLoadFactor(0);
        fineScaleStructure->NodeBuildGlobalDofs();
        fineScaleStructure->NodeExtractDofValues(activeDOF,dependentDOF);
        fineScaleStructure->NodeMergeActiveDofValues(activeDOF);

         //Get and set previous delta strain
        EngineeringStrain2D engineeringStrain;
        rDeformationGradient.GetEngineeringStrain(engineeringStrain);
        EngineeringStrain2D deltaStrain(engineeringStrain-staticData->GetPrevStrain());
        fineScaleStructure->SetDeltaTotalEngineeringStrain(deltaStrain);

        fineScaleStructure->SetNewtonRaphsonToleranceResidualForce(mToleranceResidualForce);
        fineScaleStructure->SetNewtonRaphsonAutomaticLoadStepControl(true);
        fineScaleStructure->SetNewtonRaphsonMaxDeltaLoadFactor(mMaxDeltaLoadFactor);
        fineScaleStructure->SetNewtonRaphsonMaxNumNewtonIterations(mMaxNumNewtonIterations);
        fineScaleStructure->SetNewtonRaphsonDecreaseFactor(mDecreaseFactor);
        fineScaleStructure->SetNewtonRaphsonMinNumNewtonIterations(mMinNumNewtonIterations);
        fineScaleStructure->SetNewtonRaphsonIncreaseFactor(mIncreaseFactor);
        fineScaleStructure->SetNewtonRaphsonMinDeltaLoadFactor(mMinLoadFactor);
        fineScaleStructure->SetNewtonRaphsonMinLineSearchFactor(mMinLineSearchFactor);

        std::stringstream saveStream;
        bool hasBeenSaved(false);

        fineScaleStructure->GetLogger().OpenFile();
        fineScaleStructure->GetLogger() << "\n"<< "*******************************************" << "\n";
        fineScaleStructure->GetLogger() << " GetEngineeringStressFromEngineeringStrain " << "\n";
        fineScaleStructure->GetLogger() << " engineering strain " <<  engineeringStrain.mEngineeringStrain[0] << " "
        		                            <<  engineeringStrain.mEngineeringStrain[1] << " "
        		                            <<  engineeringStrain.mEngineeringStrain[2]  << "\n";
        fineScaleStructure->GetLogger() << " deltaStrain strain " <<  deltaStrain.mEngineeringStrain[0] << " "
        		                            <<  deltaStrain.mEngineeringStrain[1] << " "
        		                            <<  deltaStrain.mEngineeringStrain[2] << "\n";
        fineScaleStructure->GetLogger() << " prevStrain strain "  <<  staticData->GetPrevStrain().mEngineeringStrain[0] << " "
        		                            <<  staticData->GetPrevStrain().mEngineeringStrain[1] << " "
        		                            <<  staticData->GetPrevStrain().mEngineeringStrain[2] << "\n";
        try
        {
            fineScaleStructure->NewtonRaphson(true, saveStream, hasBeenSaved);
        }
        catch(NuTo::MechanicsException& e)
        {
        	//restore structure
            if (hasBeenSaved)
            {
                fineScaleStructure->RestoreStructure(saveStream);
            }
            else
            {
                //set load factor to zero in order to get the same ordering of the displacements as before the routine
                fineScaleStructure->SetLoadFactor(0);
                fineScaleStructure->NodeBuildGlobalDofs();
                fineScaleStructure->NodeMergeActiveDofValues(activeDOF);
            }
            std::cout << std::string("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName() << std::endl;
            e.AddMessage(std::string("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName());
        	throw e;
        }
        catch(...)
        {
        	throw MechanicsException(std::string("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName());
        }

        //calculate average stress
        NuTo::FullMatrix<double> averageStressDamage, averageStressHomogeneous;
        fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsDamage(),fineScaleStructure->GetAreaDamage(), averageStressDamage);
        fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsHomogeneous(),fineScaleStructure->GetAreaHomogeneous(), averageStressHomogeneous);
        double scalingFactorDamage = fineScaleStructure->GetScalingFactorDamage();
        double scalingFactorHomogeneous = fineScaleStructure->GetScalingFactorHomogeneous();
        double sum(scalingFactorDamage+scalingFactorHomogeneous);
        scalingFactorDamage/=sum;
        scalingFactorHomogeneous/=sum;
        averageStressDamage*=scalingFactorDamage;
        averageStressHomogeneous*=scalingFactorHomogeneous;
        rEngineeringStress.mEngineeringStress[0] = averageStressDamage(0,0)+averageStressHomogeneous(0,0);
        rEngineeringStress.mEngineeringStress[1] = averageStressDamage(1,0)+averageStressHomogeneous(1,0);
        rEngineeringStress.mEngineeringStress[2] = averageStressDamage(3,0)+averageStressHomogeneous(3,0);

        //restore previous state (only performed if the load had to be subdivided)
        if (hasBeenSaved)
        {
            fineScaleStructure->RestoreStructure(saveStream);
        }
        else
        {
            //set load factor to zero in order to get the same ordering of the displacements as before the routine
            fineScaleStructure->SetLoadFactor(0);
            fineScaleStructure->NodeBuildGlobalDofs();
            fineScaleStructure->NodeMergeActiveDofValues(activeDOF);
        }
        fineScaleStructure->GetLogger().CloseFile();
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
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] Check the material parameters.");
    }
    EngineeringStress2D engineeringStress2D;
    GetEngineeringStressFromEngineeringStrain(rElement, rIp,rDeformationGradient, engineeringStress2D);
    //this is certainly not correct, since it is plane strain, but I do not care about the stress in thickness direction
    rEngineeringStress.mEngineeringStress[0] = engineeringStress2D.mEngineeringStress[0];
    rEngineeringStress.mEngineeringStress[1] = engineeringStress2D.mEngineeringStress[1];
    rEngineeringStress.mEngineeringStress[2] = 0.;
    rEngineeringStress.mEngineeringStress[3] = engineeringStress2D.mEngineeringStress[2];
    rEngineeringStress.mEngineeringStress[4] = 0.;
    rEngineeringStress.mEngineeringStress[5] = 0.;
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
	rDamage=0.;
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
	rDamage=0.;
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
	rDamage=0.;
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
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }
#ifdef HAVE_MUMPS
    const ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();
    if (staticData->NonlinearSolutionOn()==false)
    {
        //linear elastic solution
        ConstitutiveTangentLocal3x3 *tangent(rTangent->AsConstitutiveTangentLocal3x3());

        // store tangent at the output object
        tangent->mTangent[ 0] = mElasticStiffness(0,0);
        tangent->mTangent[ 1] = mElasticStiffness(1,0);
        tangent->mTangent[ 2] = mElasticStiffness(2,0);

        tangent->mTangent[ 3] = mElasticStiffness(0,1);
        tangent->mTangent[ 4] = mElasticStiffness(1,1);
        tangent->mTangent[ 5] = mElasticStiffness(2,1);

        tangent->mTangent[ 6] = mElasticStiffness(0,2);
        tangent->mTangent[ 7] = mElasticStiffness(1,2);
        tangent->mTangent[ 8] = mElasticStiffness(2,2);
    }
    else
    {
        //nonlinear solution
        const ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();
        StructureMultiscale *fineScaleStructure = const_cast<StructureMultiscale*>(staticData->GetFineScaleStructure());
        fineScaleStructure->GetLogger().OpenFile();

        //change the previous angle of the alpha-constraint
        fineScaleStructure->SetPrevCrackAngleElastic(staticData->GetPrevCrackAngleElastic());
        fineScaleStructure->SetPrevCrackAngle(staticData->GetPrevCrackAngle());

        //Get and set previous total strain
        EngineeringStrain2D prevStrain(staticData->GetPrevStrain());
        fineScaleStructure->SetPrevTotalEngineeringStrain(prevStrain);

        fineScaleStructure->SetLoadFactor(0);
        fineScaleStructure->NodeBuildGlobalDofs();
        NuTo::FullMatrix<double> activeDOF, dependentDOF;
        fineScaleStructure->NodeExtractDofValues(activeDOF,dependentDOF);
        fineScaleStructure->NodeMergeActiveDofValues(activeDOF);

        //Get and set previous delta strain
        EngineeringStrain2D engineeringStrain;
        rDeformationGradient.GetEngineeringStrain(engineeringStrain);
        EngineeringStrain2D deltaStrain(engineeringStrain-staticData->GetPrevStrain());
        fineScaleStructure->SetDeltaTotalEngineeringStrain(deltaStrain);

        fineScaleStructure->SetNewtonRaphsonToleranceResidualForce(mToleranceResidualForce);
        fineScaleStructure->SetNewtonRaphsonAutomaticLoadStepControl(true);
        fineScaleStructure->SetNewtonRaphsonMaxDeltaLoadFactor(mMaxDeltaLoadFactor);
        fineScaleStructure->SetNewtonRaphsonMaxNumNewtonIterations(mMaxNumNewtonIterations);
        fineScaleStructure->SetNewtonRaphsonDecreaseFactor(mDecreaseFactor);
        fineScaleStructure->SetNewtonRaphsonMinNumNewtonIterations(mMinNumNewtonIterations);
        fineScaleStructure->SetNewtonRaphsonIncreaseFactor(mIncreaseFactor);
        fineScaleStructure->SetNewtonRaphsonMinDeltaLoadFactor(mMinLoadFactor);
        fineScaleStructure->SetNewtonRaphsonMinLineSearchFactor(mMinLineSearchFactor);

        std::stringstream saveStream;
        bool hasBeenSaved(false);

        try
        {
            fineScaleStructure->GetLogger() << "\n" << "************************************************" << "\n";
            fineScaleStructure->GetLogger() << " GetTangent_EngineeringStress_EngineeringStrain " << "\n";
            fineScaleStructure->GetLogger() << " engineering strain " <<  engineeringStrain.mEngineeringStrain[0] << " "
            		                            <<  engineeringStrain.mEngineeringStrain[1] << " "
            		                            <<  engineeringStrain.mEngineeringStrain[2] << "\n";
            fineScaleStructure->GetLogger() << " deltaStrain strain " <<  deltaStrain.mEngineeringStrain[0] << " "
            		                            <<  deltaStrain.mEngineeringStrain[1] << " "
            		                            <<  deltaStrain.mEngineeringStrain[2] << "\n";
            fineScaleStructure->GetLogger() << " prevStrain strain "  <<  staticData->GetPrevStrain().mEngineeringStrain[0] << " "
            		                            <<  staticData->GetPrevStrain().mEngineeringStrain[1] << " "
            		                            <<  staticData->GetPrevStrain().mEngineeringStrain[2] << "\n"  << "\n";
            fineScaleStructure->NewtonRaphson(true, saveStream, hasBeenSaved);
        }
        catch(NuTo::MechanicsException& e)
        {
            //restore structure
            if (hasBeenSaved)
            {
                fineScaleStructure->RestoreStructure(saveStream);
            }
            else
            {
                //set load factor to zero in order to get the same ordering of the displacements as before the routine
                fineScaleStructure->SetLoadFactor(0);
                fineScaleStructure->NodeBuildGlobalDofs();
                fineScaleStructure->NodeMergeActiveDofValues(activeDOF);
            }
            std::cout << std::string("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName() << std::endl;
            e.AddMessage(std::string("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] No convergence in multiscale for ip ") + fineScaleStructure->GetIPName());
        	throw e;
        }
        catch(...)
        {
        	throw MechanicsException(std::string("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName());
        }

        // use Schur complement to calculate the stiffness

        //no change the constraints the way that the total strain is no longer fixed
        ConstraintBase* linearTotalStrainConstraintPtr = fineScaleStructure->ConstraintRelease(fineScaleStructure->GetConstraintTotalStrain()); // now this is not deleted
        //std::cout << "constraint for alpha is deleted for calculation of stiffness " << std::endl;
        //ConstraintBase* nonlinearConstraintAlphaPtr = fineScaleStructure->ConstraintRelease(fineScaleStructure->GetConstraintCrackAngle()); // now this is not deleted
        fineScaleStructure->NodeBuildGlobalDofs();
        SparseMatrixCSRVector2General<double>
            matrixJJ(fineScaleStructure->GetNumActiveDofs(), fineScaleStructure->GetNumActiveDofs());
        FullMatrix<double> rhsVector(fineScaleStructure->GetNumDofs() - fineScaleStructure->GetNumActiveDofs(),1);
        fineScaleStructure->BuildGlobalCoefficientMatrix0(matrixJJ, rhsVector);
        if (rhsVector.Abs().Max()>1e-6)
        {
        	throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] RHS vector not zero, previous equilibrium iteration was not successfull.");
        }

    /*
        const_cast<StructureMultiscale*> (fineScaleStructure)->ElementTotalUpdateTmpStaticData();
        double totEnergy = fineScaleStructure->ElementTotalGetTotalEnergy();
        const_cast<StructureMultiscale*> (fineScaleStructure)->CalculateHomogeneousEngineeringStrain();
        //calculate average strain
        NuTo::FullMatrix<double> averageStrain,averageStress;
        double mlX, mlY;
        mlX = fineScaleStructure->GetDimensionX();
        mlY = fineScaleStructure->GetDimensionY();
        fineScaleStructure->ElementTotalGetAverageStrain(mlX*mlY,averageStrain);
        std::cout<<"average strain " << std::endl;
        averageStrain.Info(12,4);
        fineScaleStructure->ElementTotalGetAverageStress(mlX*mlY,averageStress);
        std::cout<<"average stress " << std::endl;
        averageStress.Info(12,4);
        fineScaleStructure->NodeInfo(10);
        fineScaleStructure->ExportVtkDataFile("/home/unger3/develop/nuto_build/examples/c++/FineScaleConcurrentMultiscale.vtk");
    */
/*        std::cout << "matrixJJ algo "<< std::endl;
        NuTo::FullMatrix<double> matrixJJFull(matrixJJ);
        matrixJJFull.Info(12,3);

        {
        //check the Global Matrix
		NuTo::FullMatrix<double> matrixJJCDF(fineScaleStructure->GetNumActiveDofs(),fineScaleStructure->GetNumActiveDofs());
		//std::cout<<"stiffness matrix" << std::endl;
		//stiffnessMatrixCSRVector2Full.Info(10,3);
		double interval(1e-9);
		NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
		NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
		StructureMultiscale *fineScaleStructureNonConst = const_cast<StructureMultiscale*> (fineScaleStructure);

		fineScaleStructureNonConst->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
		fineScaleStructureNonConst->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
		fineScaleStructureNonConst->ElementTotalUpdateTmpStaticData();

		NuTo::FullMatrix<double> activeGrad1(displacementsActiveDOFsCheck),
								 activeGrad2(displacementsActiveDOFsCheck);
		fineScaleStructureNonConst->BuildGlobalGradientInternalPotentialVector(activeGrad1);
		for (int count=0; count<displacementsActiveDOFsCheck.GetNumRows(); count++)
		{
			displacementsActiveDOFsCheck(count,0)+=interval;
			fineScaleStructureNonConst->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
			fineScaleStructureNonConst->ElementTotalUpdateTmpStaticData();
			fineScaleStructure->BuildGlobalGradientInternalPotentialVector(activeGrad2);
			matrixJJCDF.SetColumn(count,(activeGrad2-activeGrad1)*(1./interval));

			displacementsActiveDOFsCheck(count,0)-=interval;
		}
		fineScaleStructureNonConst->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);

		std::cout << "matrixJJ cdf "<< std::endl;
		matrixJJCDF.Info(12,3);

		std::cout << "delta "<< std::endl;
		(matrixJJFull-matrixJJCDF).Info(12,3);

		std::cout << "maxError JJ " << (matrixJJCDF-NuTo::FullMatrix<double>(matrixJJ)).Abs().Max() << std::endl;
        }
*/
/*        {
        std::cout << "stiffness before solution" << std::endl;
        NuTo::FullMatrix<double> matrixJJFull(matrixJJ);
        matrixJJFull.Info(12,3);
        NuTo::FullMatrix<double> eigenvalues, eigenvectors;
        matrixJJFull.EigenValuesSymmetric(eigenvalues);
        matrixJJFull.EigenVectorsSymmetric(eigenvectors);
        std::cout << "eigenvalues " << std::endl;
        eigenvalues.Trans().Info(12,3);
        std::cout << "eigenvectors " << std::endl;
        eigenvectors.Info(12,3);
        }
*/
        //calculate schur complement
        NuTo::FullMatrix<int> schurIndicesMatrix(3,1);
        NuTo::FullMatrix<double> stiffness(3,3);
        //attention, the index is in the zero based indexing system
        schurIndicesMatrix(0,0) = fineScaleStructure->GetDofGlobalTotalStrain2D()[0];
        schurIndicesMatrix(1,0) = fineScaleStructure->GetDofGlobalTotalStrain2D()[1];
        schurIndicesMatrix(2,0) = fineScaleStructure->GetDofGlobalTotalStrain2D()[2];

        NuTo::SparseDirectSolverMUMPS mumps;
#ifdef SHOW_TIME
        mumps.SetShowTime(false);
#endif
        NuTo::SparseMatrixCSRGeneral<double> stiffnessFineScale(matrixJJ);
        stiffnessFineScale.SetOneBasedIndexing();
        mumps.SchurComplement(stiffnessFineScale,schurIndicesMatrix,stiffness);
        //scale with the dimension of the structure (area)
        double area = fineScaleStructure->GetAreaDamage()*fineScaleStructure->GetScalingFactorDamage()+
        		      fineScaleStructure->GetAreaHomogeneous()*fineScaleStructure->GetScalingFactorHomogeneous();
        // scale the stiffness in order to account for the difference between dof and epsilon
        double scalingFactorEpsilon(fineScaleStructure->GetScalingFactorEpsilon());
        stiffness*=1./(area*scalingFactorEpsilon*scalingFactorEpsilon);
        *(rTangent->AsConstitutiveTangentLocal3x3()) = stiffness;

        //restore previous state (only performed if the load had to be subdivided)
        if (hasBeenSaved)
        {
        	fineScaleStructure->RestoreStructure(saveStream);
        }
        else
        {
            //reinsert the total strain constraint set load factor to zero in order to get the same ordering of the displacements as before the routine
            fineScaleStructure->ConstraintAdd(fineScaleStructure->GetConstraintTotalStrain(),linearTotalStrainConstraintPtr);
            //fineScaleStructure->ConstraintAdd(fineScaleStructure->GetConstraintCrackAngle(),nonlinearConstraintAlphaPtr);
            fineScaleStructure->SetLoadFactor(0);
            fineScaleStructure->NodeBuildGlobalDofs();
            fineScaleStructure->NodeMergeActiveDofValues(activeDOF);
        }
		//std::cout << "crack opening "<< fineScaleStructure->GetGlobalCrackOpening2D()[0] << " " << fineScaleStructure->GetGlobalCrackOpening2D()[1] << std::endl;
		//std::cout << "crack angle "<< fineScaleStructure->GetCrackAngle()*180./M_PI << ", crack angle elastic "<< fineScaleStructure->GetCrackAngleElastic()*180./M_PI  << std::endl;
        fineScaleStructure->GetLogger()<< "stiffness algo " << "\n";
        fineScaleStructure->GetLogger().Out(stiffness,12,3,false);
		//std::cout << "stiffness" << "\n";
        //stiffness.Info(12,3);
        if (1==0)
        {
			//just for test purpose, calculate stiffness via resforces
			EngineeringStress2D stress1, stress2;
			GetEngineeringStressFromEngineeringStrain(rElement, rIp, rDeformationGradient, stress1);
			double delta(1e-8);
			NuTo::FullMatrix<double> stiffnessCDF(3,3);
			for (int count=0; count<3; count++)
			{
			   DeformationGradient2D deformationGradient(rDeformationGradient);
			   switch(count)
				{
				case 0:
					deformationGradient.mDeformationGradient[0]+=delta;
					break;
				case 1:
					deformationGradient.mDeformationGradient[3]+=delta;
					break;
				case 2:
					deformationGradient.mDeformationGradient[1]+=0.5*delta;
					deformationGradient.mDeformationGradient[2]+=0.5*delta;
					break;
				default:
					throw MechanicsException("");
				}
				GetEngineeringStressFromEngineeringStrain(rElement, rIp, deformationGradient, stress2);

				for (int count2=0; count2<3; count2++)
				{
					stiffnessCDF(count2,count) = (stress2.GetData()[count2]- stress1.GetData()[count2])/delta;
				}
			}
			std::cout << "Schur stiffness algo" << std::endl;
			stiffness.Info(12,3);
			std::cout << "Schur stiffness cdf" << std::endl;
			stiffnessCDF.Info(12,3);
			if ((stiffness-stiffnessCDF).Abs().Max()>1e-4)
				std::cout << "Multiscale stiffness not correct" << std::endl;
			else
				std::cout << "Multiscale stiffness is correct" << std::endl;
            NuTo::FullMatrix<double> eigenValues;
            stiffness.EigenValuesSymmetric(eigenValues);
            std::cout << "eigenvalues" << std::endl;
            eigenValues.Trans().Info(12,3);
            NuTo::FullMatrix<double> eigenVectors;
            stiffness.EigenVectorsSymmetric(eigenVectors);
            std::cout << "eigenvector 1" << std::endl;
            eigenVectors.GetColumn(0).Trans().Info(12,3);
            std::cout << "eigenvector 2" << std::endl;
            eigenVectors.GetColumn(1).Trans().Info(12,3);
            std::cout << "eigenvector 3" << std::endl;
            eigenVectors.GetColumn(2).Trans().Info(12,3);
        }
        fineScaleStructure->GetLogger().CloseFile();
    }//nonlinear solution
#else
    throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] MUMPS solver required to calculate Schur complement.");
#endif //HAVE_MUMPS
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
        const DeformationGradient1D& rDeformationGradient)const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient)const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }
    // this is a somehow weird situation, since for a normal material law nothing should be changed
    // since the material law is a full structure whose bc change, this can either be implemented with a cast (as I did)
    // or by removing the const flag from all material routines (which I do not consider is good)
    ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();

    // linear solution - no static data to be updated, but strain is set as the total strain of the fine scale model
    EngineeringStrain2D engineeringStrain;
    EngineeringStrain2D prevStrain(staticData->GetPrevStrain());
    EngineeringStress2D meanEngineeringStress, engineeringStress;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);
    if (staticData->NonlinearSolutionOn()==false)
    {
        // calculate Engineering stress
        engineeringStress.mEngineeringStress[0] = mElasticStiffness(0,0) * engineeringStrain.mEngineeringStrain[0] + mElasticStiffness(0,1) * engineeringStrain.mEngineeringStrain[1] + mElasticStiffness(0,2) * engineeringStrain.mEngineeringStrain[1];
        engineeringStress.mEngineeringStress[1] = mElasticStiffness(1,0) * engineeringStrain.mEngineeringStrain[0] + mElasticStiffness(1,1) * engineeringStrain.mEngineeringStrain[1] + mElasticStiffness(1,2) * engineeringStrain.mEngineeringStrain[1];
        engineeringStress.mEngineeringStress[2] = mElasticStiffness(2,0) * engineeringStrain.mEngineeringStrain[0] + mElasticStiffness(2,1) * engineeringStrain.mEngineeringStrain[1] + mElasticStiffness(2,2) * engineeringStrain.mEngineeringStrain[1];
    }
    else
    {
        StructureMultiscale *fineScaleStructure = const_cast<StructureMultiscale*>(staticData->GetFineScaleStructure());
        fineScaleStructure->GetLogger().OpenFile();

        //change the previous angle of the alpha-constraint and the previous crack angle
        fineScaleStructure->SetPrevCrackAngleElastic(staticData->GetPrevCrackAngleElastic());
        fineScaleStructure->SetPrevCrackAngle(staticData->GetPrevCrackAngle());

        //Get and set previous total strain
        fineScaleStructure->SetPrevTotalEngineeringStrain(prevStrain);

        fineScaleStructure->SetLoadFactor(0);
        fineScaleStructure->NodeBuildGlobalDofs();
        NuTo::FullMatrix<double> activeDOF, dependentDOF;
        fineScaleStructure->NodeExtractDofValues(activeDOF,dependentDOF);
        fineScaleStructure->NodeMergeActiveDofValues(activeDOF);

        // calculate engineering strain
        EngineeringStrain2D deltaStrain(engineeringStrain-prevStrain);
        fineScaleStructure->SetDeltaTotalEngineeringStrain(deltaStrain);

        fineScaleStructure->SetNewtonRaphsonToleranceResidualForce(mToleranceResidualForce);
        fineScaleStructure->SetNewtonRaphsonAutomaticLoadStepControl(true);
        fineScaleStructure->SetNewtonRaphsonMaxDeltaLoadFactor(mMaxDeltaLoadFactor);
        fineScaleStructure->SetNewtonRaphsonMaxNumNewtonIterations(mMaxNumNewtonIterations);
        fineScaleStructure->SetNewtonRaphsonDecreaseFactor(mDecreaseFactor);
        fineScaleStructure->SetNewtonRaphsonMinNumNewtonIterations(mMinNumNewtonIterations);
        fineScaleStructure->SetNewtonRaphsonIncreaseFactor(mIncreaseFactor);
        fineScaleStructure->SetNewtonRaphsonMinDeltaLoadFactor(mMinLoadFactor);
        fineScaleStructure->SetNewtonRaphsonMinLineSearchFactor(mMinLineSearchFactor);

        std::stringstream saveStream;
        bool hasBeenSaved(false);

        fineScaleStructure->GetLogger() << "\n" << "******************************************************" << "\n";
        fineScaleStructure->GetLogger() << " UpdateStaticData_EngineeringStress_EngineeringStrain " << "\n";
        fineScaleStructure->GetLogger() << " engineering strain " <<  engineeringStrain.mEngineeringStrain[0] << " "
        		                            <<  engineeringStrain.mEngineeringStrain[1] << " "
        		                            <<  engineeringStrain.mEngineeringStrain[2]  << "\n";
        fineScaleStructure->GetLogger() << " deltaStrain strain " <<  deltaStrain.mEngineeringStrain[0] << " "
        		                            <<  deltaStrain.mEngineeringStrain[1] << " "
        		                            <<  deltaStrain.mEngineeringStrain[2] <<  "\n";
        fineScaleStructure->GetLogger() << " prevStrain strain "  <<  staticData->GetPrevStrain().mEngineeringStrain[0] << " "
        		                            <<  staticData->GetPrevStrain().mEngineeringStrain[1] << " "
        		                            <<  staticData->GetPrevStrain().mEngineeringStrain[2]  << "\n"  << "\n";
        try
        {
            fineScaleStructure->NewtonRaphson(true, saveStream, hasBeenSaved);
        }
        catch(NuTo::MechanicsException& e)
        {
            //restore structure
            if (hasBeenSaved)
            {
                fineScaleStructure->RestoreStructure(saveStream);
            }
            else
            {
                //set load factor to zero in order to get the same ordering of the displacements as before the routine
                fineScaleStructure->SetLoadFactor(0);
                fineScaleStructure->NodeBuildGlobalDofs();
                fineScaleStructure->NodeMergeActiveDofValues(activeDOF);
            }
            e.AddMessage(std::string("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] Error in performing Newton-iteration on fine scale for ip") + fineScaleStructure->GetIPName());
        	throw e;
        }
        catch(...)
        {
        	throw MechanicsException(std::string("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] Error in performing Newton-iteration on fine scale for ip") + fineScaleStructure->GetIPName());
        }

        //calculate delta total energy (sigma1+sigma2)/2*delta_strain
        //calculate average stress
        NuTo::FullMatrix<double> averageStressDamage, averageStressHomogeneous;
        fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsDamage(),fineScaleStructure->GetAreaDamage(), averageStressDamage);
        fineScaleStructure->GetLogger() << "average stress in damage domain"  << "\n";
        fineScaleStructure->GetLogger().Out(averageStressDamage.Trans(),12,4);
        fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsHomogeneous(),fineScaleStructure->GetAreaHomogeneous(), averageStressHomogeneous);
        fineScaleStructure->GetLogger() << "average stress in homogeneous domain"  << "\n";
        fineScaleStructure->GetLogger().Out(averageStressHomogeneous.Trans(),12,4);
        double scalingFactorDamage = fineScaleStructure->GetScalingFactorDamage();
        double scalingFactorHomogeneous = fineScaleStructure->GetScalingFactorHomogeneous();
        double sum(scalingFactorDamage+scalingFactorHomogeneous);
        scalingFactorDamage/=sum;
        scalingFactorHomogeneous/=sum;
        averageStressDamage*=scalingFactorDamage;
        averageStressHomogeneous*=scalingFactorHomogeneous;

        engineeringStress.mEngineeringStress[0] = (averageStressDamage(0,0)+averageStressHomogeneous(0,0));
        engineeringStress.mEngineeringStress[1] = (averageStressDamage(1,0)+averageStressHomogeneous(1,0));
        engineeringStress.mEngineeringStress[2] = (averageStressDamage(3,0)+averageStressHomogeneous(3,0));

        //fineScaleStructure->GetLogger() << "Energy of fine scale exact  "  << fineScaleStructure->ElementTotalGetTotalEnergy()  << "\n";
        //fineScaleStructure->GetLogger() << "Energy of fine scale approx "  << energy * (fineScaleStructure->GetAreaDamage()*fineScaleStructure->GetScalingFactorDamage()+
        //		                                                  fineScaleStructure->GetAreaHomogeneous()*fineScaleStructure->GetScalingFactorHomogeneous())  << "\n";
        //fineScaleStructure->GetLogger() << "total area of macroscale " <<  (fineScaleStructure->GetAreaDamage()*fineScaleStructure->GetScalingFactorDamage()+
        //        fineScaleStructure->GetAreaHomogeneous()*fineScaleStructure->GetScalingFactorHomogeneous())  << "\n";
        const_cast<ConstitutiveStaticDataMultiscale2DPlaneStrain*>(staticData)->SetPrevCrackAngleElastic(fineScaleStructure->GetCrackAngleElastic());
        std::cout << "actual crack angle " << fineScaleStructure->GetCrackAngle()*180./M_PI << "prev crack angle " << staticData->GetPrevCrackAngle()*180./M_PI << "\n";
        const_cast<ConstitutiveStaticDataMultiscale2DPlaneStrain*>(staticData)->SetPrevCrackAngle(fineScaleStructure->GetCrackAngle());
        fineScaleStructure->ElementTotalUpdateStaticData();
        fineScaleStructure->GetLogger().CloseFile();
    }
    double energy = staticData->GetPrevTotalEnergy();
    const EngineeringStress2D& prevStress(staticData->GetPrevStress());
    meanEngineeringStress.mEngineeringStress[0] = (engineeringStress.mEngineeringStress[0]+prevStress.mEngineeringStress[0]);
    meanEngineeringStress.mEngineeringStress[1] = (engineeringStress.mEngineeringStress[1]+prevStress.mEngineeringStress[1]);
    meanEngineeringStress.mEngineeringStress[2] = (engineeringStress.mEngineeringStress[2]+prevStress.mEngineeringStress[2]);

    energy+=0.5*(meanEngineeringStress.mEngineeringStress[0]*(engineeringStrain.mEngineeringStrain[0]-prevStrain.mEngineeringStrain[0])+
                 meanEngineeringStress.mEngineeringStress[1]*(engineeringStrain.mEngineeringStrain[1]-prevStrain.mEngineeringStrain[1])+
                 meanEngineeringStress.mEngineeringStress[2]*(engineeringStrain.mEngineeringStrain[2]-prevStrain.mEngineeringStrain[2]));
    const_cast<ConstitutiveStaticDataMultiscale2DPlaneStrain*>(staticData)->SetPrevStrain(engineeringStrain);
    const_cast<ConstitutiveStaticDataMultiscale2DPlaneStrain*>(staticData)->SetPrevStress(engineeringStress);
    const_cast<ConstitutiveStaticDataMultiscale2DPlaneStrain*>(staticData)->SetPrevTotalEnergy(energy);
}
//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient)const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient)const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient)const
{
    return;
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient)const
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
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    const ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();
    StructureMultiscale *fineScaleStructure = const_cast<StructureMultiscale*>(staticData->GetFineScaleStructure());
    fineScaleStructure->GetLogger().OpenFile();

    //Get and set previous total strain
     EngineeringStrain2D prevStrain(staticData->GetPrevStrain());
     fineScaleStructure->SetPrevTotalEngineeringStrain(prevStrain);

    fineScaleStructure->SetLoadFactor(0);
    fineScaleStructure->NodeBuildGlobalDofs();
    NuTo::FullMatrix<double> activeDOF, dependentDOF;
    fineScaleStructure->NodeExtractDofValues(activeDOF,dependentDOF);
    fineScaleStructure->NodeMergeActiveDofValues(activeDOF);

    // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);
    EngineeringStrain2D deltaStrain(engineeringStrain-prevStrain);
    fineScaleStructure->SetDeltaTotalEngineeringStrain(deltaStrain);

    //change the previous angle of the alpha-constraint
    fineScaleStructure->SetPrevCrackAngleElastic(staticData->GetPrevCrackAngleElastic());
    fineScaleStructure->SetPrevCrackAngle(staticData->GetPrevCrackAngle());

    fineScaleStructure->SetNewtonRaphsonToleranceResidualForce(mToleranceResidualForce);
    fineScaleStructure->SetNewtonRaphsonAutomaticLoadStepControl(true);
    fineScaleStructure->SetNewtonRaphsonMaxDeltaLoadFactor(mMaxDeltaLoadFactor);
    fineScaleStructure->SetNewtonRaphsonMaxNumNewtonIterations(mMaxNumNewtonIterations);
    fineScaleStructure->SetNewtonRaphsonDecreaseFactor(mDecreaseFactor);
    fineScaleStructure->SetNewtonRaphsonMinNumNewtonIterations(mMinNumNewtonIterations);
    fineScaleStructure->SetNewtonRaphsonIncreaseFactor(mIncreaseFactor);
    fineScaleStructure->SetNewtonRaphsonMinDeltaLoadFactor(mMinLoadFactor);
    fineScaleStructure->SetNewtonRaphsonMinLineSearchFactor(mMinLineSearchFactor);

    std::stringstream saveStream;
    bool hasBeenSaved(false);
    fineScaleStructure->GetLogger() << "\n" << "****************************************************" << "\n";
    fineScaleStructure->GetLogger() << " GetTotalEnergy_EngineeringStress_EngineeringStrain "  << "\n";
    fineScaleStructure->GetLogger() << " engineering strain " <<  engineeringStrain.mEngineeringStrain[0] << " "
    		                            <<  engineeringStrain.mEngineeringStrain[1] << " "
    		                            <<  engineeringStrain.mEngineeringStrain[2]  << "\n";
    fineScaleStructure->GetLogger() << " deltaStrain strain " <<  deltaStrain.mEngineeringStrain[0] << " "
    		                            <<  deltaStrain.mEngineeringStrain[1] << " "
    		                            <<  deltaStrain.mEngineeringStrain[2] << "\n";
    fineScaleStructure->GetLogger() << " prevStrain strain "  <<  staticData->GetPrevStrain().mEngineeringStrain[0] << " "
    		                            <<  staticData->GetPrevStrain().mEngineeringStrain[1] << " "
    		                            <<  staticData->GetPrevStrain().mEngineeringStrain[2]  << "\n" << "\n";
    try
    {
        fineScaleStructure->NewtonRaphson(true, saveStream, hasBeenSaved);
    }
    catch(MechanicsException& e)
    {
        //restore structure
        if (hasBeenSaved)
        {
            fineScaleStructure->RestoreStructure(saveStream);
        }
        else
        {
            //set load factor to zero in order to get the same ordering of the displacements as before the routine
            fineScaleStructure->SetLoadFactor(0);
            fineScaleStructure->NodeBuildGlobalDofs();
            fineScaleStructure->NodeMergeActiveDofValues(activeDOF);
        }
        e.AddMessage(std::string("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName());
    	throw e;
    }
    catch(...)
    {
    	throw MechanicsException(std::string("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName());
    }

    double energy = staticData->GetPrevTotalEnergy();
    //calculate delta total energy (sigma1+sigma2)/2*delta_strain
    //calculate average stress
    NuTo::FullMatrix<double> averageStressDamage, averageStressHomogeneous;
    fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsDamage(),fineScaleStructure->GetAreaDamage(), averageStressDamage);
    fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsHomogeneous(),fineScaleStructure->GetAreaHomogeneous(), averageStressHomogeneous);
    double scalingFactorDamage = fineScaleStructure->GetScalingFactorDamage();
    double scalingFactorHomogeneous = fineScaleStructure->GetScalingFactorHomogeneous();
    double sum(scalingFactorDamage+scalingFactorHomogeneous);
    scalingFactorDamage/=sum;
    scalingFactorHomogeneous/=sum;
    averageStressDamage*=scalingFactorDamage;
    averageStressHomogeneous*=scalingFactorHomogeneous;
    EngineeringStress2D meanEngineeringStress;
    const EngineeringStress2D& prevStress(staticData->GetPrevStress());
    meanEngineeringStress.mEngineeringStress[0] = (averageStressDamage(0,0)+averageStressHomogeneous(0,0)+prevStress.mEngineeringStress[0]);
    meanEngineeringStress.mEngineeringStress[1] = (averageStressDamage(1,0)+averageStressHomogeneous(1,0)+prevStress.mEngineeringStress[1]);
    meanEngineeringStress.mEngineeringStress[2] = (averageStressDamage(3,0)+averageStressHomogeneous(3,0)+prevStress.mEngineeringStress[2]);

    energy+=0.5*(meanEngineeringStress.mEngineeringStress[0]*(engineeringStrain.mEngineeringStrain[0]-prevStrain.mEngineeringStrain[0])+
                 meanEngineeringStress.mEngineeringStress[1]*(engineeringStrain.mEngineeringStrain[1]-prevStrain.mEngineeringStrain[1])+
                 meanEngineeringStress.mEngineeringStress[2]*(engineeringStrain.mEngineeringStrain[2]-prevStrain.mEngineeringStrain[2]));
    //restore structure
    if (hasBeenSaved)
    {
        fineScaleStructure->RestoreStructure(saveStream);
    }
    else
    {
        //set load factor to zero in order to get the same ordering of the displacements as before the routine
        fineScaleStructure->SetLoadFactor(0);
        fineScaleStructure->NodeBuildGlobalDofs();
        fineScaleStructure->NodeMergeActiveDofValues(activeDOF);
    }
    fineScaleStructure->GetLogger() << "Energy of fine scale exact  "  << fineScaleStructure->ElementTotalGetTotalEnergy()  << "\n";
    fineScaleStructure->GetLogger() << "Energy of fine scale approx "  << energy * (fineScaleStructure->GetAreaDamage()*fineScaleStructure->GetScalingFactorDamage()+
    		                                                  fineScaleStructure->GetAreaHomogeneous()*fineScaleStructure->GetScalingFactorHomogeneous())  << "\n";
    fineScaleStructure->GetLogger().CloseFile();

    return energy;
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

//! @brief ... set the elastic matrix
//! @param rElasticStiffness... elastic matrix
NuTo::FullMatrix<double> NuTo::Multiscale::GetElasticStiffness()const
{
	return mElasticStiffness;
}

//! @brief ... set the elastic matrix
//! @param rElasticStiffness... elastic matrix
void NuTo::Multiscale::SetElasticStiffness(NuTo::FullMatrix<double> rElasticStiffness)
{
    this->CheckElasticStiffness(rElasticStiffness);
	mElasticStiffness=rElasticStiffness;
    this->SetParametersValid();
}

//! @brief ... check elastic stiffness
//! @param rElasticStiffness ... crack transition radius
void NuTo::Multiscale::CheckElasticStiffness(const NuTo::FullMatrix<double>& rElasticStiffness) const
{
	NuTo::FullMatrix<double> eigenValues;
	if (rElasticStiffness.GetNumRows()!=3 || rElasticStiffness.GetNumColumns()!=3)
		throw NuTo::MechanicsException("[NuTo::Multiscale::CheckElasticStiffness] elastic stiffness not yet set, needs to be a 3x3 matrix.");
	rElasticStiffness.EigenValuesSymmetric(eigenValues);
    for (int count=0; count<eigenValues.GetNumRows(); count++)
    {
		if (eigenValues(count,0) <= 0.0)
		{
			throw NuTo::MechanicsException("[NuTo::Multiscale::CheckElasticStiffness] The eigenvalues of the elastic stiffness must be positive.");
		}
    }
}

//! @brief ... return the binary file from which the fine scale model is eventually deserialized
//! @return name of the file
std::string NuTo::Multiscale::GetMultiscaleFile()const
{
	return mFileName;
}

//! @brief ... set the binary file from which the fine scale model is eventually deserialized
//! @param rFileName... name of the file
void NuTo::Multiscale::SetMultiscaleFile(std::string rFileName)
{
	mFileName = rFileName;
}

//! @brief ... return crack transition radius to smooth the Heaviside function in the multiscale model
//! @return crack transition radius
double NuTo::Multiscale::GetCrackTransitionRadius()const
{
	return mCrackTransitionRadius;
}

//! @brief ... crack transition radius to smooth the Heaviside function in the multiscale model
//! @param rCrackTransitionRadius crack transition radius
void NuTo::Multiscale::SetCrackTransitionRadius(double rCrackTransitionRadius)
{
    this->CheckCrackTransitionRadius(rCrackTransitionRadius);
	mCrackTransitionRadius = rCrackTransitionRadius;
    this->SetParametersValid();
}

//! @brief ... check if crack transition radius is positive
//! @param rCrackTransitionRadius ... crack transition radius
void NuTo::Multiscale::CheckCrackTransitionRadius(double rCrackTransitionRadius) const
{
    if (rCrackTransitionRadius <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::Multiscale::CheckCrackTransitionRadius] The crack transition radius must be a positive value.");
    }
}

//! @brief ... get penalty stiffness crack angle
//! @return ... penalty stiffness crack angle
double NuTo::Multiscale::GetPenaltyStiffnessCrackAngle() const
{
	return mPenaltyStiffnessCrackAngle;
}

//! @brief ... set PenaltyStiffnessCrackAngle
//! @param rPenaltyStiffnessCrackAngle...  penalty stiffness crack angle
void NuTo::Multiscale::SetPenaltyStiffnessCrackAngle(double rPenaltyStiffnessCrackAngle)
{
    this->CheckPenaltyStiffnessCrackAngle(rPenaltyStiffnessCrackAngle);
	mPenaltyStiffnessCrackAngle = rPenaltyStiffnessCrackAngle;
	std::cout << "set penalty stiffness crack angle to " << mPenaltyStiffnessCrackAngle << "\n";
    this->SetParametersValid();
}

//! @brief ... check if penalty stiffness crack angle is positive
//! @param rPenaltyStiffnessCrackAngle ... PenaltyStiffnessCrackAngle
void NuTo::Multiscale::CheckPenaltyStiffnessCrackAngle(double rPenaltyStiffnessCrackAngle) const
{
    if (rPenaltyStiffnessCrackAngle <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::Multiscale::CheckPenaltyStiffnessCrackAngle] The penalty stiffness for the crack angle must be a positive value.");
    }
}

//! @brief ... get tensile strength
//! @return ... tensile strength
double NuTo::Multiscale::GetTensileStrength() const
{
    return mTensileStrength;
}

//! @brief ... set tensile strength
//! @param rTensileStrength...  tensile strength
void NuTo::Multiscale::SetTensileStrength(double rTensileStrength)
{
    this->CheckTensileStrength(rTensileStrength);
    this->mTensileStrength = rTensileStrength;
    this->SetParametersValid();
}

//! @brief ... check if tensile strength is positive
//! @param rTensileStrength ... tensile strength
void NuTo::Multiscale::CheckTensileStrength(double rTensileStrength) const
{
    if (rTensileStrength <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::Multiscale::CheckTensileStrength] The tensile strength must be a positive value.");
    }
}

//! @brief ... get scaling factor for the dofs of the crack angle
//! @return ... scaling factor
double NuTo::Multiscale::GetScalingFactorCrackAngle() const
{
    return mScalingFactorCrackAngle;
}

//! @brief ... set scaling factor for the dofs of the crack angle
//! @param rScalingFactor...  scaling factor
void NuTo::Multiscale::SetScalingFactorCrackAngle(double rScalingFactorCrackAngle)
{
    this->CheckScalingFactorCrackAngle(rScalingFactorCrackAngle);
    this->mScalingFactorCrackAngle = rScalingFactorCrackAngle;
    this->SetParametersValid();
}

//! @brief ... check if ScalingFactorCrackAngle is positive
//! @param rScalingFactorCrackAngle ... ScalingFactorCrackAngle
void NuTo::Multiscale::CheckScalingFactorCrackAngle(double rScalingFactorCrackAngle) const
{
    if (rScalingFactorCrackAngle <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::Multiscale::CheckScalingFactorCrackAngle] The scaling factor for the crack angle must be a positive value.");
    }
}

//! @brief ... get scaling factor for the dofs of the crack opening
//! @return ... scaling factor
double NuTo::Multiscale::GetScalingFactorCrackOpening() const
{
    return mScalingFactorCrackOpening;
}

//! @brief ... set scaling factor for the dofs of the crack opening
//! @param rScalingFactor...  scaling factor
void NuTo::Multiscale::SetScalingFactorCrackOpening(double rScalingFactorCrackOpening)
{
    this->CheckScalingFactorCrackOpening(rScalingFactorCrackOpening);
    this->mScalingFactorCrackOpening = rScalingFactorCrackOpening;
    this->SetParametersValid();
}

//! @brief ... check if ScalingFactorCrackOpening is positive
//! @param rScalingFactorCrackOpening ... ScalingFactorCrackOpening
void NuTo::Multiscale::CheckScalingFactorCrackOpening(double rScalingFactorCrackOpening) const
{
    if (rScalingFactorCrackOpening <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::Multiscale::CheckScalingFactorCrackOpening] The scaling factor for the crack opening must be a positive value.");
    }
}

//! @brief ... get scaling factor for the dofs of total strain
//! @return ... scaling factor
double NuTo::Multiscale::GetScalingFactorEpsilon() const
{
    return mScalingFactorEpsilon;
}

//! @brief ... set scaling factor for the dofs of the total strain
//! @param rScalingFactor...  scaling factor
void NuTo::Multiscale::SetScalingFactorEpsilon(double rScalingFactorEpsilon)
{
    this->CheckScalingFactorEpsilon(rScalingFactorEpsilon);
    this->mScalingFactorEpsilon = rScalingFactorEpsilon;
    this->SetParametersValid();
}

//! @brief ... check if ScalingFactorEpsilon is positive
//! @param rScalingFactorEpsilon ... ScalingFactorEpsilon
void NuTo::Multiscale::CheckScalingFactorEpsilon(double rScalingFactorEpsilon) const
{
    if (rScalingFactorEpsilon <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::Multiscale::CheckScalingFactorEpsilon] The scaling factor for the strain must be a positive value.");
    }
}

//! @brief ... get AugmentedLagrangeStiffnessCrackOpening
//! @return ...AugmentedLagrangeStiffnessCrackOpening
double NuTo::Multiscale::GetAugmentedLagrangeStiffnessCrackOpening() const
{
    return mAugmentedLagrangeStiffnessCrackOpening;
}

//! @brief ... set AugmentedLagrangeStiffnessCrackOpening
//! @param rAugmentedLagrangeStiffnessCrackOpening...AugmentedLagrangeStiffnessCrackOpening
void NuTo::Multiscale::SetAugmentedLagrangeStiffnessCrackOpening(double rAugmentedLagrangeStiffnessCrackOpening)
{
    this->CheckAugmentedLagrangeStiffnessCrackOpening(rAugmentedLagrangeStiffnessCrackOpening);
    this->mAugmentedLagrangeStiffnessCrackOpening = rAugmentedLagrangeStiffnessCrackOpening;
    this->SetParametersValid();

}

//! @brief ... check AugmentedLagrangeStiffnessCrackOpening
//! @param rAugmentedLagrangeStiffnessCrackOpening ...AugmentedLagrangeStiffnessCrackOpening
void NuTo::Multiscale::CheckAugmentedLagrangeStiffnessCrackOpening(double rAugmentedLagrangeStiffnessCrackOpening) const
{

}

//! @brief ... get ToleranceResidualForce
//! @return ...ToleranceResidualForce
double NuTo::Multiscale::GetToleranceResidualForce() const
{
    return mToleranceResidualForce;
}

//! @brief ... set ToleranceResidualForce
//! @param rToleranceResidualForce... ToleranceResidualForce
void NuTo::Multiscale::SetToleranceResidualForce(double rToleranceResidualForce)
{
    this->CheckToleranceResidualForce(rToleranceResidualForce);
    this->mToleranceResidualForce = rToleranceResidualForce;
    this->SetParametersValid();

}

//! @brief ... check ToleranceResidualForce
//! @param r ...
void NuTo::Multiscale::CheckToleranceResidualForce(double rToleranceResidualForce) const
{

}

//! @brief ... get MaxNumNewtonIterations
//! @return ...MaxNumNewtonIterations
int NuTo::Multiscale::GetMaxNumNewtonIterations() const
{
    return mMaxNumNewtonIterations;
}

//! @brief ... MaxNumNewtonIterations
//! @param rMaxNumNewtonIterations...MaxNumNewtonIterations
void NuTo::Multiscale::SetMaxNumNewtonIterations(int rMaxNumNewtonIterations)
{
    this->CheckMaxNumNewtonIterations(rMaxNumNewtonIterations);
    this->mMaxNumNewtonIterations = rMaxNumNewtonIterations;
    this->SetParametersValid();

}

//! @brief ... check MaxNumNewtonIterations
//! @param rMaxNumNewtonIterations ...MaxNumNewtonIterations
void NuTo::Multiscale::CheckMaxNumNewtonIterations(int rMaxNumNewtonIterations) const
{

}

//! @brief ... get MaxDeltaLoadFactor
//! @return ...MaxDeltaLoadFactor
double NuTo::Multiscale::GetMaxDeltaLoadFactor() const
{
    return mMaxDeltaLoadFactor;
}

//! @brief ... set MaxDeltaLoadFactor
//! @param rMaxDeltaLoadFactor...MaxDeltaLoadFactor
void NuTo::Multiscale::SetMaxDeltaLoadFactor(double rMaxDeltaLoadFactor)
{
    this->CheckMaxDeltaLoadFactor(rMaxDeltaLoadFactor);
    this->mMaxDeltaLoadFactor = rMaxDeltaLoadFactor;
    this->SetParametersValid();

}

//! @brief ... check MaxDeltaLoadFactor
//! @param rMaxDeltaLoadFactor ...MaxDeltaLoadFactor
void NuTo::Multiscale::CheckMaxDeltaLoadFactor(double rMaxDeltaLoadFactor) const
{

}

//! @brief ... get DecreaseFactor
//! @return ...DecreaseFactor
double NuTo::Multiscale::GetDecreaseFactor() const
{
    return mDecreaseFactor;
}

//! @brief ... set DecreaseFactor
//! @param rDecreaseFactor...DecreaseFactor
void NuTo::Multiscale::SetDecreaseFactor(double rDecreaseFactor)
{
    this->CheckDecreaseFactor(rDecreaseFactor);
    this->mDecreaseFactor = rDecreaseFactor;
    this->SetParametersValid();

}

//! @brief ... check DecreaseFactor
//! @param rDecreaseFactor ...DecreaseFactor
void NuTo::Multiscale::CheckDecreaseFactor(double r) const
{
}

//! @brief ... get MinNumNewtonIterations
//! @return ...MinNumNewtonIterations
int NuTo::Multiscale::GetMinNumNewtonIterations() const
{
    return mMinNumNewtonIterations;
}

//! @brief ... set MinNumNewtonIterations
//! @param rMinNumNewtonIterations...MinNumNewtonIterations
void NuTo::Multiscale::SetMinNumNewtonIterations(double rMinNumNewtonIterations)
{
    this->CheckMinNumNewtonIterations(rMinNumNewtonIterations);
    this->mMinNumNewtonIterations = rMinNumNewtonIterations;
    this->SetParametersValid();

}

//! @brief ... check MinNumNewtonIterations
//! @param rMinNumNewtonIterations ...
void NuTo::Multiscale::CheckMinNumNewtonIterations(int rMinNumNewtonIterations) const
{

}

//! @brief ... get IncreaseFactor
//! @return ...IncreaseFactor
double NuTo::Multiscale::GetIncreaseFactor() const
{
    return mIncreaseFactor;
}

//! @brief ... set IncreaseFactor
//! @param rIncreaseFactor...
void NuTo::Multiscale::SetIncreaseFactor(double rIncreaseFactor)
{
    this->CheckIncreaseFactor(rIncreaseFactor);
    this->mIncreaseFactor = rIncreaseFactor;
    this->SetParametersValid();

}

//! @brief ... check IncreaseFactor
//! @param rIncreaseFactor ...IncreaseFactor
void NuTo::Multiscale::CheckIncreaseFactor(double rIncreaseFactor) const
{

}

//! @brief ... get MinLoadFactor
//! @return ...MinLoadFactor
double NuTo::Multiscale::GetMinLoadFactor() const
{
    return mMinLoadFactor;
}

//! @brief ... set MinLoadFactor
//! @param rMinLoadFactor...MinLoadFactor
void NuTo::Multiscale::SetMinLoadFactor(double rMinLoadFactor)
{
    this->CheckMinLoadFactor(rMinLoadFactor);
    this->mMinLoadFactor = rMinLoadFactor;
    this->SetParametersValid();

}

//! @brief ... check MinLoadFactor
//! @param rMinLoadFactor ...MinLoadFactor
void NuTo::Multiscale::CheckMinLoadFactor(double rMinLoadFactor) const
{

}

//! @brief ... get MinLineSearchFactor
//! @return ... MinLineSearchFactor
double NuTo::Multiscale::GetMinLineSearchFactor() const
{
    return mMinLineSearchFactor;
}

//! @brief ... set MinLineSearchFactor
//! @param rMinLineSearchFactorFactor...MinLineSearchFactor
void NuTo::Multiscale::SetMinLineSearchFactor(double rMinLineSearchFactor)
{
    this->CheckMinLineSearchFactor(rMinLineSearchFactor);
    this->mMinLineSearchFactor = rMinLineSearchFactor;
    this->SetParametersValid();

}

//! @brief ... check MinLineSearchFactor
//! @param rMinLineSearchFactorFactor ...MinLineSearchFactor
void NuTo::Multiscale::CheckMinLineSearchFactor(double rMinLineSearchFactor) const
{

}

//! @brief ... get result directory
//! @return ... MinLineSearchFactor
const std::string& NuTo::Multiscale::GetResultDirectory()const
{
	return mResultDirectory;
}

//! @brief ... set ResultDirectory
//! @param rResultDirectory...ResultDirectory
void NuTo::Multiscale::SetResultDirectory(const std::string& rResultDirectory)
{
	this->CheckResultDirectory(rResultDirectory);
	this->mResultDirectory = rResultDirectory;
	this->SetParametersValid();
}

//! @brief ... check ResultDirectory
//! @param ResultDirectory ...ResultDirectory
void NuTo::Multiscale::CheckResultDirectory(const std::string& rResultDirectory) const
{

}

//! @brief ... get LoadStepMacro
//! @return ...LoadStepMacro
int NuTo::Multiscale::GetLoadStepMacro() const
{
    return mLoadStepMacro;
}

//! @brief ... set MinNumNewtonIterations
//! @param rMinNumNewtonIterations...MinNumNewtonIterations
void NuTo::Multiscale::SetLoadStepMacro(int rLoadStepMacro)
{
    this->CheckLoadStepMacro(rLoadStepMacro);
    this->mLoadStepMacro = rLoadStepMacro;
    this->SetParametersValid();

}

//! @brief ... check LoadStepMacro
//! @param rLoadStepMacro ...
void NuTo::Multiscale::CheckLoadStepMacro(int rLoadStepMacro) const
{

}

// check parameters
void NuTo::Multiscale::CheckParameters()const
{
	this->CheckElasticStiffness(this->mElasticStiffness);
	this->CheckCrackTransitionRadius(this->mCrackTransitionRadius);
	this->CheckPenaltyStiffnessCrackAngle(this->mPenaltyStiffnessCrackAngle);
	this->CheckTensileStrength(this->mTensileStrength);
	this->CheckScalingFactorCrackAngle(this->mScalingFactorCrackAngle);;
	this->CheckScalingFactorCrackOpening(this->mScalingFactorCrackOpening);
	this->CheckScalingFactorEpsilon(this->mScalingFactorEpsilon);
	this->CheckAugmentedLagrangeStiffnessCrackOpening(mAugmentedLagrangeStiffnessCrackOpening);
	this->CheckMaxNumNewtonIterations(mMaxNumNewtonIterations);
	this->CheckToleranceResidualForce(mToleranceResidualForce);
	this->CheckMaxDeltaLoadFactor(mMaxDeltaLoadFactor);
	this->CheckDecreaseFactor(mDecreaseFactor);
	this->CheckMinNumNewtonIterations(mMinNumNewtonIterations);
	this->CheckIncreaseFactor(mIncreaseFactor);
	this->CheckMinLoadFactor(mMinLoadFactor);
	this->CheckResultDirectory(mResultDirectory);
	this->CheckLoadStepMacro(mLoadStepMacro);
}
//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::Multiscale::Info(unsigned short rVerboseLevel) const
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

bool NuTo::Multiscale::CheckStiffness(NuTo::StructureMultiscale* rFineScaleStructure)const
{
    std::cout << "test of stiffness still included " << std::endl;
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
    NuTo::FullMatrix<double> dispForceVector;

    //recalculate stiffness
    rFineScaleStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    rFineScaleStructure->ConstraintInfo(10);

    NuTo::FullMatrix<double> stiffnessMatrixCSRVector2Full(stiffnessMatrixCSRVector2);
    //std::cout<<"stiffness matrix" << std::endl;
    //stiffnessMatrixCSRVector2Full.Info(10,3);
    double interval(1e-9);
    NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
    NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
    NuTo::FullMatrix<double> stiffnessMatrixCSRVector2_CDF(stiffnessMatrixCSRVector2.GetNumRows(), stiffnessMatrixCSRVector2.GetNumColumns());
    NuTo::FullMatrix<double> intForceVector1, intForceVector2, intForceVectorCDF(stiffnessMatrixCSRVector2.GetNumRows(),1);
    double energy1,energy2;
    rFineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
    rFineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
    rFineScaleStructure->ElementTotalUpdateTmpStaticData();
    rFineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector1);
    std::cout << "check stiffness:: intForceVector1"<< std::endl;
    intForceVector1.Trans().Info(12,3);
    energy1 = rFineScaleStructure->ElementTotalGetTotalEnergy();
    energy1 += rFineScaleStructure->ConstraintTotalGetTotalEnergy();
    for (int count=0; count<displacementsActiveDOFsCheck.GetNumRows(); count++)
    {
        displacementsActiveDOFsCheck(count,0)+=interval;
        rFineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        rFineScaleStructure->ElementTotalUpdateTmpStaticData();
        rFineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector2);
        std::cout << "check stiffness:: intForceVector2"<< std::endl;
        intForceVector2.Trans().Info(12,3);
        rFineScaleStructure->ConstraintInfo(10);
        energy2 = rFineScaleStructure->ElementTotalGetTotalEnergy();
        energy2 += rFineScaleStructure->ConstraintTotalGetTotalEnergy();
        stiffnessMatrixCSRVector2_CDF.SetColumn(count,(intForceVector2-intForceVector1)*(1./interval));
        intForceVectorCDF(count,0) = (energy2-energy1)/interval;
        displacementsActiveDOFsCheck(count,0)-=interval;
    }
    rFineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
    if ((stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).Abs().Max()>1e-3)
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
        //throw MechanicsException("[NuTo::Multiscale::Solve] Stiffness matrix is not correct.");
        std::cout << "stiffness is wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<< std::endl;
        return false;
    }
    else
    {
        std::cout << "stiffness is OK "<< std::endl;
        return true;
    }

}


bool NuTo::Multiscale::CheckGradient(NuTo::StructureMultiscale* rFineScaleStructure)const
{
    std::cout << "test of gradient still included " << std::endl;

    NuTo::FullMatrix<double> intForceVector;
    rFineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector);
     //std::cout<<"stiffness matrix" << std::endl;
    //stiffnessMatrixCSRVector2Full.Info(10,3);
    double interval(1e-9);
    NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
    NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
    NuTo::FullMatrix<double> intForceVectorCDF(intForceVector.GetNumRows(),1);
    double energy1,energy2;
    rFineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
    rFineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
    rFineScaleStructure->ElementTotalUpdateTmpStaticData();
    energy1 = rFineScaleStructure->ElementTotalGetTotalEnergy();
    energy1 += rFineScaleStructure->ConstraintTotalGetTotalEnergy();
    for (int count=0; count<displacementsActiveDOFsCheck.GetNumRows(); count++)
    {
        displacementsActiveDOFsCheck(count,0)+=interval;
        rFineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        rFineScaleStructure->ElementTotalUpdateTmpStaticData();
        energy2 = rFineScaleStructure->ElementTotalGetTotalEnergy();
        energy2 += rFineScaleStructure->ConstraintTotalGetTotalEnergy();
        intForceVectorCDF(count,0) = (energy2-energy1)/interval;
        displacementsActiveDOFsCheck(count,0)-=interval;
    }
    rFineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
    if ((intForceVectorCDF-intForceVector).Abs().Max()>1e-3)
    {
        std::cout << "intForceVector algo" << std::endl;
        intForceVector.Trans().Info(10,3);
        std::cout<< std::endl << "intForceVector cdf" << std::endl;
        intForceVectorCDF.Trans().Info(10,3);
        std::cout<< std::endl << "error" << std::endl;
        (intForceVectorCDF-intForceVector).Trans().Info(10);
        std::cout << "maximum error is " << (intForceVectorCDF-intForceVector).Abs().Max() << std::endl;
        //throw MechanicsException("[NuTo::Multiscale::Solve] Stiffness matrix is not correct.");
        std::cout << "intforce is wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<< std::endl;
        return false;
    }
    else
    {
        std::cout << "intforce is OK "<< std::endl;
        return true;
    }

}

void NuTo::Multiscale::MultiscaleSwitchToNonlinear(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient)const
{
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
   		//throw an exception giving information related to the wrong parameter
    	CheckParameters();
    	//if there is no exception thrown there is a problem with the source code
    	//since every time a material parameter is changed, the parametes should be checked
    	throw MechanicsException("[NuTo::Multiscale::MultiscaleSwitchToNonlinear] Check the material parameters.");
    }
    // this is a somehow weird situation, since for a normal material law nothing should be changed
    // since the material law is a full structure whose bc change, this can either be implemented with a cast (as I did)
    // or by removing the const flag from all material routines (which I do not consider is good)
    ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();

    if (staticData->NonlinearSolutionOn()==true)
    	return;

    //calculate inital elastic stress
    EngineeringStress2D stressElastic;
    GetEngineeringStressFromEngineeringStrain(rElement, rIp, rDeformationGradient, stressElastic);

    //check the principal stress
    double princ_sigma = 0.5*(stressElastic.mEngineeringStress[0]+stressElastic.mEngineeringStress[1])+
    		sqrt(0.25*(stressElastic.mEngineeringStress[0]-stressElastic.mEngineeringStress[1])*(stressElastic.mEngineeringStress[0]-stressElastic.mEngineeringStress[1])+
    				stressElastic.mEngineeringStress[2]*stressElastic.mEngineeringStress[2]);

    // check if transformation has to be done
    if (princ_sigma<mTensileStrength)
    	return;

    //calculate macro length
	double macroLength = sqrt(rElement->CalculateArea());

    //calculate center
	double center[3];
	rElement->GetGlobalIntegrationPointCoordinates(rIp, center);

	//IPName
	std::stringstream elementString;
	elementString << rElement->ElementGetId();
	std::stringstream rIPString;
	rIPString << rIp;
	std::string iPName = elementString.str()+std::string("_")+rIPString.str();

	rElement->GetStructure()->GetLogger() << "transform ip " << iPName << " to nonlinear structure " << "\n";
	staticData->SetFineScaleModel(mFileName, macroLength, center, iPName);
    StructureMultiscale *fineScaleStructure = const_cast<StructureMultiscale*>(staticData->GetFineScaleStructure());
    fineScaleStructure->LoggerSetQuiet(true);
    fineScaleStructure->GetLogger().OpenFile();

    fineScaleStructure->SetResultDirectory(mResultDirectory);
    std::stringstream ssLoadStep;
    ssLoadStep << mLoadStepMacro;
    fineScaleStructure->SetResultLoadStepMacro(ssLoadStep.str());

	const_cast<StructureBase*>(rElement->GetStructure())->SetUpdateTmpStaticDataRequired();
	//set the total strain constraint and recalculate the elastic crack angle
    fineScaleStructure->SetTotalEngineeringStrain(staticData->GetPrevStrain());
	//initialize the crack angle and the previous crack angle to the current elastic solution
	double alpha = fineScaleStructure->GetCrackAngleElastic();

	alpha = 0.5*M_PI;
	fineScaleStructure->SetCrackAngle(alpha);
	staticData->SetPrevCrackAngle(alpha);
	staticData->SetPrevCrackAngleElastic(alpha);

#ifdef SHOW_TIME
	fineScaleStructure->SetShowTime(false);
#endif

	//set constraint for negative crack opening
	fineScaleStructure->CreateConstraintLagrangeCrackOpening(mAugmentedLagrangeStiffnessCrackOpening);

	//set constraint for crack angle
	//bool coupleToTotalStrain(false);
	//fineScaleStructure->CreateConstraintNonlinearCrackAngle(mPenaltyStiffnessCrackAngle,coupleToTotalStrain);
	fineScaleStructure->CreateConstraintLinearGlobalCrackAngle(alpha);

	//set constraint for crack opening
	NuTo::FullMatrix<double> direction(2,1);
	direction(0,0)=0;
	direction(1,0)=1;
	fineScaleStructure->CreateConstraintLinearGlobalCrackOpening(0,direction);
	direction(0,0)=1;
	direction(1,0)=0;
	fineScaleStructure->CreateConstraintLinearGlobalCrackOpening(0,direction);

	//set crack transition zone
	fineScaleStructure->SetCrackTransitionRadius(mCrackTransitionRadius);

	//set scaling factors
	fineScaleStructure->SetScalingFactorCrackAngle(mScalingFactorCrackAngle);
	fineScaleStructure->SetScalingFactorCrackOpening(mScalingFactorCrackOpening);
	fineScaleStructure->SetScalingFactorEpsilon(mScalingFactorEpsilon);

	//calculate maximum independent sets
	fineScaleStructure->CalculateMaximumIndependentSets();

	//set to nonlinear solution
	staticData->SetNonlinearSolutionOn(true);
    fineScaleStructure->GetLogger().OpenFile();

    //change the previous angle of the alpha-constraint
    fineScaleStructure->SetPrevCrackAngle(alpha);
    fineScaleStructure->SetPrevCrackAngleElastic(alpha);

    //Get and set previous total strain
    EngineeringStrain2D prevStrain; //set to zero for the Newton Raphson iteration from zero to previous strain
    fineScaleStructure->SetPrevTotalEngineeringStrain(prevStrain);

    fineScaleStructure->SetLoadFactor(0);
    fineScaleStructure->NodeBuildGlobalDofs();
    NuTo::FullMatrix<double> activeDOF, dependentDOF;
    fineScaleStructure->NodeExtractDofValues(activeDOF,dependentDOF);
    fineScaleStructure->NodeMergeActiveDofValues(activeDOF);

    // calculate delta engineering strain from zero to the previous strain
    fineScaleStructure->SetDeltaTotalEngineeringStrain(staticData->GetPrevStrain());

    fineScaleStructure->SetNewtonRaphsonToleranceResidualForce(mToleranceResidualForce);
    fineScaleStructure->SetNewtonRaphsonAutomaticLoadStepControl(true);
    fineScaleStructure->SetNewtonRaphsonMaxDeltaLoadFactor(mMaxDeltaLoadFactor);
    fineScaleStructure->SetNewtonRaphsonMaxNumNewtonIterations(mMaxNumNewtonIterations);
    fineScaleStructure->SetNewtonRaphsonDecreaseFactor(mDecreaseFactor);
    fineScaleStructure->SetNewtonRaphsonMinNumNewtonIterations(mMinNumNewtonIterations);
    fineScaleStructure->SetNewtonRaphsonIncreaseFactor(mIncreaseFactor);
    fineScaleStructure->SetNewtonRaphsonMinDeltaLoadFactor(mMinLoadFactor);
    fineScaleStructure->SetNewtonRaphsonMinLineSearchFactor(mMinLineSearchFactor);

    std::stringstream saveStream;
    bool hasBeenSaved(false);

    fineScaleStructure->GetLogger() << "\n" << "******************************************************" << "\n";
    fineScaleStructure->GetLogger() << " MultiscaleSwitchToNonlinear after transformation" << "\n";
    fineScaleStructure->GetLogger() << " engineering strain " <<  staticData->GetPrevStrain().mEngineeringStrain[0] << " "
    		                            <<  staticData->GetPrevStrain().mEngineeringStrain[1] << " "
    		                            <<  staticData->GetPrevStrain().mEngineeringStrain[2]  << "\n";
    fineScaleStructure->GetLogger() << " deltaStrain strain " <<  staticData->GetPrevStrain().mEngineeringStrain[0] << " "
    		                            <<  staticData->GetPrevStrain().mEngineeringStrain[1] << " "
    		                            <<  staticData->GetPrevStrain().mEngineeringStrain[2] <<  "\n";
    fineScaleStructure->GetLogger() << " prevStrain strain "  <<  0 << " " <<  0 << " " <<  0  << "\n"  << "\n";
    try
    {
        fineScaleStructure->NewtonRaphson(false, saveStream, hasBeenSaved);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage(std::string("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] Error in performing Newton-iteration on fine scale after conversion from linear model for ip") + fineScaleStructure->GetIPName());
    	throw e;
    }
    catch(...)
    {
    	throw MechanicsException(std::string("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] Error in performing Newton-iteration on fine scale after conversion from linear model for ip") + fineScaleStructure->GetIPName());
    }
    fineScaleStructure->GetLogger().CloseFile();
    //calculate new inelastic stress starting from the previous iteration
    EngineeringStress2D stressInElastic;
    GetEngineeringStressFromEngineeringStrain(rElement, rIp, rDeformationGradient, stressInElastic);

    //check the difference between elastic and inelastic solution (fine scale and homogenized stress)
    fineScaleStructure->GetLogger().OpenFile();
    fineScaleStructure->GetLogger() << "elastic solution: " << stressElastic.mEngineeringStress[0] << " " << stressElastic.mEngineeringStress[1] << " " <<  stressElastic.mEngineeringStress[2] << "\n";
    fineScaleStructure->GetLogger() << "inelastic solution: " << stressInElastic.mEngineeringStress[0] << " " << stressInElastic.mEngineeringStress[1] << " " <<  stressInElastic.mEngineeringStress[2] << "\n";
    fineScaleStructure->GetLogger().CloseFile();
    if (fabs(stressElastic.mEngineeringStress[0]-stressInElastic.mEngineeringStress[0])+fabs(stressElastic.mEngineeringStress[1]-stressInElastic.mEngineeringStress[1])+fabs(stressElastic.mEngineeringStress[2]-stressInElastic.mEngineeringStress[2])>
        0.1*fabs(stressElastic.mEngineeringStress[0]+stressInElastic.mEngineeringStress[0])+fabs(stressElastic.mEngineeringStress[1]+stressInElastic.mEngineeringStress[1])+fabs(stressElastic.mEngineeringStress[2]+stressInElastic.mEngineeringStress[2]))
    {
    	throw MechanicsException("[NuTo::Multiscale::MultiscaleSwitchToNonlinear] the difference between elastic and inelastic solution is too big.");
    }
}
