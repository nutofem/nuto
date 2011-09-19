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
#include "nuto/mechanics/elements/Plane.h"
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
    mToleranceResidualForce=1e-6;
    mMaxDeltaLoadFactor = 1;
    mMaxNumNewtonIterations=20;
    mDecreaseFactor=0.5;
    mMinNumNewtonIterations=7;
    mIncreaseFactor=1.5;
    mMinLoadFactor=1e-15;
    mMinLineSearchFactor=1e-3;
    mAugmentedLagrangeStiffnessCrackOpening = 1e3;
    mCrackTransitionRadius = 0;

    mTensileStrength = 0;

    mDamageTresholdCrackInitiation = 1e-4;

    mScalingFactorCrackAngle = 0;
    mScalingFactorCrackOpening = 0;
    mScalingFactorEpsilon = 0;

    mUseAdditionalPeriodicShapeFunctions = false;
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
          & BOOST_SERIALIZATION_NVP(mLoadStepMacro)
          & BOOST_SERIALIZATION_NVP(mDamageTresholdCrackInitiation)
          & BOOST_SERIALIZATION_NVP(mUseAdditionalPeriodicShapeFunctions);
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
NuTo::Error::eError NuTo::Multiscale::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
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
NuTo::Error::eError NuTo::Multiscale::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
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
NuTo::Error::eError NuTo::Multiscale::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
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
NuTo::Error::eError NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
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
NuTo::Error::eError NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
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
NuTo::Error::eError NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
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
    if (staticData->LinearSolution())
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
        //fineScaleStructure->SetPrevCrackAngleElastic(staticData->GetPrevCrackAngleElastic());
        //fineScaleStructure->SetPrevCrackAngle(staticData->GetPrevCrackAngle());

        //Get and set previous total strain
        EngineeringStrain2D prevStrain(staticData->GetPrevStrain());
        fineScaleStructure->SetPrevTotalEngineeringStrain(prevStrain);

        fineScaleStructure->SetLoadFactor(0);
        fineScaleStructure->NodeBuildGlobalDofs();
        fineScaleStructure->NodeExtractDofValues(activeDOF,dependentDOF);

        fineScaleStructure->NodeMergeActiveDofValues(activeDOF);
        //std::cout << "last converged dofs " << activeDOF.Norm() <<  " " << dependentDOF.Norm()<< "\n";
        //activeDOF.Trans().Info(12,5);
        //dependentDOF.Trans().Info(12,5);

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
            //this might happen due to the adaptation
            bool initialStateInEquilibrium=false;
            Error::eError error = fineScaleStructure->NewtonRaphson(true, saveStream, hasBeenSaved, initialStateInEquilibrium);

            if (error==Error::NO_CONVERGENCE)
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
                std::cout << "[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] return with no convergence for ip "<< fineScaleStructure->GetIPName() << "\n";
                return error;
            }
            else
            {
				if (error!=Error::SUCCESSFUL)
				{
					throw MechanicsException(std::string("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] Newton iteration returned with flag other than successful"));
				}
            }
        }
        catch(NuTo::MechanicsException& e)
        {
        	e.AddMessage(std::string("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] MechanicsException while performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName());
        	throw e;
        }
        catch(...)
        {
        	throw MechanicsException(std::string("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName());
        }

        //calculate average stress
        NuTo::FullMatrix<double> averageStressDamage, averageStressHomogeneous;
        if (fineScaleStructure->GetGroupElementsDamage()==-1)
        {
			fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsHomogeneous(),fineScaleStructure->GetAreaFineScale(), averageStressHomogeneous);
	        rEngineeringStress.mEngineeringStress[0] = averageStressHomogeneous(0,0);
	        rEngineeringStress.mEngineeringStress[1] = averageStressHomogeneous(1,0);
	        rEngineeringStress.mEngineeringStress[2] = averageStressHomogeneous(3,0);
        }
        else
        {
			fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsDamage(),fineScaleStructure->GetAreaFineScale(), averageStressDamage);
			fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsHomogeneous(),fineScaleStructure->GetAreaFineScale(), averageStressHomogeneous);
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
        }

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

    return Error::SUCCESSFUL;
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
NuTo::Error::eError NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
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
    NuTo::Error::eError error = GetEngineeringStressFromEngineeringStrain(rElement, rIp,rDeformationGradient, engineeringStress2D);
    //this is certainly not correct, since it is plane strain, but I do not care about the stress in thickness direction here
    rEngineeringStress.mEngineeringStress[0] = engineeringStress2D.mEngineeringStress[0];
    rEngineeringStress.mEngineeringStress[1] = engineeringStress2D.mEngineeringStress[1];
    rEngineeringStress.mEngineeringStress[2] = 0.;
    rEngineeringStress.mEngineeringStress[3] = engineeringStress2D.mEngineeringStress[2];
    rEngineeringStress.mEngineeringStress[4] = 0.;
    rEngineeringStress.mEngineeringStress[5] = 0.;

    return error;
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
NuTo::Error::eError NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
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
NuTo::Error::eError NuTo::Multiscale::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient1D& rDeformationGradient, double& rDamage) const
{
	rDamage=0.;
    return Error::SUCCESSFUL;
}

//  Damage /////////////////////////////////////
//! @brief ... calculate isotropic damage from deformation gradient in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDamage ... damage variable
NuTo::Error::eError NuTo::Multiscale::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient2D& rDeformationGradient, double& rDamage) const
{
	rDamage=0.;
    return Error::SUCCESSFUL;
}

//  Damage /////////////////////////////////////
//! @brief ... calculate isotropic damage from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDamage ... damage variable
NuTo::Error::eError NuTo::Multiscale::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient3D& rDeformationGradient, double& rDamage) const
{
	rDamage=0.;
    return Error::SUCCESSFUL;
}

//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
NuTo::Error::eError NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
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
NuTo::Error::eError NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
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
    if (staticData->LinearSolution())
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
        //fineScaleStructure->SetPrevCrackAngleElastic(staticData->GetPrevCrackAngleElastic());
        //fineScaleStructure->SetPrevCrackAngle(staticData->GetPrevCrackAngle());

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
            //this might happen due to the adaptation
            bool initialStateInEquilibrium=false;
            Error::eError error = fineScaleStructure->NewtonRaphson(true, saveStream, hasBeenSaved, initialStateInEquilibrium);

            if (error==Error::NO_CONVERGENCE)
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
                std::cout << "[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] return with no convergence for ip "<< fineScaleStructure->GetIPName() << "\n";
                return error;
            }
            if (error!=Error::SUCCESSFUL)
            {
            	throw MechanicsException(std::string("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] Newton iteration returned with flag other than successful"));
            }
        }
        catch(NuTo::MechanicsException& e)
        {
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
        if (rhsVector.Abs().Max()>mToleranceResidualForce)
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
        double area = fineScaleStructure->GetCoarseScaleArea();
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

    return Error::SUCCESSFUL;
}

//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
NuTo::Error::eError NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
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
NuTo::Error::eError NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient)const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
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
    if (staticData->LinearSolution()==true)
    {
        // calculate Engineering stress
        engineeringStress.mEngineeringStress[0] = mElasticStiffness(0,0) * engineeringStrain.mEngineeringStrain[0] + mElasticStiffness(0,1) * engineeringStrain.mEngineeringStrain[1] + mElasticStiffness(0,2) * engineeringStrain.mEngineeringStrain[2];
        engineeringStress.mEngineeringStress[1] = mElasticStiffness(1,0) * engineeringStrain.mEngineeringStrain[0] + mElasticStiffness(1,1) * engineeringStrain.mEngineeringStrain[1] + mElasticStiffness(1,2) * engineeringStrain.mEngineeringStrain[2];
        engineeringStress.mEngineeringStress[2] = mElasticStiffness(2,0) * engineeringStrain.mEngineeringStrain[0] + mElasticStiffness(2,1) * engineeringStrain.mEngineeringStrain[1] + mElasticStiffness(2,2) * engineeringStrain.mEngineeringStrain[2];
    }
    else
    {
        StructureMultiscale *fineScaleStructure = const_cast<StructureMultiscale*>(staticData->GetFineScaleStructure());
        fineScaleStructure->GetLogger().OpenFile();

        //change the previous angle of the alpha-constraint and the previous crack angle
        //fineScaleStructure->SetPrevCrackAngleElastic(staticData->GetPrevCrackAngleElastic());
        //fineScaleStructure->SetPrevCrackAngle(staticData->GetPrevCrackAngle());

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
            //this might happen due to the adaptation
            bool initialStateInEquilibrium=false;
            NuTo::Error::eError error = fineScaleStructure->NewtonRaphson(true, saveStream, hasBeenSaved, initialStateInEquilibrium);

            if (error==Error::NO_CONVERGENCE)
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
                std::cout << "[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] return with no convergence for ip "<< fineScaleStructure->GetIPName() << "\n";
                return error;
            }
            if (error!=Error::SUCCESSFUL)
            {
            	throw MechanicsException(std::string("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] Newton iteration returned with flag other than successful"));
            }
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
        if (fineScaleStructure->GetGroupElementsDamage()==-1)
        {
			fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsHomogeneous(),fineScaleStructure->GetAreaFineScale(), averageStressHomogeneous);
			engineeringStress.mEngineeringStress[0] = averageStressHomogeneous(0,0);
			engineeringStress.mEngineeringStress[1] = averageStressHomogeneous(1,0);
			engineeringStress.mEngineeringStress[2] = averageStressHomogeneous(3,0);
        }
        else
        {
			fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsDamage(),fineScaleStructure->GetAreaFineScale(), averageStressDamage);
			fineScaleStructure->GetLogger() << "average stress in damage domain"  << "\n";
			fineScaleStructure->GetLogger().Out(averageStressDamage.Trans(),12,4);
			fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsHomogeneous(),fineScaleStructure->GetAreaFineScale(), averageStressHomogeneous);
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
        }

        //fineScaleStructure->GetLogger() << "Energy of fine scale exact  "  << fineScaleStructure->ElementTotalGetTotalEnergy()  << "\n";
        //fineScaleStructure->GetLogger() << "Energy of fine scale approx "  << energy * (fineScaleStructure->GetAreaDamage()*fineScaleStructure->GetScalingFactorDamage()+
        //		                                                  fineScaleStructure->GetAreaHomogeneous()*fineScaleStructure->GetScalingFactorHomogeneous())  << "\n";
        //fineScaleStructure->GetLogger() << "total area of macroscale " <<  (fineScaleStructure->GetAreaDamage()*fineScaleStructure->GetScalingFactorDamage()+
        //        fineScaleStructure->GetAreaHomogeneous()*fineScaleStructure->GetScalingFactorHomogeneous())  << "\n";
        //std::cout << "actual crack angle " << fineScaleStructure->GetCrackAngle()*180./M_PI << "prev crack angle " << staticData->GetPrevCrackAngle()*180./M_PI << "\n";
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

    return Error::SUCCESSFUL;
}
//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient)const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient)const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient)const
{
    return Error::SUCCESSFUL;
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
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
NuTo::Error::eError NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, double& rEnergy) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, double& rEnergy) const
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
    //fineScaleStructure->SetPrevCrackAngleElastic(staticData->GetPrevCrackAngleElastic());
    //fineScaleStructure->SetPrevCrackAngle(staticData->GetPrevCrackAngle());

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
        //this might happen due to the adaptation
        bool initialStateInEquilibrium=false;
        NuTo::Error::eError error = fineScaleStructure->NewtonRaphson(true, saveStream, hasBeenSaved, initialStateInEquilibrium);

        if (error==Error::NO_CONVERGENCE)
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
            std::cout << "[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] return with no convergence " << "\n";
            return error;
        }
        if (error!=Error::SUCCESSFUL)
        {
        	throw MechanicsException(std::string("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Newton iteration returned with flag other than successful"));
        }
    }
    catch(MechanicsException& e)
    {
        e.AddMessage(std::string("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName());
    	throw e;
    }
    catch(...)
    {
    	throw MechanicsException(std::string("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName());
    }

    rEnergy = staticData->GetPrevTotalEnergy();
    //calculate delta total energy (sigma1+sigma2)/2*delta_strain
    //calculate average stress
    NuTo::FullMatrix<double> averageStressDamage, averageStressHomogeneous;
    EngineeringStress2D meanEngineeringStress;
    const EngineeringStress2D& prevStress(staticData->GetPrevStress());
    if (fineScaleStructure->GetGroupElementsDamage()==-1)
    {
        fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsHomogeneous(),fineScaleStructure->GetAreaFineScale(), averageStressHomogeneous);
        meanEngineeringStress.mEngineeringStress[0] = (averageStressHomogeneous(0,0)+prevStress.mEngineeringStress[0]);
        meanEngineeringStress.mEngineeringStress[1] = (averageStressHomogeneous(1,0)+prevStress.mEngineeringStress[1]);
        meanEngineeringStress.mEngineeringStress[2] = (averageStressHomogeneous(3,0)+prevStress.mEngineeringStress[2]);
    }
    else
    {
		fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsDamage(),fineScaleStructure->GetAreaFineScale(), averageStressDamage);
		fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsHomogeneous(),fineScaleStructure->GetAreaFineScale(), averageStressHomogeneous);
		double scalingFactorDamage = fineScaleStructure->GetScalingFactorDamage();
		double scalingFactorHomogeneous = fineScaleStructure->GetScalingFactorHomogeneous();
		double sum(scalingFactorDamage+scalingFactorHomogeneous);
		scalingFactorDamage/=sum;
		scalingFactorHomogeneous/=sum;
		averageStressDamage*=scalingFactorDamage;
		averageStressHomogeneous*=scalingFactorHomogeneous;
		meanEngineeringStress.mEngineeringStress[0] = (averageStressDamage(0,0)+averageStressHomogeneous(0,0)+prevStress.mEngineeringStress[0]);
		meanEngineeringStress.mEngineeringStress[1] = (averageStressDamage(1,0)+averageStressHomogeneous(1,0)+prevStress.mEngineeringStress[1]);
		meanEngineeringStress.mEngineeringStress[2] = (averageStressDamage(3,0)+averageStressHomogeneous(3,0)+prevStress.mEngineeringStress[2]);
    }

    rEnergy+=0.5*(meanEngineeringStress.mEngineeringStress[0]*(engineeringStrain.mEngineeringStrain[0]-prevStrain.mEngineeringStrain[0])+
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
    fineScaleStructure->GetLogger() << "Energy of fine scale approx "  << rEnergy * fineScaleStructure->GetCoarseScaleArea() << "\n";
    fineScaleStructure->GetLogger().CloseFile();

    return Error::SUCCESSFUL;
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, double& rEnergy) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, double& rEnergy) const
{
    throw MechanicsException("[NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, double& rEnergy) const
{
    throw MechanicsException("[NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
NuTo::Error::eError NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, double& rEnergy) const
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

//! @brief ... get if additional periodic shape functions are used
//! @return ... true (periodic) or false (fixed displacements)
bool NuTo::Multiscale::GetUseAdditionalPeriodicShapeFunctions()const
{
	return mUseAdditionalPeriodicShapeFunctions;

}

//! @brief ... set to use additional periodic shape functions
//! @param rUseAddPeriodicShapeFunctions...rUseAddPeriodicShapeFunctions
void NuTo::Multiscale::SetUseAdditionalPeriodicShapeFunctions(bool rUseAdditionalPeriodicShapeFunctions)
{
    this->CheckUseAdditionalPeriodicShapeFunctions(rUseAdditionalPeriodicShapeFunctions);
    this->mUseAdditionalPeriodicShapeFunctions = rUseAdditionalPeriodicShapeFunctions;
    this->SetParametersValid();

}

//! @brief ... check
//! @param rUseAddPeriodicShapeFunctions ...rUseAddPeriodicShapeFunctions
void NuTo::Multiscale::CheckUseAdditionalPeriodicShapeFunctions(bool rUseAddPeriodicShapeFunctions) const
{

}

//! @brief ... get threshold for crack initiation based on the maximum damage value within the structure
//! @return ... mDamageTresholdCrackInitiation
double NuTo::Multiscale::GetDamageTresholdCrackInitiation() const
{
    return mDamageTresholdCrackInitiation;
}

//! @brief ... set DamageTresholdCrackInitiation
//! @param rDamageTresholdCrackInitiation...DamageTresholdCrackInitiation
void NuTo::Multiscale::SetDamageTresholdCrackInitiation(double rDamageTresholdCrackInitiation)
{
	mDamageTresholdCrackInitiation = rDamageTresholdCrackInitiation;
}

//! @brief ... check DamageTresholdCrackInitiation
//! @param rDamageTresholdCrackInitiation ...DamageTresholdCrackInitiation
void NuTo::Multiscale::CheckDamageTresholdCrackInitiation(double rDamageTresholdCrackInitiation) const
{
	if (mDamageTresholdCrackInitiation<0. ||mDamageTresholdCrackInitiation>1.)
    {
        throw NuTo::MechanicsException("[NuTo::Multiscale::CheckDamageTresholdCrackInitiation] The treshold must be within the interval [0,1].");
	}
}



// check parameters
void NuTo::Multiscale::CheckParameters()const
{
	this->CheckElasticStiffness(this->mElasticStiffness);
	this->CheckCrackTransitionRadius(this->mCrackTransitionRadius);
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
	this->CheckUseAdditionalPeriodicShapeFunctions(mUseAdditionalPeriodicShapeFunctions);
	this->CheckDamageTresholdCrackInitiation(mDamageTresholdCrackInitiation);
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

NuTo::Error::eError NuTo::Multiscale::MultiscaleSwitchToNonlinear(ElementBase* rElement, int rIp,
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

    if (staticData->LinearSolution())
    {
    	//switch from linear elastic to nonlinear no crack

		//calculate inital elastic stress
		EngineeringStress2D curStressElastic;
		GetEngineeringStressFromEngineeringStrain(rElement, rIp, rDeformationGradient, curStressElastic);

		//check the principal stress
		double princ_sigma = 0.5*(curStressElastic.mEngineeringStress[0]+curStressElastic.mEngineeringStress[1])+
				sqrt(0.25*(curStressElastic.mEngineeringStress[0]-curStressElastic.mEngineeringStress[1])*(curStressElastic.mEngineeringStress[0]-curStressElastic.mEngineeringStress[1])+
						curStressElastic.mEngineeringStress[2]*curStressElastic.mEngineeringStress[2]);

		// check if transformation has to be done
		if (princ_sigma<mTensileStrength)
			return Error::SUCCESSFUL;

		//calculate center
		double center[3];
		rElement->GetGlobalIntegrationPointCoordinates(rIp, center);

		//IPName
		std::stringstream elementString;
		elementString << rElement->ElementGetId();
		std::stringstream rIPString;
		rIPString << rIp;
		std::string iPName = elementString.str()+std::string("_")+rIPString.str();

		double macroCrackLength = 0; //tmp value
		rElement->GetStructure()->GetLogger() << "transform ip " << iPName << " to nonlinear structure " << "\n";
		staticData->SetFineScaleModel(mFileName, macroCrackLength, center, iPName);
		StructureMultiscale *fineScaleStructure = const_cast<StructureMultiscale*>(staticData->GetFineScaleStructure());
		fineScaleStructure->LoggerSetQuiet(true);
		fineScaleStructure->GetLogger().OpenFile();

		fineScaleStructure->SetResultDirectory(mResultDirectory);
		std::stringstream ssLoadStep;
		ssLoadStep << mLoadStepMacro;
		fineScaleStructure->SetResultLoadStepMacro(ssLoadStep.str());

		//this is a not so good way to say, that the Newton Rapshon iteration should be continued, since an adaptation has been performed
		const_cast<StructureBase*>(rElement->GetStructure())->SetUpdateTmpStaticDataRequired();

		//calculate macro area
		double macroArea = rElement->AsPlane()->CalculateArea();
		fineScaleStructure->SetCoarseScaleArea(macroArea);

	#ifdef SHOW_TIME
		fineScaleStructure->SetShowTime(false);
	#endif

		//set constraint for crack opening
		fineScaleStructure->CreateConstraintLinearGlobalCrackOpeningTangential(0);
		fineScaleStructure->CreateConstraintLinearGlobalCrackOpeningNormal(0);

		//set constraint for fine scale fluctuations on the boundary
		bool homogeneousDomain(true);
		if (mUseAdditionalPeriodicShapeFunctions)
		    fineScaleStructure->CreateConstraintLinearFineScaleDisplacementsUsingAddShapeFunctions(homogeneousDomain);
		else
		{
		    fineScaleStructure->CreateConstraintLinearFineScaleDisplacements(0,0,homogeneousDomain);
		    fineScaleStructure->CreateConstraintLinearPeriodicBoundaryShapeFunctions(0,0);
		    fineScaleStructure->CreateConstraintLinearPeriodicBoundaryShapeFunctions(1,0);
		    fineScaleStructure->CreateConstraintLinearPeriodicBoundaryShapeFunctions(2,0);
		}

		//set crack transition zone
		fineScaleStructure->SetCrackTransitionRadius(mCrackTransitionRadius);

		//set scaling factors
		fineScaleStructure->SetScalingFactorCrackOpening(mScalingFactorCrackOpening);
		fineScaleStructure->SetScalingFactorEpsilon(mScalingFactorEpsilon);

		//calculate maximum independent sets
		fineScaleStructure->CalculateMaximumIndependentSets();

		//set to nonlinear solution
		staticData->SetSolutionPhase(Constitutive::NONLINEAR_NO_CRACK);
		//staticData->SetSolutionPhase(Constitutive::NONLINEAR_CRACKED);
		fineScaleStructure->GetLogger().OpenFile();

		//change the previous angle of the alpha-constraint
		//fineScaleStructure->SetPrevCrackAngle(alpha);
		//fineScaleStructure->SetPrevCrackAngleElastic(alpha);

		//Get and set previous total strain
		EngineeringStrain2D prevStrain; //set to zero for the Newton Raphson iteration from zero to previous strain
		fineScaleStructure->SetPrevTotalEngineeringStrain(prevStrain);

		fineScaleStructure->SetLoadFactor(0);
		//fineScaleStructure->ConstraintInfo(10);
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
            //this might happen due to the adaptation
            bool initialStateInEquilibrium=false;
            fineScaleStructure->NewtonRaphson(false, saveStream, hasBeenSaved, initialStateInEquilibrium);
	        fineScaleStructure->ElementTotalUpdateStaticData();
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
		NuTo::FullMatrix<double> stressInElastic;
		fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsHomogeneous(),fineScaleStructure->GetAreaFineScale(), stressInElastic);
		const EngineeringStress2D& stressElastic(staticData->GetPrevStress());

		//check the difference between elastic and inelastic solution (fine scale and homogenized stress)
		rElement->GetStructure()->GetLogger().OpenFile();
		rElement->GetStructure()->GetLogger() << "elastic solution: " << stressElastic.mEngineeringStress[0] << " " << stressElastic.mEngineeringStress[1] << " " <<  stressElastic.mEngineeringStress[2] << "\n";
		rElement->GetStructure()->GetLogger() << "inelastic solution: " << stressInElastic(0,0) << " " << stressInElastic(1,0) << " " <<  stressInElastic(3,0) << "\n";
		rElement->GetStructure()->GetLogger().CloseFile();
		if (fabs(stressElastic.mEngineeringStress[0]-stressInElastic(0,0))+fabs(stressElastic.mEngineeringStress[1]-stressInElastic(1,0))+fabs(stressElastic.mEngineeringStress[2]-stressInElastic(3,0))>
			0.1*fabs(stressElastic.mEngineeringStress[0]+stressInElastic(0,0))+fabs(stressElastic.mEngineeringStress[1]+stressInElastic(1,0))+fabs(stressElastic.mEngineeringStress[2]+stressInElastic(3,0)))
		{
			rElement->GetStructure()->GetLogger() << "********** reenable check for difference between elastic and nonlinear solution " << "\n";
			//throw MechanicsException("[NuTo::Multiscale::MultiscaleSwitchToNonlinear] the difference between elastic and inelastic solution is too big.",NuTo::MechanicsException::NOCONVERGENCE);
		}
    }
    else
    {
        if (staticData->NonLinearSolutionNoCrack())
        {
            //first calculate the actual state
        	//nonlinear solution
            ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();
            StructureMultiscale *fineScaleStructure = const_cast<StructureMultiscale*>(staticData->GetFineScaleStructure());
            fineScaleStructure->GetLogger().OpenFile();

            //change the previous angle of the alpha-constraint
            //fineScaleStructure->SetPrevCrackAngleElastic(staticData->GetPrevCrackAngleElastic());
            //fineScaleStructure->SetPrevCrackAngle(staticData->GetPrevCrackAngle());

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
                fineScaleStructure->GetLogger() << " Try to switch from nonlinear without crack to cracked solution" << "\n";
                fineScaleStructure->GetLogger() << " engineering strain " <<  engineeringStrain.mEngineeringStrain[0] << " "
                		                            <<  engineeringStrain.mEngineeringStrain[1] << " "
                		                            <<  engineeringStrain.mEngineeringStrain[2] << "\n";
                fineScaleStructure->GetLogger() << " deltaStrain strain " <<  deltaStrain.mEngineeringStrain[0] << " "
                		                            <<  deltaStrain.mEngineeringStrain[1] << " "
                		                            <<  deltaStrain.mEngineeringStrain[2] << "\n";
                fineScaleStructure->GetLogger() << " prevStrain strain "  <<  staticData->GetPrevStrain().mEngineeringStrain[0] << " "
                		                            <<  staticData->GetPrevStrain().mEngineeringStrain[1] << " "
                		                            <<  staticData->GetPrevStrain().mEngineeringStrain[2] << "\n"  << "\n";
                //std::clock_t start,end;
                //start=clock();
                //this might happen due to the adaptation
                bool initialStateInEquilibrium=false;
                fineScaleStructure->NewtonRaphson(true, saveStream, hasBeenSaved, initialStateInEquilibrium);
                //end=clock();
                //std::cout << "time for not enriched Newton-Raphson solution " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
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
                std::cout << std::string("[NuTo::Multiscale::MultiscaleSwitschToNonlinear] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName() << std::endl;
                e.AddMessage(std::string("[NuTo::Multiscale::MultiscaleSwitschToNonlinear] No convergence in multiscale for ip ") + fineScaleStructure->GetIPName());
            	throw e;
            }
            catch(...)
            {
            	throw MechanicsException(std::string("[NuTo::Multiscale::MultiscaleSwitschToNonlinear] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName());
            }


            //check the damage state
            double maxDamage = fineScaleStructure->ElementTotalGetMaxDamage();
			double maxShift[2] = { 0,0 };
			double maxAlpha = 0;

			if (maxDamage>=mDamageTresholdCrackInitiation)
            {
               //initial crack angle from the maximum principal strain
				double princAlpha = fineScaleStructure->CalculateCrackAnglePrincipalStrain(engineeringStrain);

				std::clock_t start,end;
				start=clock();
				int numCountShift=5;
				double rangeShiftNormal = 0.9*sqrt(fineScaleStructure->GetAreaFineScale());
				double initShiftNormal(-0.5*rangeShiftNormal);
				double deltaShiftNormal=rangeShiftNormal/(numCountShift-1);
				//double initShiftNormal(0);
				//double deltaShiftNormal = 0.;

				int numAlpha=5;
				double initAlpha = princAlpha-M_PI*0.25;
				double deltaAlpha = 0.5*M_PI/(numAlpha-1);
				//double initAlpha = M_PI*0.5;
				//double deltaAlpha = 0.;

				double shift[2];
				const boost::array<double,2>& centerDamage(fineScaleStructure->GetCenterDamage());

				//get actual coordinates and disp of all nodes
				NuTo::FullMatrix<double> coordinates;
				NuTo::FullMatrix<double> displacementsTotalReal;
				int nodeGroupHomogeneous = fineScaleStructure->GetGroupNodesHomogeneous();
				fineScaleStructure->NodeGroupGetCoordinates(nodeGroupHomogeneous,coordinates);
				//std::cout << "coordinates " << "\n";
				//coordinates.Info(12,3);
				fineScaleStructure->NodeGroupGetDisplacements(nodeGroupHomogeneous,displacementsTotalReal);
				//std::cout << "displacementsTotalReal " << "\n";
				//displacementsTotalReal.Info(12,6);
				double maxNormCrackOpeningSquare = -1;

				//calculate A
				int numNodes=coordinates.GetNumRows();
				int numDisp = displacementsTotalReal.GetNumColumns();
				NuTo::FullMatrix<double> A(numDisp*numNodes,1);
				for (int countNode=0; countNode<coordinates.GetNumRows(); countNode++)
				{
					A(countNode,0) = displacementsTotalReal(countNode,0);
					A(countNode+numNodes,0) = displacementsTotalReal(countNode,1);
				}

				for (int countAlpha = 0; countAlpha<numAlpha; countAlpha++)
				{
					double alpha = initAlpha+countAlpha*deltaAlpha;
					fineScaleStructure->SetCrackAngle(alpha);
					for (int countShift=0; countShift<numCountShift; countShift++)
					{
						double normalShift(initShiftNormal+deltaShiftNormal*countShift);
						fineScaleStructure->SetCrackAngle(alpha);
						shift[0] = normalShift*sin(alpha);
						shift[1] = -normalShift*cos(alpha);
						fineScaleStructure->SetShiftCenterDamage(shift);

						//calculate matrices according to A=B*([exx_hom, eyy_hom, gxy_hom, ut, un]^T)
						assert(numDisp==2);
						NuTo::FullMatrix<double> B(numDisp*numNodes,5);

						//set crack opening to zero
						NuTo::FullMatrix<double> crackOpening(2,1);
						crackOpening(0,0) = 0;
						crackOpening(1,0) = 0;
						fineScaleStructure->SetCrackOpening(crackOpening);

						NuTo::EngineeringStrain2D homStrain;
						for (int countEpsilon=0; countEpsilon<3; countEpsilon++)
						{
							//Set homogeneous strain
							homStrain.mEngineeringStrain[countEpsilon]=1;
							fineScaleStructure->SetHomogeneousEngineeringStrain(homStrain);

							for (int countNode=0; countNode<coordinates.GetNumRows(); countNode++)
							{
								double coordinatesNode[2];
								coordinatesNode[0] = coordinates(countNode,0);
								coordinatesNode[1] = coordinates(countNode,1);
								double displacementsNodeHom[2];
								fineScaleStructure->GetDisplacementsEpsilonHom2D(coordinatesNode, displacementsNodeHom, centerDamage);
								B(countNode,countEpsilon) = displacementsNodeHom[0];
								B(countNode+numNodes,countEpsilon) = displacementsNodeHom[1];
							}
							homStrain.mEngineeringStrain[countEpsilon]=0;
						}
						fineScaleStructure->SetHomogeneousEngineeringStrain(homStrain);
						//std::cout << "displacements hom " << "\n";
						//A.Info(12,6);

						//set tangential crackopening to 1
						crackOpening(0,0) = 1;
						crackOpening(1,0) = 0;
						fineScaleStructure->SetCrackOpening(crackOpening);
						for (int countNode=0; countNode<coordinates.GetNumRows(); countNode++)
						{
							double coordinatesNode[2];
							coordinatesNode[0] = coordinates(countNode,0);
							coordinatesNode[1] = coordinates(countNode,1);
							double displacementsNodeCrack[2];
							//the homogeneous displacements change as well, since the crack opening changes the hom strain due to constant total strain
							fineScaleStructure->GetDisplacementsCrack2D(coordinatesNode, displacementsNodeCrack);
							//fineScaleStructure->GetDisplacementsEpsilonHom2D(coordinatesNode, displacementsNodeHom, centerDamage);
							B(countNode,3) = displacementsNodeCrack[0];// + displacementsNodeHom[0]-A(countNode,0);
							B(countNode+numNodes,3) = displacementsNodeCrack[1];// + displacementsNodeHom[1]-A(countNode+numNodes,0);
						}

						//set normal crackopening to 1
						crackOpening(0,0) = 0;
						crackOpening(1,0) = 1;
						//dont forget to recalculate the homogenous part of the strain
						fineScaleStructure->SetCrackOpening(crackOpening);
						for (int countNode=0; countNode<coordinates.GetNumRows(); countNode++)
						{
							double coordinatesNode[2];
							coordinatesNode[0] = coordinates(countNode,0);
							coordinatesNode[1] = coordinates(countNode,1);
							double displacementsNodeCrack[2];
							//the homogeneous displacements change as well, since the crack opening changes the hom strain due to constant total strain
							fineScaleStructure->GetDisplacementsCrack2D(coordinatesNode, displacementsNodeCrack);
							//fineScaleStructure->GetDisplacementsEpsilonHom2D(coordinatesNode, displacementsNodeHom, centerDamage);
							B(countNode,4) = displacementsNodeCrack[0];// + displacementsNodeHom[0]-A(countNode,0);
							B(countNode+numNodes,4) = displacementsNodeCrack[1];// + displacementsNodeHom[1]-A(countNode+numNodes,0);
						}
						//std::cout << "displacements crack " << "\n";
						//B.Info(12,3);

						//solve the least squares problem to find un and ut
						NuTo::FullMatrix<double> sol;
						//std::cout << "A " << "\n";
						//A.Info(12,5);
						//std::cout << "B " << "\n";
						//B.Info(12,5);
						sol = (B.Trans()*B).Inverse()*B.Trans()*A;

						//no negativ crack opening in normal direction
						if (sol(4,0)<0)
							sol(4,0) = 0;

						//std::cout << "normalShift " << normalShift << " x " << shift[0] << " y " << shift[1] << " ut " << sol(0,0) << " un " << sol(1,0) << "\n";
						double normCrackOpeningSquare(sol(3,0)*sol(3,0)+sol(4,0)*sol(4,0));
						if (normCrackOpeningSquare>maxNormCrackOpeningSquare)
						{
							maxNormCrackOpeningSquare = normCrackOpeningSquare;
							maxShift[0] = shift[0];
							maxShift[1] = shift[1];
							maxAlpha = alpha;
						}
					}
				}
				fineScaleStructure->GetLogger() << "max crack opening for alpha=" << maxAlpha*180/M_PI << " and shift " <<  maxShift[0] << " " << maxShift[1] << "\n";
				end=clock();
				fineScaleStructure->GetLogger() << "time for solving least squares problem " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
            }
    	    //restore structure
    	    if (hasBeenSaved)
    	    {
    	        fineScaleStructure->RestoreStructure(saveStream);
    	    }
    	    else
    	    {
    			//set crack opening to zero
    			NuTo::FullMatrix<double> crackOpening(2,1);
    			crackOpening(0,0) = 0;
    			crackOpening(1,0) = 0;
    			fineScaleStructure->SetCrackOpening(crackOpening);
    	        //set load factor to zero in order to get the same ordering of the displacements as before the routine
    	        fineScaleStructure->SetLoadFactor(0);
    	        fineScaleStructure->NodeBuildGlobalDofs();
    	        fineScaleStructure->NodeMergeActiveDofValues(activeDOF);
    	    }

    	    //if (engineeringStrain.mEngineeringStrain[0]>0.0001)
	    	//if ((averageStressWithCrack-averageStressNoCrack).Norm()>0.01*2)
    		//if (sqrt(crackOpening[0]*crackOpening[0]+crackOpening[1]*crackOpening[1])>0.0001)
            if (maxDamage>=mDamageTresholdCrackInitiation)
            {
				//test the difference in average stress between the initial solution (no crack enrichment) and the enriched solution (crack enrichment at pos maxShift)
                rElement->GetStructure()->GetLogger() << "Add a crack with angle " << maxAlpha*180/M_PI << " with shift " << maxShift[0] << " " << maxShift[1] << "\n";
                //rElement->GetStructure()->GetLogger() << "Add a crack with angle " << maxAlpha*180/M_PI << " with shift " << maxShift[0] << " " << maxShift[1] << "\n";
				//std::cin.getline (title,256);

				//set the crack shift
				fineScaleStructure->SetShiftCenterDamage(maxShift);
				//set crack angle
				fineScaleStructure->SetCrackAngle(maxAlpha);

				//create the damage domain by coyping and shifting the homogeneous domain
				fineScaleStructure->CreateDamageDomainFromHomogeneousDomain();

				//calculate macro length
				double macroCrackLength = rElement->AsPlane()->CalculateCrackLength2D(maxAlpha);
				fineScaleStructure->SetlCoarseScaleCrack(macroCrackLength);

				//set the fine scale crack length
				fineScaleStructure->CalculateAndSetCrackLengthFineScale();

				if (fineScaleStructure->GetScalingFactorDamage()<=1e-10 || fineScaleStructure->GetScalingFactorHomogeneous()<=1e-10)
				{
					std::cout << "scaling factor hom " << fineScaleStructure->GetScalingFactorHomogeneous() << " damage " << fineScaleStructure->GetScalingFactorDamage() << "\n";
					std::cout << "mCoarseScaleArea " << fineScaleStructure->GetCoarseScaleArea() << "\n";
					std::cout << "mlCoarseScale " << fineScaleStructure->GetlCoarseScaleCrack() << "\n";
					std::cout << "mFineScaleArea " << fineScaleStructure->GetAreaFineScale() << "\n";
					std::cout << "mlFineScaleDamage " << fineScaleStructure->GetlFineScaleCrack() << "\n";
					throw MechanicsException("[NuTo::Multiscale::MultiscaleSwitchToNonlinear] scaling factor is less than 0, probably your macro element is smaller than the fine scale model." );
				}

				//delete constraints for crack opening
				fineScaleStructure->ConstraintDeleteTangentialCrackOpening();
				fineScaleStructure->ConstraintDeleteNormalCrackOpening();
				//delete linear constraint for normal crack opening and create constraint to avoid negative crack opening
				//std::cout << "Constraint Lagrange for crack opening removed " << "\n";
				fineScaleStructure->CreateConstraintLagrangeGlobalCrackOpeningNormal(mAugmentedLagrangeStiffnessCrackOpening);

				//set constraint for fine scale fluctuations on the boundary of the damaged domain
				if (mUseAdditionalPeriodicShapeFunctions)
					fineScaleStructure->CreateConstraintLinearFineScaleDisplacementsUsingAddShapeFunctions(false);

				else
				{
					fineScaleStructure->CreateConstraintLinearFineScaleDisplacements(0,0,false);
				}

				//set to nonlinear solution with crack
				staticData->SetSolutionPhase(Constitutive::NONLINEAR_CRACKED);
				//this is a not so good way to say, that the Newton Rapshon iteration should be continued, since an adaptation has been performed
				const_cast<StructureBase*>(rElement->GetStructure())->SetUpdateTmpStaticDataRequired();

				//test the difference in average stress between the initial solution (no crack enrichment) and the enriched solution (crack enrichment at pos maxShift)
				//Get and set previous total strain
				fineScaleStructure->SetPrevTotalEngineeringStrain(prevStrain);

				fineScaleStructure->SetLoadFactor(0);
				fineScaleStructure->NodeBuildGlobalDofs();
				fineScaleStructure->NodeExtractDofValues(activeDOF,dependentDOF);
				fineScaleStructure->NodeMergeActiveDofValues(activeDOF);

				//Get and set previous delta strain
				EngineeringStrain2D zeroDeltaStrain;
				fineScaleStructure->SetDeltaTotalEngineeringStrain(zeroDeltaStrain);

				hasBeenSaved = false;
				try
				{
					fineScaleStructure->GetLogger() << "\n" << "************************************************" << "\n";
					fineScaleStructure->GetLogger() << " Switch from nonlinear without crack to cracked solution with crack" << "\n";
					fineScaleStructure->GetLogger() << " engineering strain " <<  engineeringStrain.mEngineeringStrain[0] << " "
														<<  engineeringStrain.mEngineeringStrain[1] << " "
														<<  engineeringStrain.mEngineeringStrain[2] << "\n";
					fineScaleStructure->GetLogger() << " deltaStrain strain " <<  deltaStrain.mEngineeringStrain[0] << " "
														<<  deltaStrain.mEngineeringStrain[1] << " "
														<<  deltaStrain.mEngineeringStrain[2] << "\n";
					fineScaleStructure->GetLogger() << " prevStrain strain "  <<  staticData->GetPrevStrain().mEngineeringStrain[0] << " "
														<<  staticData->GetPrevStrain().mEngineeringStrain[1] << " "
														<<  staticData->GetPrevStrain().mEngineeringStrain[2] << "\n"  << "\n";
					//std::clock_t start,end;
					//start=clock();
					bool saveStructureBeforeUpdate(false);
					bool initialStateInEquilibrium(false);
					fineScaleStructure->NewtonRaphson(saveStructureBeforeUpdate, saveStream, hasBeenSaved, initialStateInEquilibrium);
			        fineScaleStructure->ElementTotalUpdateStaticData();
					//end=clock();
					//std::cout << "time for enriched Newton-Raphson solution " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
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
					std::cout << std::string("[NuTo::Multiscale::MultiscaleSwitschToNonlinear] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName() << std::endl;
					e.AddMessage(std::string("[NuTo::Multiscale::MultiscaleSwitschToNonlinear] No convergence in multiscale for ip ") + fineScaleStructure->GetIPName());
					throw e;
				}
				catch(...)
				{
					throw MechanicsException(std::string("[NuTo::Multiscale::MultiscaleSwitschToNonlinear] Error in performing Newton-iteration on fine scale for ip ") + fineScaleStructure->GetIPName());
				}

				//calculate average stress
				NuTo::FullMatrix<double> averageStressWithCrack, averageStressDamage, averageStressHomogeneous;
				fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsDamage(),fineScaleStructure->GetAreaFineScale(), averageStressDamage);
				fineScaleStructure->ElementGroupGetAverageStress(fineScaleStructure->GetGroupElementsHomogeneous(),fineScaleStructure->GetAreaFineScale(), averageStressHomogeneous);
				double scalingFactorDamage = fineScaleStructure->GetScalingFactorDamage();
				double scalingFactorHomogeneous = fineScaleStructure->GetScalingFactorHomogeneous();
				double sum = scalingFactorDamage + scalingFactorHomogeneous;
				scalingFactorDamage/=sum;
				scalingFactorHomogeneous/=sum;
				averageStressWithCrack=averageStressDamage*scalingFactorDamage+averageStressHomogeneous*scalingFactorHomogeneous;
				boost::array<double,2> crackOpening=fineScaleStructure->GetGlobalCrackOpening2D();

				const EngineeringStress2D& stressWithoutCrack(staticData->GetPrevStress());
				//check the difference between elastic and inelastic solution (fine scale and homogenized stress)
				rElement->GetStructure()->GetLogger().OpenFile();
				rElement->GetStructure()->GetLogger() << "average stress without crack: " << stressWithoutCrack.mEngineeringStress[0] << " " << stressWithoutCrack.mEngineeringStress[1] << " " <<  stressWithoutCrack.mEngineeringStress[2] << "\n";
				rElement->GetStructure()->GetLogger() << "average stress with crack: " << averageStressWithCrack(0,0) << " " << averageStressWithCrack(1,0) << " " <<  averageStressWithCrack(3,0) << "\n";
				rElement->GetStructure()->GetLogger() << "crackopening after adaptation: " << crackOpening[0] << "[t]," << crackOpening[1] << "[n] \n";
				rElement->GetStructure()->GetLogger().CloseFile();
				if (fabs(stressWithoutCrack.mEngineeringStress[0]-averageStressWithCrack(0,0))+fabs(stressWithoutCrack.mEngineeringStress[1]-averageStressWithCrack(1,0))+fabs(stressWithoutCrack.mEngineeringStress[2]-averageStressWithCrack(3,0))>
					0.1*fabs(stressWithoutCrack.mEngineeringStress[0]+averageStressWithCrack(0,0))+fabs(stressWithoutCrack.mEngineeringStress[1]+averageStressWithCrack(1,0))+fabs(stressWithoutCrack.mEngineeringStress[2]+averageStressWithCrack(3,0)))
				{
					rElement->GetStructure()->GetLogger() << "********** reenable check for difference between cracked and noncracked solution " << "\n";
					//throw MechanicsException("[NuTo::Multiscale::MultiscaleSwitchToNonlinear] the difference between elastic and inelastic solution is too big.",NuTo::MechanicsException::NOCONVERGENCE);
				}
            }
    	    else
    	    {
    			//std::cout << "continue without crack enrichment "<< (averageStressWithCrack-averageStressNoCrack).Norm() << "<"<< 0.01*2 << "\n";
    			//std::cin.getline (title,256);
    	    }
        }
    }
    return Error::NOT_IMPLEMENTED;
}
