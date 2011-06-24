// $Id: ConstraintGlobalCrackAngle.cpp 314 2010-09-27 16:31:43Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <iostream>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constraints/ConstraintNonlinearGlobalCrackAngle2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"

// constructor
NuTo::ConstraintNonlinearGlobalCrackAngle2D::ConstraintNonlinearGlobalCrackAngle2D(const StructureMultiscale* rStructure,
        double rPenaltyStiffness, bool rCoupleToTotalStrain):
        ConstraintNonlinear()
{
    mStructure = rStructure;
    mPenaltyStiffness = rPenaltyStiffness;
    mCoupleToTotalStrain = rCoupleToTotalStrain;
}

NuTo::ConstraintNonlinearGlobalCrackAngle2D* NuTo::ConstraintNonlinearGlobalCrackAngle2D::AsConstraintNonlinearGlobalCrackAngle2D()
{
    return this;
}

const NuTo::ConstraintNonlinearGlobalCrackAngle2D* NuTo::ConstraintNonlinearGlobalCrackAngle2D::AsConstraintNonlinearGlobalCrackAngle2D()const
{
    return this;
}

void NuTo::ConstraintNonlinearGlobalCrackAngle2D::SetPenaltyStiffness(double rPenaltyStiffness)
{
    mPenaltyStiffness = rPenaltyStiffness;
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::ConstraintNonlinearGlobalCrackAngle2D::Info(unsigned short rVerboseLevel) const
{
	mStructure->GetLogger() << "ConstraintNonLinearGlobalCrackAngle : delta angle " << mStructure->CalculateDeltaCrackAngleElastic() << "\n";
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix (Hessian) of the total potential
//! @param rResult ... coefficient matrix
//! @param rGlobalDofs ... row and column numbers in global system
void NuTo::ConstraintNonlinearGlobalCrackAngle2D::CalculateCoefficientMatrix_0(NuTo::SparseMatrixCSRVector2Symmetric<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
	if (mCoupleToTotalStrain)
	{
		//crack angle
		int dof(4);
		rResult.Resize(dof,dof);
		rGlobalDofs.resize(dof,1);

		rGlobalDofs[0] = mStructure->GetDofCrackAngle();
		rGlobalDofs[1] = mStructure->GetDofGlobalTotalStrain2D()[0];
		rGlobalDofs[2] = mStructure->GetDofGlobalTotalStrain2D()[1];
		rGlobalDofs[3] = mStructure->GetDofGlobalTotalStrain2D()[2];

		//std::cout << "global dofs " << rGlobalDofs[0] << " " << rGlobalDofs[1] << " " << rGlobalDofs[2] << " " << rGlobalDofs[3] << std::endl;

		//calculate angle orthogonal to second principal stress
		double delta_alpha = mStructure->CalculateDeltaCrackAngleElastic();

		NuTo::FullMatrix<double> dDeltaAlphadEpsilon = mStructure->CalculateDDeltaCrackAngleElastic();
		NuTo::FullMatrix<double> d2DeltaAlphadEpsilon = mStructure->CalculateD2DeltaCrackAngleElastic();

		//std::cout << "dDeltaAlphadEpsilon " << std::endl;
		//dDeltaAlphadEpsilon.Info(12,3);

		//std::cout << "d2DeltaAlphadEpsilon " << std::endl;
		//d2DeltaAlphadEpsilon.Info(12,3);

		//std::cout << "scalingFactorAlpha "<< mStructure->GetScalingFactorAlpha() << std::endl;
		//std::cout << "scalingFactorEpsilon "<< mStructure->GetScalingFactorEpsilon() << std::endl;
		rResult.AddEntry(0,0, mPenaltyStiffness*mStructure->GetScalingFactorCrackAngle()*mStructure->GetScalingFactorCrackAngle());
		rResult.AddEntry(0,1, -mPenaltyStiffness*mStructure->GetScalingFactorEpsilon()*mStructure->GetScalingFactorCrackAngle()*dDeltaAlphadEpsilon(0,0));
		rResult.AddEntry(0,2, -mPenaltyStiffness*mStructure->GetScalingFactorEpsilon()*mStructure->GetScalingFactorCrackAngle()*dDeltaAlphadEpsilon(1,0));
		rResult.AddEntry(0,3, -mPenaltyStiffness*mStructure->GetScalingFactorEpsilon()*mStructure->GetScalingFactorCrackAngle()*dDeltaAlphadEpsilon(2,0));
		rResult.AddEntry(1,1, mPenaltyStiffness*mStructure->GetScalingFactorEpsilon()*mStructure->GetScalingFactorEpsilon()*(dDeltaAlphadEpsilon(0,0)*dDeltaAlphadEpsilon(0,0)-delta_alpha*d2DeltaAlphadEpsilon(0,0)));
		rResult.AddEntry(1,2, mPenaltyStiffness*mStructure->GetScalingFactorEpsilon()*mStructure->GetScalingFactorEpsilon()*(dDeltaAlphadEpsilon(0,0)*dDeltaAlphadEpsilon(1,0)-delta_alpha*d2DeltaAlphadEpsilon(0,1)));
		rResult.AddEntry(1,3, mPenaltyStiffness*mStructure->GetScalingFactorEpsilon()*mStructure->GetScalingFactorEpsilon()*(dDeltaAlphadEpsilon(0,0)*dDeltaAlphadEpsilon(2,0)-delta_alpha*d2DeltaAlphadEpsilon(0,2)));
		rResult.AddEntry(2,2, mPenaltyStiffness*mStructure->GetScalingFactorEpsilon()*mStructure->GetScalingFactorEpsilon()*(dDeltaAlphadEpsilon(1,0)*dDeltaAlphadEpsilon(1,0)-delta_alpha*d2DeltaAlphadEpsilon(1,1)));
		rResult.AddEntry(2,3, mPenaltyStiffness*mStructure->GetScalingFactorEpsilon()*mStructure->GetScalingFactorEpsilon()*(dDeltaAlphadEpsilon(1,0)*dDeltaAlphadEpsilon(2,0)-delta_alpha*d2DeltaAlphadEpsilon(1,2)));
		rResult.AddEntry(3,3, mPenaltyStiffness*mStructure->GetScalingFactorEpsilon()*mStructure->GetScalingFactorEpsilon()*(dDeltaAlphadEpsilon(2,0)*dDeltaAlphadEpsilon(2,0)-delta_alpha*d2DeltaAlphadEpsilon(2,2)));

	/*    std::cout << "constraint stiffness alpha analytical " << std::endl;
		NuTo::FullMatrix<double> stiffness_analytic(rResult);
		stiffness_analytic.Info(12,7);

		//check stiffness
		{
		double delta(1e-10);
		NuTo::FullMatrix<double> stiffness_analytic(rResult);
		NuTo::FullMatrix<double> stiffness_cdf(rResult);
		NuTo::FullMatrix<double> gradient1, gradient2;
		CalculateGradientInternalPotential(gradient1, rGlobalDofs);
		double alpha = mStructure->GetCrackAngle();
		alpha+=delta;
		mStructure->SetCrackAngle(alpha);
		CalculateGradientInternalPotential(gradient2, rGlobalDofs);
		alpha-=delta;
		mStructure->SetCrackAngle(alpha);
		stiffness_cdf.SetColumn(0,(gradient2-gradient1)*(mStructure->GetScalingFactorCrackAngle()/delta));

		EngineeringStrain2D strainTensor(mStructure->GetTotalEngineeringStrain());
		const double* dataPtr = strainTensor.GetData();

		for (int count=0; count<3; count++)
		{
			const_cast<double*>(dataPtr)[count]+=delta;
			const_cast<StructureMultiscale*>(mStructure)->SetTotalEngineeringStrain(strainTensor);
			CalculateGradientInternalPotential(gradient2, rGlobalDofs);
			const_cast<double*>(dataPtr)[count]-=delta;
			const_cast<StructureMultiscale*>(mStructure)->SetTotalEngineeringStrain(strainTensor);
			stiffness_cdf.SetColumn(count+1,(gradient2-gradient1)*(mStructure->GetScalingFactorEpsilon()/delta));
		}
		std::cout << "constraint stiffness alpha cdf        " << std::endl;
		stiffness_cdf.Info(12,5);
		}
	*/
	}
	else
	{
		//crack angle
		int dof(1);
		rResult.Resize(dof,dof);
		rGlobalDofs.resize(dof,1);

		rGlobalDofs[0] = mStructure->GetDofCrackAngle();

		rResult.AddEntry(0,0, mPenaltyStiffness*mStructure->GetScalingFactorCrackAngle()*mStructure->GetScalingFactorCrackAngle());
	}
}

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::ConstraintNonlinearGlobalCrackAngle2D::CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
	if (mCoupleToTotalStrain)
	{
		//crack angle
		int dof(4);
		rResult.Resize(dof,1);
		rGlobalDofs.resize(dof,1);

		rGlobalDofs[0] = mStructure->GetDofCrackAngle();
		rGlobalDofs[1] = mStructure->GetDofGlobalTotalStrain2D()[0];
		rGlobalDofs[2] = mStructure->GetDofGlobalTotalStrain2D()[1];
		rGlobalDofs[3] = mStructure->GetDofGlobalTotalStrain2D()[2];

		//calculate angle orthogonal to second principal stress
		double delta_alpha = mStructure->CalculateDeltaCrackAngleElastic();

		NuTo::FullMatrix<double> dDeltaAlphadEpsilon = mStructure->CalculateDDeltaCrackAngleElastic();

		//std::cout << "dDeltaAlphadEpsilon" << std::endl;
		//dDeltaAlphadEpsilon.Trans().Info(12,5);
		//std::cout << "delta alpha " << delta_alpha << std::endl;

		rResult(0,0) =  mPenaltyStiffness * delta_alpha * mStructure->GetScalingFactorCrackAngle();
		rResult(1,0) = -mPenaltyStiffness * delta_alpha * dDeltaAlphadEpsilon(0,0)*mStructure->GetScalingFactorEpsilon();
		rResult(2,0) = -mPenaltyStiffness * delta_alpha * dDeltaAlphadEpsilon(1,0)*mStructure->GetScalingFactorEpsilon();
		rResult(3,0) = -mPenaltyStiffness * delta_alpha * dDeltaAlphadEpsilon(2,0)*mStructure->GetScalingFactorEpsilon();

		//std::cout << "constraint gradient alpha analytical " << std::endl;
		//rResult.Info(12,8);

		/*
		//check the result
		double delta(1e-9);
		NuTo::FullMatrix<double> gradient_cdf(rResult);
		double pot1, pot2;
		pot1 = CalculateTotalPotential();
		double alpha = mStructure->GetCrackAngle();
		alpha+=delta;
		mStructure->SetCrackAngle(alpha);
		pot2 = CalculateTotalPotential();
		alpha-=delta;
		mStructure->SetCrackAngle(alpha);
		gradient_cdf(0,0) = (pot2-pot1)*(mStructure->GetScalingFactorAlpha()/delta);

		EngineeringStrain2D strainTensor(mStructure->GetTotalEngineeringStrain());
		const double* dataPtr = strainTensor.GetData();
		delta=1e-9*sqrt(dataPtr[0]*dataPtr[0]+dataPtr[1]*dataPtr[1]+0.5*dataPtr[2]*dataPtr[2]);

		for (int count=0; count<3; count++)
		{
			const_cast<double*>(dataPtr)[count]+=delta;
			const_cast<StructureMultiscale*>(mStructure)->SetTotalEngineeringStrain(strainTensor);
			pot2 = CalculateTotalPotential();
			const_cast<double*>(dataPtr)[count]-=delta;
			const_cast<StructureMultiscale*>(mStructure)->SetTotalEngineeringStrain(strainTensor);
			gradient_cdf(count+1,0) = (pot2-pot1)*(mStructure->GetScalingFactorEpsilon()/delta);
		}
		std::cout << "constraint gradient alpha cdf        " << std::endl;
		gradient_cdf.Info(12,8);
	*/
	}
	else
	{
		//crack angle
		int dof(1);
		rResult.Resize(dof,1);
		rGlobalDofs.resize(dof,1);

		rGlobalDofs[0] = mStructure->GetDofCrackAngle();

		//calculate delta of crack angle with respect to previous crack angle
		double delta_alpha = mStructure->CalculateDeltaCrackAnglePrev();

		rResult(0,0) =  mPenaltyStiffness * delta_alpha * mStructure->GetScalingFactorCrackAngle();
	}
}

//! @brief calculates the internal potential
double NuTo::ConstraintNonlinearGlobalCrackAngle2D::CalculateTotalPotential()const
{
	if (mCoupleToTotalStrain)
	{
		//calculate angle orthogonal to second principal stress
		double delta_alpha = mStructure->CalculateDeltaCrackAngleElastic();
		double delta_alpha_square=delta_alpha*delta_alpha;
		return 0.5*mPenaltyStiffness*delta_alpha_square;
	}
	else
	{
		double delta_alpha = mStructure->CalculateDeltaCrackAnglePrev();
		double delta_alpha_square=delta_alpha*delta_alpha;
		return 0.5*mPenaltyStiffness*delta_alpha_square;
	}
}


#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintNonLinearGlobalCrackAngle2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNonlinear)
       & BOOST_SERIALIZATION_NVP(mPenaltyStiffness)
       & BOOST_SERIALIZATION_NVP(mStructure);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintNonLinearGlobalCrackAngle2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintNonlinearGlobalCrackAngle2D)
#endif // ENABLE_SERIALIZATION
