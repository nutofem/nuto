// $Id: JumpDirect.cpp 2014-11-06 12:55:55Z vkindrac $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#endif // ENABLE_SERIALIZATION

# ifdef _OPENMP
#include <omp.h>
# endif

#include <iostream>
#include <fstream>

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverMKLPardiso.h"

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/JumpDirect.h"
#include "nuto/mechanics/timeIntegration/ResultGroupNodeForce.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataDamageViscoPlasticity3D.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::JumpDirect::JumpDirect (StructureBase* rStructure)  : NewmarkDirect (rStructure)
{
    mMinLineSearchStep = 0.01;
    mHarmonicIncrementation = 16;
    mHarmonicExtrapolation = false;
    mHarmonicExcitation = false;
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::JumpDirect::Info()const
{
    NewmarkBase::Info();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::JumpDirect::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::JumpDirect::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::JumpDirect::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::JumpDirect::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::JumpDirect::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::JumpDirect::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::JumpDirect::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of JumpDirect" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NewmarkBase)
           & BOOST_SERIALIZATION_NVP(mMinLineSearchStep)
           & BOOST_SERIALIZATION_NVP(mHarmonicIncrementation)
           & BOOST_SERIALIZATION_NVP(mHarmonicExtrapolation)
           & BOOST_SERIALIZATION_NVP(mHarmonicExcitation);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of JumpDirect" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


NuTo::Error::eError NuTo::JumpDirect::Solve(double rTimeDelta)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    //allocate solver
    NuTo::SparseDirectSolverMUMPS mySolver;
    //NuTo::SparseDirectSolverMKLPardiso mySolver;
#ifdef SHOW_TIME
        mySolver.SetShowTime(mStructure->GetShowTime());
#endif
    try
    {
    	this->SetNewmarkBeta(15.);
    	//calculate the end time of monotonic loading which is exactly the beginning of cyclic loading
    	NuTo::Error::eError Error;
    	double BeginHarmonicExcitation(mTimeDependentConstraintFactor(mTimeDependentConstraintFactor.GetNumRows()-1,0));

    	//calculate till the end of the first three cycle
    	if (mHarmonicExtrapolation) {
    		//calculate the first three cycles
    		Error = Solve(BeginHarmonicExcitation + 3./mHarmonicConstraintFactor(0,1));
    		if (Error != NuTo::Error::SUCCESSFUL) {
    			return Error;
    		}
		} else {
			//calculate the whole loading history
			return NuTo::NewmarkDirect::Solve(rTimeDelta);
		}

    	//calculate matrices
        //calculate constraint matrix
        NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
        mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
        NuTo::SparseMatrixCSRVector2General<double> Cmat(CmatTmp);
        SparseMatrixCSRVector2General<double> CmatT (Cmat.Transpose());
        FullVector<double,Eigen::Dynamic> bRHSprev, bRHSend;

        // allocate space for stiffness matrix
        SparseMatrixCSRVector2General<double> stiffMatrix_jj(mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
        SparseMatrixCSRVector2General<double> stiffMatrix_jk(mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
        SparseMatrixCSRVector2General<double> stiffMatrix_kj(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
        SparseMatrixCSRVector2General<double> stiffMatrix_kk(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());

        //calculate individual mass matrix
        NuTo::SparseMatrixCSRVector2General<double> massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk;

		massMatrix_jj.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
		massMatrix_jk.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
		massMatrix_kj.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
		massMatrix_kk.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
		mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);

        NuTo::FullVector<double,Eigen::Dynamic> extForce_j, extForce_k, prevExtForce_j, prevExtForce_k;
        NuTo::FullVector<double,Eigen::Dynamic> intForce_j(mStructure->GetNumActiveDofs()), intForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs()),
                                 prevIntForce_j(mStructure->GetNumActiveDofs()), prevIntForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
        NuTo::FullVector<double,Eigen::Dynamic> residual_j, residual_k, residual_mod, prevResidual_j, prevResidual_k;

        NuTo::FullVector<double,Eigen::Dynamic> lastConverged_disp_j,lastConverged_disp_k, lastConverged_vel_j, lastConverged_vel_k, lastConverged_acc_j, lastConverged_acc_k;

//    	//Example get number of elements
//    	std::cout << "Number of Elements = " << mStructure->GetNumElements() << std::endl;

//    	//Example get stress
//    	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rEngineeringStress;
//    	mStructure->ElementGetEngineeringStress(0, rEngineeringStress);
//    	std::cout << "Return stress " << std::endl;
//    	std::cout << rEngineeringStress << std::endl;

//    	//Example get damage
//    	if (mStructure->ElementGetElementPtr(0)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS){
//    		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(0)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    		std::cout << "Damage = " << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl;
//    	}

        // save static data and displacements after the 3th cycle
    	NuTo::FullVector<double,Eigen::Dynamic> save_disp_j, save_disp_k, delta_disp_j, delta_disp_k;;
    	double SaveTime(mStructure->GetTime());
    	mStructure->NodeExtractDofValues(0,save_disp_j, save_disp_k);
    	mStructure->ElementFatigueSaveStaticData();

//    	// calculate 10 cycles
//    	std::cout << "Calculate 10 Cycles from here, mTime = " << mTime << std::endl;
//    	Error = NuTo::NewmarkDirect::Solve(SaveTime + 10./mHarmonicConstraintFactor(0,1));

//    	//Example get damage
//    	if (mStructure->ElementGetElementPtr(0)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS){
//    		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(0)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    		std::cout << "Damage after 10 cycles = " << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl;
//    	}

//    	// repeat with restoring
//    	mStructure->SetPrevTime(SaveTime);
//    	mStructure->SetTime(SaveTime);
//    	// postprocessing time
//    	mTime = SaveTime;
//    	mStructure->NodeMergeActiveDofValues(0,save_disp_j);
//    	mStructure->ElementFatigueRestoreStaticData();
//    	std::cout << " SaveTime " << SaveTime << ", mStructure->GetTime() = " << mStructure->GetTime() << std::endl;
//
//    	//Example get damage
//    	if (mStructure->ElementGetElementPtr(0)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS){
//    		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(0)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    		std::cout << "Damage after restoring " << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl;
//    	}
//
//    	std::cout << "Repeat 10 Cycles from here, mTime = " << mTime << std::endl;
//    	Error = NuTo::NewmarkDirect::Solve(SaveTime + 10./mHarmonicConstraintFactor(0,1));

    	// perform the 4th cycle
    	// initialize mean displacement and displacement amplitude
    	NuTo::FullVector<double,Eigen::Dynamic> disp_Max_j, disp_Max_k, disp_Min_j, disp_Min_k;

    	Error = NuTo::NewmarkDirect::Solve(SaveTime + 0.25/mHarmonicConstraintFactor(0,1));
    	mStructure->NodeExtractDofValues(0,disp_Max_j, disp_Max_k);

    	Error = NuTo::NewmarkDirect::Solve(SaveTime + 0.75/mHarmonicConstraintFactor(0,1));
    	mStructure->NodeExtractDofValues(0,disp_Min_j, disp_Min_k);

    	NuTo::FullVector<double,Eigen::Dynamic> disp_Mean_j, disp_Mean_k, disp_Ampl_j, disp_Ampl_k;

    	Error = NuTo::NewmarkDirect::Solve(SaveTime + 1./mHarmonicConstraintFactor(0,1));
    	mStructure->NodeExtractDofValues(0,disp_Mean_j, disp_Mean_k);

    	disp_Ampl_j = 0.5*(disp_Max_j - disp_Min_j);
    	disp_Ampl_k = 0.5*(disp_Max_k - disp_Min_k);

        if (mStructure->GetNumTimeDerivatives()>1)
	    {
	        mStructure->NodeExtractDofValues(1,lastConverged_vel_j,lastConverged_vel_k);
	        mStructure->NodeExtractDofValues(2,lastConverged_acc_j,lastConverged_acc_k);
			std::cout << "=== acc_j = " << lastConverged_acc_j[0] << std::endl;

	    }

    	// extrapolate
    		// calculate the number of cycles to be extrapolated
	    	// initiate number of cycles
    	int Njump(0);
    	Njump = 10;
    	do
    	{
			mStructure->SetNumExtrapolatedCycles(Njump);
				// extrapolate time starting from the end of the 4th cycle (SaveTime is still the end of the 3d cycle)
			mStructure->SetTime(SaveTime + (Njump + 1)/mHarmonicConstraintFactor(0,1));
			mStructure->SetPrevTime(SaveTime + (Njump + 1)/mHarmonicConstraintFactor(0,1));
	//    	mTime += Njump/mHarmonicConstraintFactor(0,1);

				// extrapolate state variables
			mStructure->ElementFatigueExtrapolateStaticData();

            if (mStructure->GetNumTimeDerivatives()>1)
			{
				mStructure->NodeExtractDofValues(1,lastConverged_vel_j,lastConverged_vel_k);
				mStructure->NodeExtractDofValues(2,lastConverged_acc_j,lastConverged_acc_k);
				std::cout << "=== acc_j = " << lastConverged_acc_j[0] << std::endl;
			}

			// find equilibrium
	//    		// find for the mean value (rFourier = 0)
	//    		// set BC for the mean displacement, which is the displacement at the beginning of the harmonic excitation
	//		double timeDependentConstraintFactor(this->CalculateTimeDependentConstraintFactor(BeginHarmonicExcitation));
	//		mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
	//		mStructure->ConstraintGetRHSAfterGaussElimination(bRHSend);
	//			// set mean displacement field
	//		mStructure->NodeMergeActiveDofValues(0,disp_Mean_j);
	//		mStructure->ElementTotalUpdateTmpStaticData();
	//			// calculate internal forces and external forces
	//		mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(intForce_j,intForce_k);
	//		CalculateExternalLoad(*mStructure, BeginHarmonicExcitation, extForce_j, extForce_k);
	//
	//		residual_j = intForce_j - extForce_j;
	//		residual_k = intForce_k - extForce_k;
	//
	//        if (Cmat.GetNumEntries()>0)
	//        {
	//            residual_mod = residual_j - CmatT*residual_k;
	//        } else {
	//            residual_mod = residual_j;
	//        }
	//
	//        this->CalculateGlobalModifiedStiffness(&stiffMatrix_jj,0);
	//
	//		// solve for trial state
	//		NuTo::SparseMatrixCSRGeneral<double> hessianModSolver(stiffMatrix_jj);
	//
	//		hessianModSolver.SetOneBasedIndexing();
	//		mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
	//		delta_disp_j*=-1;
	//
	//        // calculate mean displacement
	//		disp_Mean_j += delta_disp_j;
	//		disp_Mean_k = bRHSend - Cmat*disp_Mean_j;

			this->CalculateFourierCofficients(&disp_Mean_j,&disp_Mean_k,&disp_Ampl_j,&disp_Ampl_k);

            if (mStructure->GetNumTimeDerivatives()>1)
			{
				mStructure->NodeExtractDofValues(1,lastConverged_vel_j,lastConverged_vel_k);
				mStructure->NodeExtractDofValues(2,lastConverged_acc_j,lastConverged_acc_k);
				std::cout << "=== acc_j = " << lastConverged_acc_j[0] << std::endl;

			}

			// store static data for the next jump
			mStructure->ElementFatigueSaveStaticData();

			// straight forward integration after extrapolation
			mStructure->NodeMergeActiveDofValues(0,disp_Mean_j);
			mStructure->ElementTotalUpdateTmpStaticData();
	//		std::cout << "mStress before USD" << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetPrevStress() << std::endl;
	//		std::cout << "mOmegaCompr before USD" << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl;
				// the state variables are already extrapolated, except of mPrevStrain and mPrevSigma
				// call UpdateStaticDatat in order to update mPrevStrain and mPrevSigma too (this has to be done after the equilibrium has been found)
			mStructure->ElementTotalUpdateStaticData();
	//		std::cout << "mStress after USD " << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetPrevStress() << std::endl;
	//		std::cout << "mOmegaCompr after USD" << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl;
			SaveTime = mStructure->GetTime();
			this->IntegrateSingleCycle(&disp_Mean_j,&disp_Mean_k,&disp_Ampl_j,&disp_Ampl_k,true);
	//	    if (this->IsDynamic())
	//	    {
	//	        mStructure->NodeExtractDofValues(1,lastConverged_vel_j,lastConverged_vel_k);
	//	        mStructure->NodeExtractDofValues(2,lastConverged_acc_j,lastConverged_acc_k);
	//			std::cout << "=== acc_j = " << lastConverged_acc_j[0] << std::endl;
	//
	//	    }
	//		Error = NuTo::NewmarkDirect::Solve(mStructure->GetTime() + 1./mHarmonicConstraintFactor(0,1));
    	}
    	while (SaveTime + (Njump + 1)/mHarmonicConstraintFactor(0,1) <= rTimeDelta);
    }
//    try
//    {
//    	double timeDependentConstraintTime(mTimeDependentConstraintFactor(mTimeDependentConstraintFactor.GetNumRows()-1,0));
//    	std::cout << "... applied harmonic = " << mHarmonicExcitation << std::endl;
//    	std::cout << "... displecement amplitude = " << mHarmonicConstraintFactor(0,0) << std::endl;
//    	std::cout << "... displecement frequency = " << mHarmonicConstraintFactor(0,1) << std::endl;
//    	std::cout << "... number of cycles = " << int (mHarmonicConstraintFactor(0,2)) << std::endl;
//    	std::cout << "... timeDependentConstraintTime = " << timeDependentConstraintTime << std::endl;
//
//        if (mMaxTimeStep==0)
//            throw MechanicsException("[NuTo::JumpDirect::Solve] max time step is set to zero.");
//
//        //renumber dofs and build constraint matrix
//        mStructure->NodeBuildGlobalDofs();
//
//        //calculate constraint matrix
//        NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
//        mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
//        NuTo::SparseMatrixCSRVector2General<double> Cmat(CmatTmp);
//        SparseMatrixCSRVector2General<double> CmatT (Cmat.Transpose());
//        FullVector<double,Eigen::Dynamic> bRHSprev, bRHShalf, bRHSend, bRHSdot, bRHSddot;
//
//        // allocate space for stiffness matrix
//        SparseMatrixCSRVector2General<double> stiffMatrix_jj(mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
//        SparseMatrixCSRVector2General<double> stiffMatrix_jk(mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
//        SparseMatrixCSRVector2General<double> stiffMatrix_kj(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
//        SparseMatrixCSRVector2General<double> stiffMatrix_kk(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
//
//        //calculate individual mass matrix (I was lazy here and have just used general matrices, but symmetric is certainly better here)
//        NuTo::SparseMatrixCSRVector2General<double> massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk;
//        NuTo::FullVector<double, Eigen::Dynamic> lumped_massMatrix_j(mStructure->GetNumDofs()),
//        		                                 lumped_massMatrix_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
//        if (this->IsDynamic())
//        {
//        	if (mUseLumpedMass)
//        	{
//        		mStructure->BuildGlobalLumpedHession2(lumped_massMatrix_j,lumped_massMatrix_k);
//        		std::cout << "use lumped mass matrix " << std::endl;
//        	}
//        	else
//        	{
//				massMatrix_jj.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
//				massMatrix_jk.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
//				massMatrix_kj.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
//				massMatrix_kk.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
//				mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);
//        		std::cout << "use full mass matrix " << std::endl;
//        	}
//            //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> m11Full(massMatrix_jj);
//            //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenValues;
//            //m11Full.EigenValuesSymmetric(eigenValues);
//            //std::cout << "eigenvalues mass" << "\n";
//            //std::cout << eigenValues.Trans() << "\n";
//        }
//
//        //extract displacements, velocities and accelerations
//        NuTo::FullVector<double,Eigen::Dynamic> disp_j,disp_k, vel_j, vel_k, acc_j, acc_k, delta_disp_j, delta_disp_k;
//        NuTo::FullVector<double,Eigen::Dynamic> trial_disp_j,trial_disp_k, trial_vel_j, trial_vel_k, trial_acc_j, trial_acc_k;
//        NuTo::FullVector<double,Eigen::Dynamic> extForce_j, extForce_k, prevExtForce_j, prevExtForce_k;
//        NuTo::FullVector<double,Eigen::Dynamic> intForce_j(mStructure->GetNumActiveDofs()), intForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs()),
//                                 prevIntForce_j(mStructure->GetNumActiveDofs()), prevIntForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
//        NuTo::FullVector<double,Eigen::Dynamic> residual_j, residual_k, residual_mod, prevResidual_j, prevResidual_k;
//
//        //store last converged displacements, velocities and accelerations
//        NuTo::FullVector<double,Eigen::Dynamic> lastConverged_disp_j,lastConverged_disp_k, lastConverged_vel_j, lastConverged_vel_k, lastConverged_acc_j, lastConverged_acc_k;
//        mStructure->NodeExtractDofValues(0,lastConverged_disp_j, lastConverged_disp_k);
//
//        if (this->IsDynamic())
//        {
//            mStructure->NodeExtractDofValues(1,lastConverged_vel_j,lastConverged_vel_k);
//            mStructure->NodeExtractDofValues(2,lastConverged_acc_j,lastConverged_acc_k);
//        }
//
//        double curTime  = 0;
//
//        // initialize the structure times
//        mStructure->SetPrevTime(curTime);
//        mStructure->SetTime(curTime);
//
//        //apply constraints for last converged time step
//        double timeDependentConstraintFactor(0);
//        if (mTimeDependentConstraint!=-1)
//        {
//            timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
//            mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
//        }
//        mStructure->ConstraintGetRHSAfterGaussElimination(bRHSprev);
//        mStructure->NodeMergeActiveDofValues(0,lastConverged_disp_j); //disp_k is internally calculated from the previous constraint matrix
//        mStructure->ElementTotalUpdateTmpStaticData();
//
//        //calculate internal force, update of history variables=false
//        mStructure->BuildGlobalGradientInternalPotentialSubVectors(prevIntForce_j,prevIntForce_k,false);
//        prevResidual_j = prevIntForce_j;
//        prevResidual_k = prevIntForce_k;
//
//        if (this->IsDynamic())
//        {
//            if (mMuDampingMass>0)
//            {
//                //add damping terme
//            	if (mUseLumpedMass)
//            	{
//            		prevResidual_j += lumped_massMatrix_j.asDiagonal()*lastConverged_vel_j * mMuDampingMass;
//            		prevResidual_k += lumped_massMatrix_k.asDiagonal()*lastConverged_vel_k * mMuDampingMass;
//            	}
//            	else
//            	{
//					prevResidual_j += (massMatrix_jj*lastConverged_vel_j+massMatrix_jk*lastConverged_vel_k) *mMuDampingMass;
//					prevResidual_k += (massMatrix_kj*lastConverged_vel_j+massMatrix_kk*lastConverged_vel_k) *mMuDampingMass;
//            	}
//            }
//        	if (mUseLumpedMass)
//        	{
//        		prevResidual_j += lumped_massMatrix_j.asDiagonal()*lastConverged_acc_j;
//        		prevResidual_k += lumped_massMatrix_k.asDiagonal()*lastConverged_acc_k;
//        	}
//        	else
//        	{
//				//add mass terme
//				prevResidual_j += (massMatrix_jj*lastConverged_acc_j+massMatrix_jk*lastConverged_acc_k);
//				prevResidual_k += (massMatrix_kj*lastConverged_acc_j+massMatrix_kk*lastConverged_acc_k);
//        	}
//        }
//        //add external force
//        CalculateExternalLoad(*mStructure, curTime, prevExtForce_j, prevExtForce_k);
//        prevResidual_j -= prevExtForce_j;
//        prevResidual_k -= prevExtForce_k;
//
//        if (Cmat.GetNumEntries()>0)
//        {
//            residual_mod=prevResidual_j - CmatT*prevResidual_k;
//        }
//        else
//        {
//            residual_mod=prevResidual_j;
//        }
//
//        if (residual_mod.Norm()>mToleranceForce)
//            throw MechanicsException("[NuTo::JumpDirect::Solve] Initial configuration is not in (dynamic) equilibrium.");
//        std::cout << "residual in initial configuration " << residual_mod.Norm() << std::endl;
//
///*        FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> plotVector0(plotHistory);
//        for (int countGroup=0; countGroup<mVecGroupNodesReactionForces.GetNumRows(); countGroup++)
//        {
//            FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> reactionForce(1,mStructure->GetDimension());
//
//            boost::ptr_map<int,GroupBase>::iterator itGroup = mStructure->mGroupMap.find(mVecGroupNodesReactionForces(countGroup,0));
//            if (itGroup==mStructure->mGroupMap.end())
//                throw MechanicsException("[NuTo::NewmarkDirect::Solve] node group with the given identifier for the reaction forces does not exist.");
//            if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
//                throw MechanicsException("[NuTo::NewmarkDirect::Solve] Group is not a node group (reaction forces).");
//            Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();
//            assert(nodeGroup!=0);
//
//            //all nodes have to have the same dimension
//            if(nodeGroup->GetNumMembers()<1)
//                throw MechanicsException("[NuTo::StructureBase::NodeGroupGetCoordinates] Group has no members.");
//
//            //std::cout << "residual k " << residual_k.transpose() << std::endl;
//            for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++)
//            {
//                for (int countDimension=0; countDimension<mStructure->GetDimension(); countDimension++)
//                {
//                    int theDof = itNode->second->GetDofDisplacement(countDimension);
//                    if (theDof<mStructure->GetNumActiveDofs())
//                    {
//                        reactionForce(0,countDimension)+=prevResidual_j(theDof);
//                    }
//                    else
//                    {
//                        reactionForce(0,countDimension)+=prevResidual_k(theDof-mStructure->GetNumActiveDofs());
//                    }
//                }
//            }
//            //std::cout << "reaction force for node group " <<  countGroup << " : " << reactionForce.Trans() << std::endl;
//            plotVector0.AppendColumns(reactionForce);
//        }
//*/
//        PostProcess(prevResidual_j, prevResidual_k);
//
//        double timeStep = mTimeStep;
//        while (curTime < rTimeDelta)
//        {
//            //apply constraints for last converged time step
//            if (mTimeDependentConstraint!=-1)
//            {
//				timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
//				mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
//            }
//            mStructure->ConstraintGetRHSAfterGaussElimination(bRHSprev);
//            mStructure->NodeMergeActiveDofValues(0,lastConverged_disp_j); //disp_k is internally calculated from the previous constraint matrix
//            mStructure->ElementTotalUpdateTmpStaticData();
//
//            //add external force
//            CalculateExternalLoad(*mStructure, curTime, prevExtForce_j, prevExtForce_k);
//
//            //increase time step
//            curTime += timeStep;
//
//            if (timeStep<mMinTimeStep)
//                throw MechanicsException("[NuTo::JumpDirect::Solve] time step is smaller than minimum - no convergence is obtained.");
//
//            // check wether a harmonic time step and check wether curTime is too close to the time data
//            this->SetTimeAndTimeStep(curTime, timeStep, rTimeDelta);
//
//            // set new structure time at the end of the time increment
//            mStructure->SetTime(curTime);
//
//            if (this->IsDynamic())
//            {
//                if (1.-mGamma*mMuDampingMass*timeStep<=0)
//                    throw MechanicsException("[NuTo::JumpDirect::Solve] Factor for the mass matrix is negative - reduce your time step.");
//            }
//
//            //calculate the initial out-of-balance force
//            disp_j = lastConverged_disp_j;
//            disp_k = lastConverged_disp_k;
//
//            //calculate internal force, update of history variables=false
//            mStructure->BuildGlobalGradientInternalPotentialSubVectors(prevIntForce_j,prevIntForce_k,false);
//            residual_j = prevIntForce_j;
//            residual_k = prevIntForce_k;
//
//            //increase time step by half and calculate the constraint matrix (used to approximate the velocity and acceleration of the rhs)
//            if (this->IsDynamic())
//            {
//                if (mTimeDependentConstraint!=-1)
//                {
//					timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
//					mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
//                }
//                mStructure->ConstraintGetRHSAfterGaussElimination(bRHShalf);
//            }
//
//            //apply constraints for the new time step (modified bRHS)
//            if (mTimeDependentConstraint!=-1)
//            {
//				timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
//				mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
//            }
//            mStructure->ConstraintGetRHSAfterGaussElimination(bRHSend);
//            FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deltaBRHS(bRHSend-bRHSprev);
//
//            //std::cout << "bRHSprev " << bRHSprev.Trans() << "\n";
//            //std::cout << "bRHShalf " << bRHShalf.Trans() << "\n";
//            //std::cout << "bRHSend " << bRHSend.Trans() << "\n";
//            //std::cout << "bRHSdot " << bRHSdot.Trans() << "\n";
//            //std::cout << "bRHSddot " << bRHSddot.Trans() << "\n";
//
//            if (this->IsDynamic())
//            {
//                //calculate approximations to the time derivates of the rhs of the constraint matrix
//                bRHSdot = (bRHSprev*5.-bRHShalf*8.+bRHSend*3.)*(-1./(timeStep));
//                bRHSddot = (bRHSprev-bRHShalf*2+bRHSend)*(4./(timeStep*timeStep));
//
//                //calculate new accelerations and velocities of independent dofs
//                acc_j = (lastConverged_vel_j+lastConverged_acc_j*(timeStep*(0.5-mBeta)))*(-1./(timeStep*mBeta));
//                vel_j = lastConverged_vel_j+lastConverged_acc_j*((1.-mGamma)*timeStep)+acc_j*(mGamma*timeStep);
//
//                //calculate new accelerations and velocities of dependent dofs
//                acc_k = bRHSddot - (Cmat*acc_j);
//                vel_k = bRHSdot - (Cmat*vel_j);
//
//				if (mMuDampingMass>0)
//				{
//					//add damping terme
//					if (mUseLumpedMass)
//					{
//						residual_j += lumped_massMatrix_j.asDiagonal()*vel_j * mMuDampingMass;
//						residual_k += lumped_massMatrix_k.asDiagonal()*vel_k * mMuDampingMass;
//					}
//					else
//					{
//	                    //add damping terme
//	                    residual_j += (massMatrix_jj*vel_j+massMatrix_jk*vel_k) *mMuDampingMass;
//	                    residual_k += (massMatrix_kj*vel_j+massMatrix_kk*vel_k) *mMuDampingMass;
//					}
//				}
//				if (mUseLumpedMass)
//				{
//					residual_j += lumped_massMatrix_j.asDiagonal()*acc_j;
//					residual_k += lumped_massMatrix_k.asDiagonal()*acc_k;
//				}
//				else
//				{
//	                //add mass terme
//	                residual_j += (massMatrix_jj*acc_j+massMatrix_jk*acc_k);
//	                residual_k += (massMatrix_kj*acc_j+massMatrix_kk*acc_k);
//				}
//            }
//            //add external force
//            CalculateExternalLoad(*mStructure, curTime, extForce_j, extForce_k);
//            residual_j -= extForce_j;
//            residual_k -= extForce_k;
//
//            //calculate stiffness and disp force vector (the displacements of the dependent dofs do not take into account the modified rhs)
//            if (Cmat.GetNumEntries()>0)
//            {
//                stiffMatrix_jj.SetZeroEntries();stiffMatrix_jk.SetZeroEntries();stiffMatrix_kj.SetZeroEntries();stiffMatrix_kk.SetZeroEntries();
//                mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::STIFFNESS, stiffMatrix_jj, stiffMatrix_jk, stiffMatrix_kj, stiffMatrix_kk);
//                if (this->IsDynamic())
//                {
//                    double factor((1.+mGamma*mMuDampingMass*timeStep)/(timeStep*timeStep*mBeta));
//    				if (mUseLumpedMass)
//    				{
//                        stiffMatrix_jj.AddScalDiag(lumped_massMatrix_j,factor);
//                        stiffMatrix_kk.AddScalDiag(lumped_massMatrix_k,factor);
//    				}
//    				else
//    				{
//                        stiffMatrix_jj.AddScal(massMatrix_jj,factor);
//                        stiffMatrix_jk.AddScal(massMatrix_jk,factor);
//                        stiffMatrix_kj.AddScal(massMatrix_kj,factor);
//                        stiffMatrix_kk.AddScal(massMatrix_kk,factor);
//    				}
//                }
//                //calculate the residual for the zero time step
//                residual_j += stiffMatrix_jk * deltaBRHS;
//                residual_k += stiffMatrix_kk * deltaBRHS;
//                residual_mod=residual_j - CmatT*residual_k;
//
//                //add damping and mass to full hessian
//                stiffMatrix_jj -= CmatT * stiffMatrix_kj + stiffMatrix_jk * Cmat;
//                stiffMatrix_jj += CmatT * stiffMatrix_kk * Cmat;
//            }
//            else
//            {
//                stiffMatrix_jj.SetZeroEntries();stiffMatrix_jk.SetZeroEntries();
//                mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::STIFFNESS, stiffMatrix_jj, stiffMatrix_jk);
//
//                //calculate the residual for the zero time step
//                residual_j += stiffMatrix_jk * deltaBRHS;
//                residual_mod=residual_j;
//
//                //add damping and mass to full hessian
//                if (this->IsDynamic())
//                {
//                    double factor((1.+mGamma*mMuDampingMass*timeStep)/(timeStep*timeStep*mBeta));
//    				if (mUseLumpedMass)
//    				{
//    					stiffMatrix_jj.AddScalDiag(lumped_massMatrix_j,factor);
//    				}
//    				else
//    				{
//                        stiffMatrix_jj.AddScal(massMatrix_jj,factor);
//                        stiffMatrix_jk.AddScal(massMatrix_jk,factor);
//    				}
//                }
//           }
//
//			//solve for trial state
//			NuTo::SparseMatrixCSRGeneral<double> hessianModSolver(stiffMatrix_jj);
////FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> hessianModFull(stiffMatrix_jj);
////std::cout << "stiffness\n" << hessianModFull << std::endl;
//
//			hessianModSolver.SetOneBasedIndexing();
//			mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
//			delta_disp_j*=-1;
//
//            //calculate trial state
//            disp_j += delta_disp_j;
//            disp_k = bRHSend - Cmat*disp_j;
//            if (this->IsDynamic())
//            {
//                acc_j  += delta_disp_j * (1./(timeStep*timeStep*mBeta));
//                vel_j  += delta_disp_j * (mGamma/(timeStep*mBeta));
//
//                acc_k = bRHSddot - (Cmat*acc_j);
//                vel_k = bRHSdot - (Cmat*vel_j);
//            }
//
//            //apply displacements
//            //std::cout << "norm of delta disp " << delta_disp_j.Norm() << std::endl;
//            mStructure->NodeMergeActiveDofValues(0,disp_j);
//            mStructure->ElementTotalUpdateTmpStaticData();
//
//            //calculate internal force, update of history variables=false
//            mStructure->BuildGlobalGradientInternalPotentialSubVectors(intForce_j, intForce_k,false);
//            residual_j = intForce_j;
//            residual_k = intForce_k;
//            if (this->IsDynamic())
//            {
//                if (mMuDampingMass>0)
//                {
//                    if (mUseLumpedMass)
//                    {
//                    	residual_j += (lumped_massMatrix_j.asDiagonal()*vel_j) *mMuDampingMass;
//                    }
//                    else
//                    {
//						//add damping terme
//						residual_j += (massMatrix_jj*vel_j+massMatrix_jk*vel_k) *mMuDampingMass;
//						residual_k += (massMatrix_kj*vel_j+massMatrix_kk*vel_k) *mMuDampingMass;
//                    }
//                }
//
//                if (mUseLumpedMass)
//                {
//					residual_j += lumped_massMatrix_j.asDiagonal()*acc_j;
//
//                }
//                else
//                {
//					//add mass terme
//					residual_j += massMatrix_jj*acc_j+massMatrix_jk*acc_k;
//					residual_k += massMatrix_kj*acc_j+massMatrix_kk*acc_k;
//                }
//            }
//
//            //add external force
//            residual_j -= extForce_j;
//            residual_k -= extForce_k;
//
//            //calculate the initial residual
//            if (Cmat.GetNumEntries()>0)
//            {
//                residual_mod=residual_j - CmatT*residual_k;
//            }
//            else
//            {
//                residual_mod=residual_j;
//            }
//
//            //calculate norm of initial residual (first time step)
//            double normResidual(residual_mod.Norm());
//            //std::cout << "max Residual " << residual_mod.Max() << std::endl;
//
//            int iteration(0);
//            while(normResidual>mToleranceForce && iteration<mMaxNumIterations)
//            {
//                //calculate stiffness
//                if (Cmat.GetNumEntries()>0)
//                {
//                    stiffMatrix_jj.SetZeroEntries();stiffMatrix_jk.SetZeroEntries();stiffMatrix_kj.SetZeroEntries();stiffMatrix_kk.SetZeroEntries();
//                    mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::STIFFNESS, stiffMatrix_jj, stiffMatrix_jk, stiffMatrix_kj, stiffMatrix_kk);
//                    if (this->IsDynamic())
//                    {
//                        //add damping and mass to full hessian
//                        double factor((1.+mGamma*mMuDampingMass*timeStep)/(timeStep*timeStep*mBeta));
//                        if (mUseLumpedMass)
//                        {
//                        	stiffMatrix_jj.AddScalDiag(lumped_massMatrix_j,factor);
//                        }
//                        else
//                        {
//                            stiffMatrix_jj.AddScal(massMatrix_jj,factor);
//                            stiffMatrix_jk.AddScal(massMatrix_jk,factor);
//                        }
//                    }
//
//                    stiffMatrix_jj -= CmatT * stiffMatrix_kj + stiffMatrix_jk * Cmat;
//                    stiffMatrix_jj += CmatT * stiffMatrix_kk * Cmat;
//                }
//                else
//                {
//                    stiffMatrix_jj.SetZeroEntries();stiffMatrix_jk.SetZeroEntries();
//                    mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::STIFFNESS, stiffMatrix_jj, stiffMatrix_jk);
//
//                    if (this->IsDynamic())
//                    {
//                        //add damping and mass to full hessian
//                        double factor((1.+mGamma*mMuDampingMass*timeStep)/(timeStep*timeStep*mBeta));
//                        if (mUseLumpedMass)
//                        {
//                        	stiffMatrix_jj.AddScalDiag(lumped_massMatrix_j,factor);
//                        }
//                        else
//                        {
//							stiffMatrix_jj.AddScal(massMatrix_jj,factor);
//							stiffMatrix_jk.AddScal(massMatrix_jk,factor);
//                        }
//                    }
//                }
//
//                //solve for new state
//                NuTo::SparseMatrixCSRGeneral<double> hessianModSolver(stiffMatrix_jj);
//                hessianModSolver.SetOneBasedIndexing();
//                mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
//                delta_disp_j*=-1;
//
//                //std::cout << "norm of delta disp " << delta_disp_j.Norm() <<  std::endl;
//                //std::cout << "delta disp " << std::endl;
//                //std::cout << delta_disp_j.Trans() << std::endl;
//
//                //perform a line search
//                double alpha(1);
//                double trialNormResidual(0);
//                do
//                {
//                     //calculate trial state
//                    trial_disp_j = disp_j + delta_disp_j * alpha;
//                    delta_disp_k = (Cmat*delta_disp_j)*(-1.);
//                    trial_disp_k = disp_k + delta_disp_k * alpha;
//
//                    if (this->IsDynamic())
//                    {
//                        trial_acc_j  = acc_j+delta_disp_j * (alpha/(timeStep*timeStep*mBeta));
//                        trial_acc_k  = acc_k+delta_disp_k * (alpha/(timeStep*timeStep*mBeta));
//
//                        trial_vel_j  = vel_j+delta_disp_j * (alpha*mGamma/(timeStep*mBeta));
//                        trial_vel_k  = vel_k+delta_disp_k * (alpha*mGamma/(timeStep*mBeta));
//                    }
//
//                    //apply displacements
//                    mStructure->NodeMergeActiveDofValues(0,trial_disp_j);
//                    mStructure->ElementTotalUpdateTmpStaticData();
//
//                    //calculate internal force, update of history variables=false
//                    mStructure->BuildGlobalGradientInternalPotentialSubVectors(intForce_j, intForce_k,false);
//                    residual_j = intForce_j;
//                    residual_k = intForce_k;
//
//                    if (this->IsDynamic())
//                    {
//                        if (mMuDampingMass>0)
//                        {
//                            if (mUseLumpedMass)
//                            {
//								residual_j += (lumped_massMatrix_j.asDiagonal()*trial_vel_j) *mMuDampingMass;
//								residual_k += (lumped_massMatrix_k.asDiagonal()*trial_vel_k) *mMuDampingMass;
//                            }
//                            else
//                            {
//								//add damping terme
//								residual_j += (massMatrix_jj*trial_vel_j+massMatrix_jk*trial_vel_k) *mMuDampingMass;
//								residual_k += (massMatrix_jk.TransMult(trial_vel_j)+massMatrix_kk*trial_vel_k) *mMuDampingMass;
//                            }
//                        }
//
//                        //add mass terme
//                        if (mUseLumpedMass)
//                        {
//							residual_j += lumped_massMatrix_j.asDiagonal()*trial_acc_j;
//							residual_k += lumped_massMatrix_k.asDiagonal()*trial_acc_k;
//                        }
//                        else
//                        {
//                            residual_j += (massMatrix_jj*trial_acc_j+massMatrix_jk*trial_acc_k);
//                            residual_k += (massMatrix_jk.TransMult(trial_acc_j)+massMatrix_kk*trial_acc_k);
//                        }
//                    }
//
//                    //add external force (assumed to be independent of the deformation)
//                    residual_j -=extForce_j;
//                    residual_k -=extForce_k;
//
//                    //std::cout << "residual_j\n " << residual_j << std::endl;
//                    //std::cout << "residual_k\n " << residual_k << std::endl;
//
//                    //calculate the initial residual
//                    if (Cmat.GetNumEntries()>0)
//                    {
//                        residual_mod = residual_j - CmatT*residual_k;
//                    }
//                    else
//                    {
//                        residual_mod = residual_j;
//                    }
//                    //std::cout << "residual_mod\n " << residual_mod << std::endl;
//
//                    trialNormResidual=residual_mod.Norm();
//
//                    std::cout << "  linesearch alpha " << alpha <<" previous residual " << normResidual << " trial residual " << trialNormResidual <<  "\n";
//
//                    alpha*=0.5;
//                }
//                while(alpha>mMinLineSearchStep && trialNormResidual>(1.-alpha)*normResidual);
//                   //in the first iteration, there is no line search
//
//                if (alpha>mMinLineSearchStep)
//                {
//                    //improvement is achieved, go to next Newton step
//                    disp_j = trial_disp_j;
//                    disp_k = trial_disp_k;
//
//                    if (this->IsDynamic())
//                    {
//                        acc_j  = trial_acc_j;
//                        acc_k  = trial_acc_k;
//
//                        vel_j  = trial_vel_j;
//                        vel_k  = trial_vel_k;
//                    }
//
//                    normResidual = trialNormResidual;
//
//                    iteration++;
//                }
//                else
//                {
//                    //and leave
//                    iteration = mMaxNumIterations;
//                }
//            }
//            // end of while(normResidual<mToleranceForce && iteration<mMaxNumIterations)
//            if (normResidual<=mToleranceForce)
//            {
//                //converged solution
//            	std::cout << " *** UpdateStaticData *** from JumpDirect" << std::endl;
//                mStructure->ElementTotalUpdateStaticData();
//
//                //store converged step
//                lastConverged_disp_j = disp_j;
//                lastConverged_disp_k = disp_k;
//                if (this->IsDynamic())
//                {
//                    lastConverged_vel_j = vel_j;
//                    lastConverged_vel_k = vel_k;
//                    lastConverged_acc_j = acc_j;
//                    lastConverged_acc_k = acc_k;
//                }
//
//                prevResidual_j = residual_j;
//                prevResidual_k = residual_k;
//
//                //update nodal data
//                mStructure->NodeMergeActiveDofValues(0,disp_j);
//                if (this->IsDynamic())
//                {
//                    mStructure->NodeMergeDofValues(1,vel_j,vel_k);
//                    mStructure->NodeMergeDofValues(2,acc_j,acc_k);
//                }
//                mStructure->ElementTotalUpdateTmpStaticData();
//
//                //update structure time
//                mStructure->SetPrevTime(curTime);
//
//                //postprocess
//                mTime+=timeStep;
//
//                std::cout << "Convergence after " << iteration << " iterations at time " << mTime << "(timestep " << timeStep << ").\n";
//                //std::cout << "plot Vector " << plotVector << std::endl;
//
//				//perform Postprocessing
//            	std::cout << " *** PostProcess *** from JumpDirect" << std::endl;
//				PostProcess(prevResidual_j, prevResidual_k);
//
//                //eventually increase next time step
//                if (mAutomaticTimeStepping && iteration<0.25*mMaxNumIterations)
//                {
//                    timeStep*=1.5;
//                    if (timeStep>mMaxTimeStep)
//                        timeStep = mMaxTimeStep;
//                }
//            }
//            else
//            {
//                std::cout << "No convergence with timestep " << timeStep << "\n";
//                //no convergence
//                if (mAutomaticTimeStepping)
//                {
//                    //no convergence, reduce the time step and start from scratch
//                    curTime -= timeStep;
//                    timeStep*=0.5;
//                }
//                else
//                {
//                    throw MechanicsException("[NuTo::JumpDirect::Solve] no convergence with the current maximum number of iterations, either use automatic time stepping, reduce the time step or the minimal line search cut back factor.");
//                }
//            }
//        }
//    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::JumpDirect::Solve] performing Newton-Raphson iteration.");
        throw e;
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mStructure->GetLogger()<<"[NuTo::JumpDirect::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mStructure->GetLogger()<< "[NuTo::JumpDirect::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return NuTo::Error::SUCCESSFUL;

}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::JumpDirect::GetTypeId()const
{
    return std::string("JumpDirect");
}

//!@brief calculate the global stiffness matrix for the BVP with mean displacement (rFourier = 0) or with displacement amplitude (rFourier = 1)
void NuTo::JumpDirect::CalculateGlobalModifiedStiffness(NuTo::SparseMatrixCSRVector2General<double>* rStiffnessMod, int rFourierMode)
{
	switch (rFourierMode)
	{
		case 0:

			break;
		case 1:

			break;
		default:
			throw MechanicsException("[NuTo::JumpDirect::CalculateModifiedStiffness] Implemented for monoharmonic excitation: 0 - mean displacement, 1 - displacement amplitude.");
	}

	//calculate matrices
    //calculate constraint matrix
    NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
    mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
    NuTo::SparseMatrixCSRVector2General<double> Cmat(CmatTmp);
    SparseMatrixCSRVector2General<double> CmatT (Cmat.Transpose());

    // allocate space for stiffness matrix
    SparseMatrixCSRVector2General<double> stiffMatrix_jj(mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
    SparseMatrixCSRVector2General<double> stiffMatrix_jk(mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
    SparseMatrixCSRVector2General<double> stiffMatrix_kj(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
    SparseMatrixCSRVector2General<double> stiffMatrix_kk(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
    mStructure->BuildGlobalElasticStiffnessSubMatricesGeneral(stiffMatrix_jj, stiffMatrix_jk, stiffMatrix_kj, stiffMatrix_kk);

    //calculate individual mass matrix
    NuTo::SparseMatrixCSRVector2General<double> massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk;

	massMatrix_jj.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
	massMatrix_jk.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
	massMatrix_kj.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
	massMatrix_kk.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
	mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);

	if (rFourierMode > 0 && mStructure->GetNumTimeDerivatives()>1)
		{
		double frequency(mHarmonicConstraintFactor(0,1));
		const double pi = boost::math::constants::pi<double>();
		double factor(2*rFourierMode*pi*frequency);

        stiffMatrix_jj.AddScal(massMatrix_jj,-factor*factor);
        stiffMatrix_jk.AddScal(massMatrix_jk,-factor*factor);
        stiffMatrix_kj.AddScal(massMatrix_kj,-factor*factor);
        stiffMatrix_kk.AddScal(massMatrix_kk,-factor*factor);
	}

	if (Cmat.GetNumEntries()>0) {
        stiffMatrix_jj -= CmatT * stiffMatrix_kj + stiffMatrix_jk * Cmat;
        stiffMatrix_jj += CmatT * stiffMatrix_kk * Cmat;
	}

	(*rStiffnessMod) = stiffMatrix_jj;

}

//!@brief find equilibrium by calculating the Fourier coefficients (which are the displacement fields)
void NuTo::JumpDirect::CalculateFourierCofficients(NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Mean_j, NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Mean_k,
		NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Ampl_j, NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Ampl_k)
{
    NuTo::FullVector<double,Eigen::Dynamic> extForce_j, extForce_k;
    NuTo::FullVector<double,Eigen::Dynamic> intForce_j(mStructure->GetNumActiveDofs()), intForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs()),
    		intForceMax_j(mStructure->GetNumActiveDofs()), intForceMax_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
    NuTo::FullVector<double,Eigen::Dynamic> residual_j, residual_k, residual_mod;
	NuTo::FullVector<double,Eigen::Dynamic> delta_disp_j, delta_disp_k;;
    FullVector<double,Eigen::Dynamic> bRHS, bRHSmax;

    //calculate constraint matrix
    NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
    mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
    NuTo::SparseMatrixCSRVector2General<double> Cmat(CmatTmp);
    SparseMatrixCSRVector2General<double> CmatT (Cmat.Transpose());

    // allocate space for stiffness matrix
    SparseMatrixCSRVector2General<double> stiffMatrix_jj(mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());


	double BeginHarmonicExcitation(mTimeDependentConstraintFactor(mTimeDependentConstraintFactor.GetNumRows()-1,0));

	// find the mean value (rFourier = 0)

	// set BC for the mean displacement, which is the displacement at the beginning of the harmonic excitation
	double timeDependentConstraintFactor(this->CalculateTimeDependentConstraintFactor(BeginHarmonicExcitation));
	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
	mStructure->ConstraintGetRHSAfterGaussElimination(bRHS);
	// set mean displacement field
	mStructure->NodeMergeActiveDofValues(0,*rDisp_Mean_j);
	mStructure->ElementTotalUpdateTmpStaticData();
	// calculate internal forces and external forces
	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(intForce_j,intForce_k);
	CalculateExternalLoad(*mStructure, BeginHarmonicExcitation, extForce_j, extForce_k);

	residual_j = intForce_j - extForce_j;
	residual_k = intForce_k - extForce_k;



	////////**************

    //calculate individual mass matrix
    NuTo::SparseMatrixCSRVector2General<double> massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk;
    NuTo::FullVector<double, Eigen::Dynamic> lumped_massMatrix_j(mStructure->GetNumDofs()),
    		                                 lumped_massMatrix_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
	//store last converged velocities and accelerations
	NuTo::FullVector<double,Eigen::Dynamic> lastConverged_vel_j, lastConverged_vel_k, lastConverged_acc_j, lastConverged_acc_k;

	bool FirstEqulibrium(false);
    if (mStructure->GetNumTimeDerivatives()>1 && FirstEqulibrium)
    {
    	FirstEqulibrium = false;
    	if (mUseLumpedMass)
    	{
    		mStructure->BuildGlobalLumpedHession2(lumped_massMatrix_j,lumped_massMatrix_k);
    		std::cout << "use lumped mass matrix " << std::endl;
    	}
    	else
    	{
			massMatrix_jj.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
			massMatrix_jk.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
			massMatrix_kj.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
			massMatrix_kk.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
			mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);
    		std::cout << "use full mass matrix " << std::endl;
    	}

        mStructure->NodeExtractDofValues(1,lastConverged_vel_j,lastConverged_vel_k);
        mStructure->NodeExtractDofValues(2,lastConverged_acc_j,lastConverged_acc_k);

        if (mMuDampingMass>0)
        {
            //add damping terms
        	if (mUseLumpedMass)
        	{
        		residual_j += lumped_massMatrix_j.asDiagonal()*lastConverged_vel_j * mMuDampingMass;
        		residual_k += lumped_massMatrix_k.asDiagonal()*lastConverged_vel_k * mMuDampingMass;
        	}
        	else
        	{
				residual_j += (massMatrix_jj*lastConverged_vel_j+massMatrix_jk*lastConverged_vel_k) *mMuDampingMass;
				residual_k += (massMatrix_kj*lastConverged_vel_j+massMatrix_kk*lastConverged_vel_k) *mMuDampingMass;
        	}
        }

    	if (mUseLumpedMass)
    	{
    		residual_j += lumped_massMatrix_j.asDiagonal()*lastConverged_acc_j;
    		residual_k += lumped_massMatrix_k.asDiagonal()*lastConverged_acc_k;
    	}
    	else
    	{
			//add mass terms
			residual_j += (massMatrix_jj*lastConverged_acc_j+massMatrix_jk*lastConverged_acc_k);
			residual_k += (massMatrix_kj*lastConverged_acc_j+massMatrix_kk*lastConverged_acc_k);
    	}

    }

	////////**************

	if (Cmat.GetNumEntries()>0)
	{
		residual_mod = residual_j - CmatT*residual_k;
	} else {
		residual_mod = residual_j;
	}

	this->CalculateGlobalModifiedStiffness(&stiffMatrix_jj,0);

	// solve for the mean displacement
	NuTo::SparseMatrixCSRGeneral<double> hessianModSolver(stiffMatrix_jj);

	hessianModSolver.SetOneBasedIndexing();
    // allocate solver
    NuTo::SparseDirectSolverMUMPS mySolver;
	mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
	delta_disp_j*=-1;

	// calculate mean displacement
	*rDisp_Mean_j += delta_disp_j;
	*rDisp_Mean_k = bRHS - Cmat*(*rDisp_Mean_j);

	////////**************
	// the previous calculation of disp_j is based on the update scheme
	// the following calculation is based on direct calculation, that is:
	// Kmod * rDisp_mean_j = -( (K12 - CmatT*K22)*bRHS - Fext ), Fint = 0
	bool DirectCalculation(false);
    // allocate space for stiffness matrix
    SparseMatrixCSRVector2General<double> K_jj(mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
    SparseMatrixCSRVector2General<double> K_jk(mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
    SparseMatrixCSRVector2General<double> K_kj(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
    SparseMatrixCSRVector2General<double> K_kk(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
    mStructure->BuildGlobalElasticStiffnessSubMatricesGeneral(K_jj, K_jk, K_kj, K_kk);
    if (DirectCalculation)
    {
    	// contribution of external forces to residual
    	residual_j = - extForce_j;
    	residual_k = - extForce_k;

    	// contribution of bRHS to residual
        residual_j += K_jk * bRHS;
        residual_k += K_kk * bRHS;

    	if (Cmat.GetNumEntries()>0)
    	{
    		residual_mod = residual_j - CmatT*residual_k;
    	} else {
    		residual_mod = residual_j;
    	}

    	// solve
    	mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
    	delta_disp_j*=-1;

    	// calculate mean displacement
    	*rDisp_Mean_j = delta_disp_j;
    	*rDisp_Mean_k = bRHS - Cmat*(*rDisp_Mean_j);

    	// calculate internal force
    	mStructure->NodeMergeActiveDofValues(0,*rDisp_Mean_j);
    	mStructure->ElementTotalUpdateTmpStaticData();
    	// calculate internal force
    	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(intForce_j,intForce_k);
    }

	////////**************

	// find the amplitude (rFourier = 1)

	// we need to calculate intForce(maximal displacement) - intForce(mean displacement)
	// calculate time at which the harmonic displacement is maximal, this is a quarter of the period
	double MaximumHarmonicExcitation(BeginHarmonicExcitation + 0.25/mHarmonicConstraintFactor(0,1));
	// set BC for the maximal displacement
	double timeDependentConstraintFactorMax(this->CalculateTimeDependentConstraintFactor(MaximumHarmonicExcitation));
	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactorMax);
	mStructure->ConstraintGetRHSAfterGaussElimination(bRHSmax);
	// set maximal displacement
	mStructure->NodeMergeActiveDofValues(0,(*rDisp_Mean_j) + (*rDisp_Ampl_j));
	mStructure->ElementTotalUpdateTmpStaticData();

	// calculate intForce(maximal displacement)
//	intForce_j.Zero(mStructure->GetNumActiveDofs()); intForce_k.Zero(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
//	intForceMax_j.Zero(mStructure->GetNumActiveDofs()); intForceMax_k.Zero(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(intForceMax_j,intForceMax_k);
	intForce_j = intForceMax_j - intForce_j; // intForce is currently the intForce(meanDisplacement)
	intForce_k = intForceMax_k - intForce_k;

//	double timeDependentConstraintFactorAmpl(this->CalculateTimeDependentConstraintFactor(MaximumHarmonicExcitation));
//	timeDependentConstraintFactorAmpl -= timeDependentConstraintFactor;
//	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactorAmpl);
//	bRHS.Zero(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
//	mStructure->ConstraintGetRHSAfterGaussElimination(bRHS);
	bRHS = bRHSmax - bRHS;
	// set ampl displacement field
//	mStructure->NodeMergeActiveDofValues(0,*rDisp_Ampl_j);
//	mStructure->ElementTotalUpdateTmpStaticData();
	// calculate internal forces, external forces harmonic amplitude should be provided here
	// currently for the displacement controlled excitation, the external force is included into the mean equilibrium only
//	intForce_j.Zero(mStructure->GetNumActiveDofs()); intForce_k.Zero(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
//	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(intForce_j,intForce_k);

	// calculate the residuals and stiffness
	residual_j.Zero(mStructure->GetNumActiveDofs()); residual_k.Zero(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());

	residual_j = intForce_j;
	residual_k = intForce_k;

    //calculate individual mass matrix
//    NuTo::SparseMatrixCSRVector2General<double> massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk;
//    NuTo::FullVector<double, Eigen::Dynamic> lumped_massMatrix_j(mStructure->GetNumDofs()),
//    		                                 lumped_massMatrix_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());

	if (mUseLumpedMass)
	{
		mStructure->BuildGlobalLumpedHession2(lumped_massMatrix_j,lumped_massMatrix_k);
		std::cout << "use lumped mass matrix " << std::endl;
	}
	else
	{
		massMatrix_jj.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
		massMatrix_jk.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
		massMatrix_kj.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
		massMatrix_kk.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
		mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);
		std::cout << "use full mass matrix " << std::endl;
	}

	double frequency(mHarmonicConstraintFactor(0,1));
	const double pi = boost::math::constants::pi<double>();
	double factor(2*pi*frequency);

	if (mStructure->GetNumTimeDerivatives()>1)
	{
		if (mUseLumpedMass)
		{
			residual_j -= lumped_massMatrix_j.asDiagonal()*(*rDisp_Ampl_j)*factor*factor;
			residual_k -= lumped_massMatrix_k.asDiagonal()*(*rDisp_Ampl_k)*factor*factor;
		}
		else
		{
			//add mass terms
			residual_j -= (massMatrix_jj*(*rDisp_Ampl_j)+massMatrix_jk*(*rDisp_Ampl_k))*factor*factor;
			residual_k -= (massMatrix_kj*(*rDisp_Ampl_j)+massMatrix_kk*(*rDisp_Ampl_k))*factor*factor;
		}
	}

	if (Cmat.GetNumEntries()>0)
	{
		residual_mod = residual_j - CmatT*residual_k;
	} else {
		residual_mod = residual_j;
	}

	stiffMatrix_jj.SetZeroEntries();

	this->CalculateGlobalModifiedStiffness(&stiffMatrix_jj,1);

	// solve for the ampl displacement
	hessianModSolver = stiffMatrix_jj;

	hessianModSolver.SetOneBasedIndexing();
	delta_disp_j.Zero(mStructure->GetNumActiveDofs());
	mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
	delta_disp_j*=-1;

	// calculate ampl displacement
	*rDisp_Ampl_j += delta_disp_j;
	*rDisp_Ampl_k = bRHS - Cmat*(*rDisp_Ampl_j);

	////////**************
	// the previous calculation of rDisp_Ampl_j is based on the update scheme
	// the following calculation is based on direct calculation, that is:
	// Kmod * rDisp_mean_j = -( (K12 - CmatT*K22)*delta_bRHS ), Fint = 0
	DirectCalculation = true;
    if (DirectCalculation)
    {
    	// since the external load is constant, there is no contribution to the Fourier amplitude of external forces

        K_jj.AddScal(massMatrix_jj,-factor*factor);
        K_jk.AddScal(massMatrix_jk,-factor*factor);
        K_kj.AddScal(massMatrix_kj,-factor*factor);
        K_kk.AddScal(massMatrix_kk,-factor*factor);
    	// contribution of bRHS to residual
        residual_j = K_jk * bRHS;
        residual_k = K_kk * bRHS;

    	if (Cmat.GetNumEntries()>0)
    	{
    		residual_mod = residual_j - CmatT*residual_k;
    	} else {
    		residual_mod = residual_j;
    	}

    	// solve
    	mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
    	delta_disp_j*=-1;

    	// calculate mean displacement
    	*rDisp_Ampl_j = delta_disp_j;
    	*rDisp_Ampl_k = bRHS - Cmat*(*rDisp_Ampl_j);
    }

	////////**************
}

//!@brief straight-forward integration of a single cycle with a prescribed Fourier coefficients
// the displacement fields have the same meaning as above, see CalculateFourierCoefficients
// rIncludePostProcess ...false, if no postprocessing should be done during integration
// rIncludePostProcess ...true, postprocessing will be done
// Postprocessing should be done if the jump is acceptable; in this case the IntegrateSinleCycle should be repeated with a true option
void NuTo::JumpDirect::IntegrateSingleCycle(NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Mean_j, NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Mean_k,
		NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Ampl_j, NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Ampl_k, bool rIncludePostProcess)
{
    //calculate constraint matrix
    NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
    mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
    NuTo::SparseMatrixCSRVector2General<double> Cmat(CmatTmp);
    SparseMatrixCSRVector2General<double> CmatT (Cmat.Transpose());
    FullVector<double,Eigen::Dynamic> bRHSprev, bRHShalf, bRHSend, bRHSdot, bRHSddot;

    // allocate space for stiffness matrix
    SparseMatrixCSRVector2General<double> stiffMatrix_jj(mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
    SparseMatrixCSRVector2General<double> stiffMatrix_jk(mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
    SparseMatrixCSRVector2General<double> stiffMatrix_kj(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
    SparseMatrixCSRVector2General<double> stiffMatrix_kk(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());

    //calculate individual mass matrix (I was lazy here and have just used general matrices, but symmetric is certainly better here)
    NuTo::SparseMatrixCSRVector2General<double> massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk;
    NuTo::FullVector<double, Eigen::Dynamic> lumped_massMatrix_j(mStructure->GetNumDofs()),
    		                                 lumped_massMatrix_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());

    // check equilibrium
    if (mStructure->GetNumTimeDerivatives()>1)
    {
    	if (mUseLumpedMass)
    	{
    		mStructure->BuildGlobalLumpedHession2(lumped_massMatrix_j,lumped_massMatrix_k);
    		std::cout << "use lumped mass matrix " << std::endl;
    	}
    	else
    	{
			massMatrix_jj.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
			massMatrix_jk.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
			massMatrix_kj.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
			massMatrix_kk.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
			mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);
    		std::cout << "use full mass matrix " << std::endl;
    	}
    }

    //extract displacements, velocities and accelerations
    NuTo::FullVector<double,Eigen::Dynamic> disp_j,disp_k, vel_j, vel_k, acc_j, acc_k, delta_disp_j, delta_disp_k;
    NuTo::FullVector<double,Eigen::Dynamic> trial_disp_j,trial_disp_k, trial_vel_j, trial_vel_k, trial_acc_j, trial_acc_k;
    NuTo::FullVector<double,Eigen::Dynamic> extForce_j, extForce_k, prevExtForce_j, prevExtForce_k;
    NuTo::FullVector<double,Eigen::Dynamic> intForce_j(mStructure->GetNumActiveDofs()), intForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs()),
                             prevIntForce_j(mStructure->GetNumActiveDofs()), prevIntForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
    NuTo::FullVector<double,Eigen::Dynamic> residual_j, residual_k, residual_mod, prevResidual_j, prevResidual_k;

    //store last converged displacements, velocities and accelerations
    NuTo::FullVector<double,Eigen::Dynamic> lastConverged_disp_j,lastConverged_disp_k, lastConverged_vel_j, lastConverged_vel_k, lastConverged_acc_j, lastConverged_acc_k;
//  mStructure->NodeExtractDofValues(0,lastConverged_disp_j, lastConverged_disp_k);

    if (mStructure->GetNumTimeDerivatives()>1)
    {
        mStructure->NodeExtractDofValues(1,lastConverged_vel_j,lastConverged_vel_k);
        mStructure->NodeExtractDofValues(2,lastConverged_acc_j,lastConverged_acc_k);
    }

    double curTime(mStructure->GetTime());
    double startingCycleTime(curTime);
    // initialize the structure previous times
    mStructure->SetPrevTime(curTime);

    //apply constraints for last converged time step
    double timeDependentConstraintFactor(0);
    if (mTimeDependentConstraint!=-1)
    {
        timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
        mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
    }

    lastConverged_disp_j = *rDisp_Mean_j;

    mStructure->ConstraintGetRHSAfterGaussElimination(bRHSprev);
    mStructure->NodeMergeActiveDofValues(0,lastConverged_disp_j); //disp_k is internally calculated from the previous constraint matrix
    mStructure->ElementTotalUpdateTmpStaticData();

    //calculate internal force, update of history variables=false
    mStructure->BuildGlobalGradientInternalPotentialSubVectors(prevIntForce_j,prevIntForce_k,false);
    prevResidual_j = prevIntForce_j;
    prevResidual_k = prevIntForce_k;
    std::cout << "     RES.: prevIntForce_j = " << prevIntForce_j.Norm() << std::endl;

    if (mStructure->GetNumTimeDerivatives()>1)
    {
        if (mMuDampingMass>0)
        {
            //add damping terms
        	if (mUseLumpedMass)
        	{
        		prevResidual_j += lumped_massMatrix_j.asDiagonal()*lastConverged_vel_j * mMuDampingMass;
        		prevResidual_k += lumped_massMatrix_k.asDiagonal()*lastConverged_vel_k * mMuDampingMass;
        	}
        	else
        	{
				prevResidual_j += (massMatrix_jj*lastConverged_vel_j+massMatrix_jk*lastConverged_vel_k) *mMuDampingMass;
				prevResidual_k += (massMatrix_kj*lastConverged_vel_j+massMatrix_kk*lastConverged_vel_k) *mMuDampingMass;
        	}
        }

    	if (mUseLumpedMass)
    	{
    		prevResidual_j += lumped_massMatrix_j.asDiagonal()*lastConverged_acc_j;
    		prevResidual_k += lumped_massMatrix_k.asDiagonal()*lastConverged_acc_k;
    	}
    	else
    	{
			//add mass terms
			prevResidual_j += (massMatrix_jj*lastConverged_acc_j+massMatrix_jk*lastConverged_acc_k);
			prevResidual_k += (massMatrix_kj*lastConverged_acc_j+massMatrix_kk*lastConverged_acc_k);
		    std::cout << "     RES.: prevIntForce_j + M*acc_j = " << prevResidual_j.Norm() << std::endl;
    	}
    }
    //add external force
    CalculateExternalLoad(*mStructure, curTime, prevExtForce_j, prevExtForce_k);
    prevResidual_j -= prevExtForce_j;
    prevResidual_k -= prevExtForce_k;
    std::cout << "     RES.: prevExtForce_j = " << prevExtForce_j.Norm() << std::endl;


    if (Cmat.GetNumEntries()>0)
    {
        residual_mod=prevResidual_j - CmatT*prevResidual_k;
    }
    else
    {
        residual_mod=prevResidual_j;
    }

    if (mStructure->GetNumTimeDerivatives()>1) {
    	std::cout << "=== vel_j = " << lastConverged_vel_j[0] << std::endl;
    }

//    if (residual_mod.Norm()>mToleranceForce){
//        std::cout << "residual in initial configuration " << residual_mod.Norm() << std::endl;
//        throw MechanicsException("[NuTo::JumpDirect::Solve] Configuration after the extrapolation jump is not in (dynamic) equilibrium.");
//    }

    std::ofstream DamageFile;
    DamageFile.open("Damage.txt", std::ios::app);

    if (rIncludePostProcess){
    	mTime = curTime;
    	NuTo::NewmarkDirect::PostProcess(prevResidual_j, prevResidual_k);
    	DamageFile << mTime << " " << mStructure->ElementGetElementPtr(141)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl; // for Brick8N was ElementGetElementPtr(9)
	}

    // iterate over the increments within the cycle
	const double pi = boost::math::constants::pi<double>();
    double frequency(mHarmonicConstraintFactor(0,1));
    double timeStep(1./(mHarmonicIncrementation*frequency));
    for (int incr = 1; incr <= mHarmonicIncrementation; ++incr) {
    	// set new time
    	curTime += timeStep;
    	mStructure->SetTime(curTime);

    	// set BC
    	timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
    	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
    	bRHSend.Zero(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
    	mStructure->ConstraintGetRHSAfterGaussElimination(bRHSend);

		// add harmonic excitation
		disp_j = (*rDisp_Mean_j) + (*rDisp_Ampl_j)*sin(2*pi*frequency*(curTime - startingCycleTime));
		mStructure->NodeMergeActiveDofValues(0,disp_j);
		mStructure->ElementTotalUpdateTmpStaticData();

		// IP integration
//		mStructure->ElementTotalUpdateStaticData();
	    //calculate internal force, update of history variables=true
	    mStructure->BuildGlobalGradientInternalPotentialSubVectors(prevIntForce_j,prevIntForce_k,true);
	    prevResidual_j = prevIntForce_j;
	    prevResidual_k = prevIntForce_k;

	    //add external force
	    CalculateExternalLoad(*mStructure, curTime, prevExtForce_j, prevExtForce_k);
	    prevResidual_j -= prevExtForce_j;
	    prevResidual_k -= prevExtForce_k;
	    std::cout << "     RES.: prevIntForce_j = " << prevIntForce_j.Norm() << std::endl;

        // postprocessing
        if (rIncludePostProcess) {
        	// set postprocessing time
        	mTime = curTime;

			// calculate residual forces

            // calculate the velocity and acceleration of the rhs
            if (mStructure->GetNumTimeDerivatives()>1)
            {
                if (mTimeDependentConstraint!=-1)
                {
					timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime - 0.5*timeStep);
					mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
                }
                mStructure->ConstraintGetRHSAfterGaussElimination(bRHShalf);

                // get last velocities
                mStructure->NodeExtractDofValues(1,lastConverged_vel_j,lastConverged_vel_k);
                mStructure->NodeExtractDofValues(2,lastConverged_acc_j,lastConverged_acc_k);

                // calculate approximations to the time derivates of the rhs of the constraint matrix
                bRHSdot =  (bRHSend - bRHSprev)*(1./(timeStep));
                bRHSddot = (bRHSprev-bRHShalf*2+bRHSend)*(4./(timeStep*timeStep));

                // calculate new accelerations and velocities of independent dofs
                acc_j = (lastConverged_vel_j+lastConverged_acc_j*(timeStep*(0.5-mBeta)))*(-1./(timeStep*mBeta));
                acc_j += (disp_j - lastConverged_disp_j) * (1./(timeStep*timeStep*mBeta));
                vel_j = lastConverged_vel_j+lastConverged_acc_j*((1.-mGamma)*timeStep)+acc_j*(mGamma*timeStep);

                // calculate new accelerations and velocities of dependent dofs
                acc_k = bRHSddot - (Cmat*acc_j);
                vel_k = bRHSdot - (Cmat*vel_j);

                // calculate the dynamic contribution to the residuals
                if (mMuDampingMass>0)
                {
                	//add damping terms
                    if (mUseLumpedMass)
                    {
                    	prevResidual_j += lumped_massMatrix_j.asDiagonal()*lastConverged_vel_j * mMuDampingMass;
                    	prevResidual_k += lumped_massMatrix_k.asDiagonal()*lastConverged_vel_k * mMuDampingMass;
                    }
                    else
                    {
                		prevResidual_j += (massMatrix_jj*lastConverged_vel_j+massMatrix_jk*lastConverged_vel_k) *mMuDampingMass;
                		prevResidual_k += (massMatrix_kj*lastConverged_vel_j+massMatrix_kk*lastConverged_vel_k) *mMuDampingMass;
                    }
                }
                if (mUseLumpedMass)
                {
                	prevResidual_j += lumped_massMatrix_j.asDiagonal()*lastConverged_acc_j;
                	prevResidual_k += lumped_massMatrix_k.asDiagonal()*lastConverged_acc_k;
                }
                else
                {
                //add mass terms
                	prevResidual_j += (massMatrix_jj*lastConverged_acc_j+massMatrix_jk*lastConverged_acc_k);
                	prevResidual_k += (massMatrix_kj*lastConverged_acc_j+massMatrix_kk*lastConverged_acc_k);
                    std::cout << "     RES.: prevIntForce_j + M*acc_j = " << prevResidual_j.Norm() << std::endl;
                }

            }

        	std::cout << " Cyclic increment " << incr << std::endl;
        	NuTo::NewmarkDirect::PostProcess(prevResidual_j, prevResidual_k);
        	DamageFile << mTime << " " << mStructure->ElementGetElementPtr(141)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl; // for Brick8N was ElementGetElementPtr(9)

            // update previous BC and displacements
        	bRHSprev = bRHSend;
        	lastConverged_disp_j = disp_j;

        	// update nodal kinematic data
            if (mStructure->GetNumTimeDerivatives()>1)
            {
                mStructure->NodeMergeDofValues(1,vel_j,vel_k);
                mStructure->NodeMergeDofValues(2,acc_j,acc_k);
            }

		}

        // update structure time
        mStructure->SetPrevTime(curTime);
	}
    DamageFile.close();
}


#ifdef ENABLE_SERIALIZATION
//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
//! @brief ... save the object to a file
void NuTo::JumpDirect::Restore (const std::string &filename, std::string rType )
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        std::string tmpString;
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[JumpDirect::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[JumpDirect::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[JumpDirect::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            ota & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else
        {
            throw MathException ( "[Matrix::Restore]File type not implemented" );
        }
    }
    catch ( MechanicsException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
    catch ( ... )
    {
        throw MechanicsException ( "[JumpDirect::Restore]Unhandled exception." );
    }
}

//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::JumpDirect::Save (const std::string &filename, std::string rType )const
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        std::string tmpStr ( GetTypeId() );
        std::string baseClassStr = tmpStr.substr ( 4,100 );
        if (rType=="BINARY")
        {
            boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", tmpStr );
            oba & boost::serialization::make_nvp(tmpStr.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpStr );
            oxa & boost::serialization::make_nvp(tmpStr.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpStr );
            ota & boost::serialization::make_nvp(tmpStr.c_str(), *this);
        }
        else
        {
            throw MechanicsException ( "[JumpDirect::Save]File type not implemented." );
        }
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[JumpDirect::Save]File save exception in boost - " ) +std::string ( e.what() ) );
        std::cout << s << "\n";
        throw MathException ( s );
    }
    catch ( MechanicsException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
    catch ( ... )
    {
        throw MechanicsException ( "[JumpDirect::Save]Unhandled exception." );
    }
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::JumpDirect)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
