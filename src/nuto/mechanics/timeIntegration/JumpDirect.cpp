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
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage2DFatigue.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::JumpDirect::JumpDirect (StructureBase* rStructure)  : NewmarkDirect (rStructure)
{
    mMinLineSearchStep = 0.01;
    mHarmonicIncrementation = 16;
    mHarmonicExtrapolation = false;
    mHarmonicExtrapolationTolerance = 0.25;
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
    	this->SetNewmarkBeta(3./10.);
    	this->SetNewmarkGamma(11./20.);

    	//this->SetNewmarkBeta(0.5);   //for displac. controlled

    	//this->SetNewmarkBeta(15.);   //for displac. controlled
    	//this->SetNewmarkBeta(10.);   //for load control, 10%  incr 32
    	//this->SetNewmarkBeta(15.);   //for load control, 5% incr 32
    	//this->SetNewmarkBeta(50.);   //for load control, 10% incr 24

    	//calculate the end time of monotonic loading which is exactly the beginning of cyclic loading
    	NuTo::Error::eError Error;

    	//straightforward integration over the loading history, if extrapolation is not required
    	if (!mHarmonicExtrapolation) {
    		return NuTo::NewmarkDirect::Solve(rTimeDelta);
		}

    	//extrapolation of the harmonic response
    	double BeginHarmonicExcitation;
    	if (mTimeDependentConstraint != -1) {
    		BeginHarmonicExcitation = mTimeDependentConstraintFactor(mTimeDependentConstraintFactor.GetNumRows()-1,0);
    	} else if (mTimeDependentLoadCase != -1) {
    		BeginHarmonicExcitation = mTimeDependentLoadFactor(mTimeDependentLoadFactor.GetNumRows()-1,0);
    	} else {
    		throw MechanicsException("[NuTo::JumpDirect::Solve] the harmonic excitation can be only applied to the time dependent constraint or load case.");
    	}

    	//integrate till the end of the third cycle
    	Error = NuTo::NewmarkDirect::Solve(BeginHarmonicExcitation + 3./mHarmonicFactor(0,1));
    	if (Error != NuTo::Error::SUCCESSFUL) {
    		return Error;
    	}

    	std::cout << " * * * After Newmark " << std::endl;

    	//calculate matrices
        //calculate constraint matrix

        /*
        // ConstraintGetConstraintMatrixAfterGaussElimination(const NuTo::SparseMatrixCSRGeneral<double>& ... ) replaced by
        // const NuTo::SparseMatrixCSRGeneral<double> ConstraintGetConstraintMatrixAfterGaussElimination() const;
        // ---> No need for CmatTmp anymore! --- vhirtham

        NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
        mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
        */
        NuTo::SparseMatrixCSRVector2General<double> Cmat(mStructure->ConstraintGetConstraintMatrixAfterGaussElimination());
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
        mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);

        NuTo::FullVector<double,Eigen::Dynamic> extForce_j, extForce_k, prevExtForce_j, prevExtForce_k;
        NuTo::FullVector<double,Eigen::Dynamic> cyclicPreJumpIntForce_j(mStructure->GetNumActiveDofs()), cyclicPreJumpIntForceNext_j(mStructure->GetNumActiveDofs()),
        		cyclicPreJumpIntForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs()), cyclicPreJumpIntForceNext_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
        NuTo::FullVector<double,Eigen::Dynamic> cyclicAfterJumpIntForce_j(mStructure->GetNumActiveDofs()),
        		cyclicAfterJumpIntForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
        NuTo::FullVector<double,Eigen::Dynamic> cyclicPreJumpIntForce0_j(mStructure->GetNumActiveDofs()),
        		cyclicPreJumpIntForce0_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
        NuTo::FullVector<double,Eigen::Dynamic> cyclicPreJumpIntForceMax_j(mStructure->GetNumActiveDofs()),
                		cyclicPreJumpIntForceMax_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
        NuTo::FullVector<double,Eigen::Dynamic> residual_j, residual_k, residual_mod, prevResidual_j, prevResidual_k;
        NuTo::FullVector<double,Eigen::Dynamic> vel_j, vel_k, acc_j, acc_k;
        //preJump values are stored for situation if the jump has to be repeated
        NuTo::FullVector<double,Eigen::Dynamic> preJump_disp_Mean_j,preJump_disp_Mean_k, preJump_disp_Ampl_j, preJump_disp_Ampl_k, preJump_vel_j, preJump_vel_k,
				preJump_acc_j, preJump_acc_k;

        bool JumpCriterion;
        double JumpCriterionDouble;
        const double pi = boost::math::constants::pi<double>();
        const int NjumpMin(2);

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

        //save static data, time and displacements after the 3th cycle for situation if the jump has to be repeated
        NuTo::FullVector<double,Eigen::Dynamic> save_disp_j, save_disp_k, save_vel_j, save_vel_k, save_acc_j, save_acc_k;
    	double SaveTime(mStructure->GetTime());
    	mStructure->ElementFatigueSaveStaticData();
//    	std::cout << "Saved statevs =======================" <<std::endl;		// delete me
//    	if (mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE){
//    	 		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    	 		std::cout << "Damage  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << std::endl;
//    	 		std::cout << "DamageF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmegaFatigue() << std::endl;
//    	 		std::cout << "Kappa   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << std::endl;
//    	 		std::cout << "KappaF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappaFatigue() << std::endl;
//    	 		std::cout << "NonLEq  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() << std::endl;
//    	 		std::cout << "NonLEqF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrainFatigue() << std::endl;
//    	 		std::cout << "PrevEp  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrain() << std::endl;
//    	 		std::cout << "PrevEpF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrainFatigue() << std::endl;
//    	 		std::cout << "PrevS   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStress() << std::endl;
//    	 		std::cout << "PrevSF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStressFatigue() << std::endl;
//    	}
    	mStructure->NodeExtractDofValues(0,save_disp_j, save_disp_k);
        if (mStructure->GetNumTimeDerivatives()>1)
	    {
	        mStructure->NodeExtractDofValues(1,save_vel_j,save_vel_k);
	        mStructure->NodeExtractDofValues(2,save_acc_j,save_acc_k);
	    }

//      //save static data and time after the 3th cycle for situation if the jump has to be repeated
//    	NuTo::FullVector<double,Eigen::Dynamic> save_disp_j, save_disp_k, delta_disp_j, delta_disp_k;;
//    	double SaveTime(mStructure->GetTime());
//    	mStructure->NodeExtractDofValues(0,save_disp_j, save_disp_k);
//    	mStructure->ElementFatigueSaveStaticData();
//    	//calculate 10 cycles
//    	std::cout << "Calculate 10 Cycles from here, mTime = " << mTime << std::endl;
//    	Error = NuTo::NewmarkDirect::Solve(SaveTime + 10./mHarmonicFactor(0,1));

//    	//Example get damage
//    	if (mStructure->ElementGetElementPtr(0)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS){
//    		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(0)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    		std::cout << "Damage after 10 cycles = " << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl;
//    	}

//    	//repeat with restoring
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
//    	Error = NuTo::NewmarkDirect::Solve(SaveTime + 10./mHarmonicFactor(0,1));

    	// perform the 4th cycle
    	// initialize mean displacement and displacement amplitude
    	NuTo::FullVector<double,Eigen::Dynamic> disp_Max_j, disp_Max_k, disp_Min_j, disp_Min_k;

    	Error = NuTo::NewmarkDirect::Solve(SaveTime + 0.25/mHarmonicFactor(0,1));
//    	std::cout << "0.25 READY ==================" <<std::endl;		// delete me
    	mStructure->NodeExtractDofValues(0,disp_Max_j, disp_Max_k);
        if (mStructure->GetNumTimeDerivatives()>1)
	    {
	        mStructure->NodeExtractDofValues(1,preJump_vel_j,preJump_vel_k);
	        mStructure->NodeExtractDofValues(2,preJump_acc_j,preJump_acc_k);
//			std::cout << " i | PreJumpAccMax_k " << std::endl;
//			for (int i = 0; i < mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(); ++i) {
//				std::cout << "i = " << i << " | " << preJump_acc_k[i] << std::endl;
//			}
	    }

    	Error = NuTo::NewmarkDirect::Solve(SaveTime + 0.75/mHarmonicFactor(0,1));
//    	std::cout << "0.75 READY ==================" <<std::endl;		// delete me
    	mStructure->NodeExtractDofValues(0,disp_Min_j, disp_Min_k);

    	NuTo::FullVector<double,Eigen::Dynamic> disp_Mean_j, disp_Mean_k, disp_Ampl_j, disp_Ampl_k;

    	Error = NuTo::NewmarkDirect::Solve(SaveTime + 1./mHarmonicFactor(0,1));
//    	std::cout << "1 READY ==================" <<std::endl;		// delete me
//    	if (mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE){
//    	 		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    	 		std::cout << "Damage  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << std::endl;
//    	 		std::cout << "DamageF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmegaFatigue() << std::endl;
//    	 		std::cout << "Kappa   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << std::endl;
//    	 		std::cout << "KappaF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappaFatigue() << std::endl;
//    	 		std::cout << "NonLEq  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() << std::endl;
//    	 		std::cout << "NonLEqF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrainFatigue() << std::endl;
//    	 		std::cout << "PrevEp  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrain() << std::endl;
//    	 		std::cout << "PrevEpF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrainFatigue() << std::endl;
//    	 		std::cout << "PrevS   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStress() << std::endl;
//    	 		std::cout << "PrevSF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStressFatigue() << std::endl;
//    	}

    	mStructure->NodeExtractDofValues(0,disp_Mean_j, disp_Mean_k);

    	disp_Ampl_j = 0.5*(disp_Max_j - disp_Min_j);
    	disp_Ampl_k = 0.5*(disp_Max_k - disp_Min_k);

    	//save node data for situation, if the jump has to be repeated
        if (mStructure->GetNumTimeDerivatives()>1)
	    {
			mStructure->BuildGlobalGradientInternalPotentialSubVectors(cyclicPreJumpIntForce0_j,cyclicPreJumpIntForce0_k,false);

	        mStructure->NodeExtractDofValues(1,preJump_vel_j,preJump_vel_k);
	        mStructure->NodeExtractDofValues(2,preJump_acc_j,preJump_acc_k);
//			std::cout << "=== acc_j = " << preJump_acc_j[0] << std::endl;
	    }
        preJump_disp_Mean_j = disp_Mean_j; preJump_disp_Mean_k = disp_Mean_k;
        preJump_disp_Ampl_j = disp_Ampl_j; preJump_disp_Ampl_k = disp_Ampl_k;

        //evaluate cyclic change of the out-of-balance intForce => integrate the subsequent single cycles (ISC)
        //this is necessary for initialization of Njump
        //no postprocessing, therefore=>false
        //when calculating the internal force, update of history variables=false
        this->IntegrateSingleCycle(&disp_Mean_j,&disp_Mean_k,&disp_Ampl_j,&disp_Ampl_k,false);						    	//the statevs and structure time evolve through the cycle
//    	std::cout << "IntSingCycleReady ==================" <<std::endl;		// delete me
//    	if (mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE){
//    	 		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    	 		std::cout << "Damage  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << std::endl;
//    	 		std::cout << "DamageF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmegaFatigue() << std::endl;
//    	 		std::cout << "Kappa   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << std::endl;
//    	 		std::cout << "KappaF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappaFatigue() << std::endl;
//    	 		std::cout << "NonLEq  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() << std::endl;
//    	 		std::cout << "NonLEqF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrainFatigue() << std::endl;
//    	 		std::cout << "PrevEp  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrain() << std::endl;
//    	 		std::cout << "PrevEpF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrainFatigue() << std::endl;
//    	 		std::cout << "PrevS   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStress() << std::endl;
//    	 		std::cout << "PrevSF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStressFatigue() << std::endl;
//    	}
        mStructure->BuildGlobalGradientInternalPotentialSubVectors(cyclicPreJumpIntForce_j,cyclicPreJumpIntForce_k,false);	//calculate cyclic change of the internal forces within the 5th cycle
//    	std::cout << "Build Global Gradient Int Pot READY ==================" <<std::endl;		// delete me
//    	if (mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE){
//    	 		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    	 		std::cout << "Damage  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << std::endl;
//    	 		std::cout << "DamageF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmegaFatigue() << std::endl;
//    	 		std::cout << "Kappa   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << std::endl;
//    	 		std::cout << "KappaF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappaFatigue() << std::endl;
//    	 		std::cout << "NonLEq  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() << std::endl;
//    	 		std::cout << "NonLEqF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrainFatigue() << std::endl;
//    	 		std::cout << "PrevEp  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrain() << std::endl;
//    	 		std::cout << "PrevEpF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrainFatigue() << std::endl;
//    	 		std::cout << "PrevS   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStress() << std::endl;
//    	 		std::cout << "PrevSF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStressFatigue() << std::endl;
//    	}
        this->IntegrateSingleCycle(&disp_Mean_j,&disp_Mean_k,&disp_Ampl_j,&disp_Ampl_k,false);
//    	std::cout << "IntSingCycleReady2 ==================" <<std::endl;		// delete me
//    	if (mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE){
//    	 		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    	 		std::cout << "Damage  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << std::endl;
//    	 		std::cout << "DamageF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmegaFatigue() << std::endl;
//    	 		std::cout << "Kappa   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << std::endl;
//    	 		std::cout << "KappaF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappaFatigue() << std::endl;
//    	 		std::cout << "NonLEq  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() << std::endl;
//    	 		std::cout << "NonLEqF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrainFatigue() << std::endl;
//    	 		std::cout << "PrevEp  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrain() << std::endl;
//    	 		std::cout << "PrevEpF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrainFatigue() << std::endl;
//    	 		std::cout << "PrevS   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStress() << std::endl;
//    	 		std::cout << "PrevSF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStressFatigue() << std::endl;
//    	}
        mStructure->BuildGlobalGradientInternalPotentialSubVectors(cyclicPreJumpIntForceNext_j,cyclicPreJumpIntForceNext_k,false);	//calculate cyclic change within the 6th cycle
//    	std::cout << "Build Global Gradient Int Pot READY2 ==================" <<std::endl;		// delete me
//    	if (mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE){
//    	 		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    	 		std::cout << "Damage  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << std::endl;
//    	 		std::cout << "DamageF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmegaFatigue() << std::endl;
//    	 		std::cout << "Kappa   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << std::endl;
//    	 		std::cout << "KappaF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappaFatigue() << std::endl;
//    	 		std::cout << "NonLEq  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() << std::endl;
//    	 		std::cout << "NonLEqF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrainFatigue() << std::endl;
//    	 		std::cout << "PrevEp  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrain() << std::endl;
//    	 		std::cout << "PrevEpF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrainFatigue() << std::endl;
//    	 		std::cout << "PrevS   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStress() << std::endl;
//    	 		std::cout << "PrevSF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStressFatigue() << std::endl;
//    	}
        cyclicPreJumpIntForceNext_j -= cyclicPreJumpIntForce_j;
        cyclicPreJumpIntForceNext_k -= cyclicPreJumpIntForce_k;

        if (mStructure->GetNumTimeDerivatives()>1)
	    {
			cyclicPreJumpIntForce_j -= cyclicPreJumpIntForce0_j;
			cyclicPreJumpIntForce_k -= cyclicPreJumpIntForce0_k;
	    }
//        for (int i = 0; i <= sizeof(cyclicPreJumpIntForce0_j); ++i){
//        	std::cout << "i = " << i << " , cyclicPreJump0 = " << cyclicPreJumpIntForce0_j[i] << std::endl;
//        }

        //restore the structure after the 3th cycle and integrate the 4th cycle
        //this is necessary for extrapolation of statevs
        mStructure->SetPrevTime(SaveTime);
        mStructure->SetTime(SaveTime);
        mTime = SaveTime; 																								//postprocessing time
//    	std::cout << "Before starting restoring =======================" <<std::endl;		// delete me
//    	if (mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE){
//    	 		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    	 		std::cout << "Damage  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << std::endl;
//    	 		std::cout << "DamageF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmegaFatigue() << std::endl;
//    	 		std::cout << "Kappa   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << std::endl;
//    	 		std::cout << "KappaF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappaFatigue() << std::endl;
//    	 		std::cout << "NonLEq  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() << std::endl;
//    	 		std::cout << "NonLEqF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrainFatigue() << std::endl;
//    	 		std::cout << "PrevEp  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrain() << std::endl;
//    	 		std::cout << "PrevEpF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrainFatigue() << std::endl;
//    	 		std::cout << "PrevS   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStress() << std::endl;
//    	 		std::cout << "PrevSF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStressFatigue() << std::endl;
//    	}
        mStructure->ElementFatigueRestoreStaticData();
//    	std::cout << "Restored to the Third cycle =======================" <<std::endl;		// delete me
//    	if (mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE){
//    	 		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() << std::endl;
//    	 		std::cout << "Damage  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << std::endl;
//    	 		std::cout << "DamageF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmegaFatigue() << std::endl;
//    	 		std::cout << "Kappa   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << std::endl;
//    	 		std::cout << "KappaF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappaFatigue() << std::endl;
//    	 		std::cout << "NonLEq  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() << std::endl;
//    	 		std::cout << "NonLEqF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrainFatigue() << std::endl;
//    	 		std::cout << "PrevEp  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrain() << std::endl;
//    	 		std::cout << "PrevEpF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrainFatigue() << std::endl;
//    	 		std::cout << "PrevS   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStress() << std::endl;
//    	 		std::cout << "PrevSF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStressFatigue() << std::endl;
//    	}
        if (mStructure->GetNumTimeDerivatives()>1)
        {
        	mStructure->NodeMergeDofValues(1,save_vel_j,save_vel_k);
        	mStructure->NodeMergeDofValues(2,save_acc_j,save_acc_k);
        }
        //capture the value of the internal force at the maximal displacement, prior to extrapolation
    	double timeDependentConstraintFactor(0);
        if (mTimeDependentConstraint!=-1)
        {
        	timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(SaveTime + 0.25/mHarmonicFactor(0,1));
        	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
        }
        mStructure->NodeMergeActiveDofValues(0,disp_Mean_j + disp_Ampl_j);																			//new place
        mStructure->ElementTotalUpdateTmpStaticData();																				//new place
    	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(cyclicPreJumpIntForceMax_j,cyclicPreJumpIntForceMax_k);	//new place
    	//the structure after the 3th cycle is already restored, except of displacements. Now restore the displacements and
    	//capture the respective internal force
        if (mTimeDependentConstraint!=-1)
        {
        	timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(SaveTime);
        	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
        }
        mStructure->NodeMergeActiveDofValues(0,save_disp_j);
        mStructure->ElementTotalUpdateTmpStaticData();
    	//the acceleration is close to zero after a full cycle, but not exactly equal to it, so the internal force too.
    	//This is because the real cycle, especially at the beginning of the fatigue history, is not perfectly harmonic (because statevs != const).
    	//We are interesting in the change of the internal force during extrapolation, thus we have to subtract
    	//the initial value of the internal force prior to extrapolation. In a static case this value is a zero and can be ignored.
    	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(cyclicPreJumpIntForce0_j,cyclicPreJumpIntForce0_k);	//new place

        std::cout<< " Starting Newmark " << std::endl;

    	Error = NuTo::NewmarkDirect::Solve(SaveTime + 1./mHarmonicFactor(0,1));											//integrate over the 4th cycle

    	//counter of iterations
    	int NumberOfIterations(0);
    	//initialize the number of cycles to be extrapolated
        int Njump(0);
        NuTo::FullVector<double,Eigen::Dynamic> NjumpVector;
        // set the number of cycles to be extrapolated in the cycle jump routine
        // ... NjumpVector[0] is the number of extrapolated cycles itself Njump
        // ... NjumpVector[1] is the weighting coefficient of the implicit term
        // ... NjumpVector[2] is the weighting coefficient of the explicit term
        // ... NjumpVector[3] and higher are the weighting coefficients of the terms for a higher-order extrapolation
        // the first three components are mandatory
        double NjumpDouble(0.);
//        NjumpDouble = mHarmonicExtrapolationTolerance*( cyclicPreJumpIntForceNext_j.array().abs() / (cyclicPreJumpIntForceNext_j - cyclicPreJumpIntForce_j).array().abs().max(mToleranceForce) ).minCoeff();
        NjumpDouble = mHarmonicExtrapolationTolerance*( cyclicPreJumpIntForceNext_j.array().abs() / (cyclicPreJumpIntForceNext_j - cyclicPreJumpIntForce_j).array().abs().maxCoeff() ).maxCoeff();
        Njump = std::max(NjumpMin,1 + boost::numeric_cast<int>(NjumpDouble));													//convert to integer, but not less than 3
        std::cout<<"INITIALIZATION NjumpDouble = " << NjumpDouble << ", Njump = " << Njump <<std::endl;

        // extrapolation methods
        bool Explicit(true);
        bool Theta(false); double theta(0.2);
        bool Ralston(false);
        bool Heun(false);
        bool Midpoint(false);

        //cycle jump loop
    	do {
    		//check whether Njump moves to outside of the time range rTimeDelta
    		NjumpDouble = (rTimeDelta - SaveTime)*mHarmonicFactor(0,1);													//number of cycles till the end of loading rTimeDelta
    		Njump = std::min(Njump,1 + boost::numeric_cast<int>(NjumpDouble));
    		std::cout<< "DO LOOP NjumpDuble = " << NjumpDouble << ", Njump = " << Njump << std::endl;

    		// set number of cycles
    		NjumpVector.Resize(3);
    		NjumpVector[0] = Njump; 				// number of cycles

    		//set weighting coefficients for explicit and implicit terms
    		if (NumberOfIterations >= 1 && Explicit == false) {
        		std::cout << "place impl" << std::endl;
    			//one explicit iteration has already been done, switch to implicit
        		if (Ralston) {
        			Njump *= 3./2.;
        			NjumpVector[0] = Njump;
        			NjumpVector[1] = 3./4.;			// implicit weighting coefficient
        			NjumpVector[2] = 1./4.;			// explicit weighting coefficient

        		} else if (Midpoint) {
        			Njump *= 2.;
        			NjumpVector[0] = Njump;
        			NjumpVector[1] = 1./2.;			// implicit weighting coefficient
        			NjumpVector[2] = 1./2.;			// explicit weighting coefficient

				} else if (Heun) {
        			NjumpVector[1] = 1./2.;			// implicit weighting coefficient
        			NjumpVector[2] = 1./2.;			// explicit weighting coefficient

				} else if (Theta) {
        			NjumpVector[1] = theta;			// implicit weighting coefficient
        			NjumpVector[2] = 1. - theta;	// explicit weighting coefficient
				}

	    		//save the values of the previous iteration, required for the termination criterion of the implicit iteration
        		//the values of two sequental iterations are compared
        		// this is only for implicit iterations
	  //  		cyclicPreJumpIntForce_j = cyclicAfterJumpIntForce_j;
	  //  		cyclicPreJumpIntForce_k = cyclicAfterJumpIntForce_k;

			} else {
				//the starting iteration is always explicit
	    		std::cout << "place expl" << std::endl;
	    		NjumpVector[1] = 0.;				// implicit weighting coefficient
	    		NjumpVector[2] = 1.;				// explicit weighting coefficient
			}
    		std::cout<< "Time " << SaveTime << ", Iteration " << NumberOfIterations << ", Njump " << NjumpVector[0] << std::endl;


    		//set NjumpVector to the structure
    		mStructure->SetNumExtrapolatedCycles(NjumpVector);

    		//extrapolate and set time (SaveTime is prior to the 4th cycle or prior to the ISC, whereas
    		//the actual time is the end of the 4th or the end of the ISC;
    		//the statevs are stored at the SaveTime, whereas the actual statevs correspond to the actual time)
    		mStructure->SetTime(SaveTime + (Njump + 1)/mHarmonicFactor(0,1));
    		mStructure->SetPrevTime(SaveTime + (Njump + 1)/mHarmonicFactor(0,1));

			//extrapolate state variables: statev += (statev - saved_statev) * Njump;
			mStructure->ElementFatigueExtrapolateStaticData();
//	    	std::cout << "Extrapolation READY =======================" <<std::endl;		// delete me
//	    	if (mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE){
//	    	 		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() << std::endl;
//	    	 		std::cout << "Damage  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << std::endl;
//	    	 		std::cout << "DamageF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmegaFatigue() << std::endl;
//	    	 		std::cout << "Kappa   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << std::endl;
//	    	 		std::cout << "KappaF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappaFatigue() << std::endl;
//	    	 		std::cout << "NonLEq  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() << std::endl;
//	    	 		std::cout << "NonLEqF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrainFatigue() << std::endl;
//	    	 		std::cout << "PrevEp  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrain() << std::endl;
//	    	 		std::cout << "PrevEpF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrainFatigue() << std::endl;
//	    	 		std::cout << "PrevS   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStress() << std::endl;
//	    	 		std::cout << "PrevSF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStressFatigue() << std::endl;
//	    	}

            //find equilibrium: update mean and amplitude displacements																					//new place
			this->CalculateFourierCoefficients(&disp_Mean_j,&disp_Mean_k,&disp_Ampl_j,&disp_Ampl_k,
					&cyclicPreJumpIntForce0_j,&cyclicPreJumpIntForce0_k,&cyclicPreJumpIntForceMax_j,&cyclicPreJumpIntForceMax_k);						//the structure is set to
																																						//the updated mean displacements

            //find equilibrium: update mean and amplitude displacements																					//new place
			this->CalculateFourierCoefficientsCoupledDofs(&disp_Mean_j,&disp_Mean_k,&disp_Ampl_j,&disp_Ampl_k);											//the structure is set to
																																						//the updated mean displacements
			//call UpdateStaticData in order to update mPrevStrain and mPrevSigma too (this has to
			//be done after the equilibrium has been found)
	        mStructure->ElementTotalUpdateStaticData();
//	    	std::cout << "After update static data =======================" <<std::endl;		// delete me
//	    	if (mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() == NuTo::Constitutive::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE){
//	    	 		std::cout << "Constitutive law = " << mStructure->ElementGetElementPtr(189)->GetConstitutiveLaw(0)->GetType() << std::endl;
//	    	 		std::cout << "Damage  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << std::endl;
//	    	 		std::cout << "DamageF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmegaFatigue() << std::endl;
//	    	 		std::cout << "Kappa   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << std::endl;
//	    	 		std::cout << "KappaF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappaFatigue() << std::endl;
//	    	 		std::cout << "NonLEq  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() << std::endl;
//	    	 		std::cout << "NonLEqF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrainFatigue() << std::endl;
//	    	 		std::cout << "PrevEp  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrain() << std::endl;
//	    	 		std::cout << "PrevEpF = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStrainFatigue() << std::endl;
//	    	 		std::cout << "PrevS   = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStress() << std::endl;
//	    	 		std::cout << "PrevSF  = " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevStressFatigue() << std::endl;
//	    	}

	        //evaluate cyclic increase of the out-of-balance intForce => integrate the subsequent single cycle (ISC)
	        //no postprocessing, therefore=>false
	        this->IntegrateSingleCycle(&disp_Mean_j,&disp_Mean_k,&disp_Ampl_j,&disp_Ampl_k,false);						//the statevs evolve through the cycle
	        mStructure->BuildGlobalGradientInternalPotentialSubVectors(cyclicAfterJumpIntForce_j,cyclicAfterJumpIntForce_k,false);

	        //compare the cyclic change of intForce prior and after the extrapolation (jump)
//	        JumpCriterionDouble = ( (cyclicAfterJumpIntForce_j - cyclicPreJumpIntForce_j).array().abs() /
//	        		cyclicAfterJumpIntForce_j.array().abs().max(mToleranceForce) ).maxCoeff();
	        JumpCriterionDouble = ( (cyclicAfterJumpIntForce_j - cyclicPreJumpIntForce_j).array().abs() /
	        		cyclicAfterJumpIntForce_j.array().abs().maxCoeff() ).maxCoeff();
	        JumpCriterion = ( JumpCriterionDouble <= mHarmonicExtrapolationTolerance );

	        std::cout << "Maximal component of cyclic change of Fint = "<<((cyclicAfterJumpIntForce_j - cyclicPreJumpIntForce_j).array().abs()).maxCoeff() << std::endl;

//	        for (int i = 0; i <= sizeof(cyclicAfterJumpIntForce_j); ++i){
//	        	std::cout << "i = " << i << " , cyclicPreJump = " << cyclicPreJumpIntForce_j[i] <<" , cycliAfterJump = " << cyclicAfterJumpIntForce_j[i] << std::endl;
//	        }
	        std::cout << "JumpCriterionDouble = " << JumpCriterionDouble << std::endl;

	        if ((JumpCriterion && ( Explicit || NumberOfIterations >= 1)) || Njump == NjumpMin || ((Ralston || Midpoint) && NumberOfIterations == 1)) {
//	        if ((JumpCriterion && ( Explicit || NumberOfIterations >= 1)) || Njump == NjumpMin ) {

	        	//jump is successful, continue with the next one

	        	//restart the counter of iterations for the next jump
	        	NumberOfIterations = 0;

//	        	//update velocities and accelerations due to slow evolution of the mean displacement and amplitude
//	        	//Explanation:
//	        	//Let u(t) = (u0 + Du0 t/ N T) + (u1 + Du1 t / N T) sin wt
//	        	//w = 2pi/T = 2pi*frequency, N = Njump, Du0, Du1 = changes of the mean displacement and amplitude after the jump
//	        	//The jump change of displacement due to the small changes Du0 and Du1 is equal to:
//	        	//Du = Du0 t / N T + (Du1 t / N T) sin wt
//	        	//The jump change of velocity Ddu and acceleration Dddu due to the small changes Du0 and Du1 are equal to:
//	        	//Ddu  = Du0 / N T + (Du1 / N T) sin wt  + (Du1 t/N T) d sin wt
//	        	//Dddu = (2 Du1 / N T) d sin wt + (Du1 t /N T) dd sin wt
//	        	//at  t = N T  (after the jump) we get
//	        	//Ddu  = Du0 / N T + Du1 w = frequency (Du0 / N + 2pi Du1)
//	        	//Dddu = 2 Du1 w / N T     = 4pi Du1 frequency^2/ N
	            if (mStructure->GetNumTimeDerivatives()>1)
	            {
	            	preJump_vel_j += ( (disp_Mean_j - preJump_disp_Mean_j)/Njump +											//extrapolate velocity
	            			2*pi*(disp_Ampl_j - preJump_disp_Ampl_j) )*mHarmonicFactor(0,1);
	            	preJump_vel_k += ( (disp_Mean_k - preJump_disp_Mean_k)/Njump +
	            			2*pi*(disp_Ampl_k - preJump_disp_Ampl_k) )*mHarmonicFactor(0,1);

	            	preJump_acc_j += 4*pi*(disp_Ampl_j - preJump_disp_Ampl_j)*												//extrapolate acceleration
	            			mHarmonicFactor(0,1)*mHarmonicFactor(0,1)/Njump;
	            	preJump_acc_k += 4*pi*(disp_Ampl_k - preJump_disp_Ampl_k)*
	            			mHarmonicFactor(0,1)*mHarmonicFactor(0,1)/Njump;
	            	mStructure->NodeMergeDofValues(1,preJump_vel_j,preJump_vel_k);
	            	mStructure->NodeMergeDofValues(2,preJump_acc_j,preJump_acc_k);
	            }
	            //save new static data and actual time for the next jump
	        	SaveTime = mStructure->GetTime();
	        	mStructure->ElementFatigueSaveStaticData();

	        	//update Njump
	        	NjumpDouble = Njump*
	        			std::sqrt(mHarmonicExtrapolationTolerance/std::max(JumpCriterionDouble,0.25*mHarmonicExtrapolationTolerance));
	    		if (Ralston) {
	    			NjumpDouble *= 2./3.;
	    		//	if (NjumpDouble < 3) NjumpDouble = 1.;			// to omit the loop, Njump = 3 -> Njump = 6 -> Njump = 3 ...
	    		} else if (Midpoint) {
	    			NjumpDouble *= 1./2.;
	    		//	if (NjumpDouble < 3) NjumpDouble = 1.;
	    		}
	    		Njump = std::max(NjumpMin,1 + boost::numeric_cast<int>(NjumpDouble));
	    		std::cout<< "NEW STEP, Njump = " << Njump << std::endl;

	    		//save deltaPreJumpIntForce and maximal internal force for the next jump
	    		cyclicPreJumpIntForce_j = cyclicAfterJumpIntForce_j;
	    		cyclicPreJumpIntForce_k = cyclicAfterJumpIntForce_k;

	            if (mTimeDependentConstraint!=-1)
	            {
	            	timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(SaveTime + 0.25/mHarmonicFactor(0,1));		//new place
	            	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);									//new place
	            }																															//new place
	            mStructure->NodeMergeActiveDofValues(0,disp_Mean_j + disp_Ampl_j);															//new place
	            mStructure->ElementTotalUpdateTmpStaticData();																				//new place
	        	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(cyclicPreJumpIntForceMax_j,cyclicPreJumpIntForceMax_k);	//new place

	            if (mTimeDependentConstraint!=-1)																							//new place
	            {																															//new place
	            	timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(SaveTime);									//new place
	            	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);									//new place
	            }																															//new place
	            mStructure->NodeMergeActiveDofValues(0,disp_Mean_j);																		//new place
	        	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(cyclicPreJumpIntForce0_j,cyclicPreJumpIntForce0_k);		//new place

	    		//integrate single cycle for postprocessing=true and for the next extrapolation of statevs
		        this->IntegrateSingleCycle(&disp_Mean_j,&disp_Mean_k,&disp_Ampl_j,&disp_Ampl_k,true);						//the statevs, vel and acc evolve through the postprocessing

		        //save node data for situation, if the jump has to be repeated
		        if (mStructure->GetNumTimeDerivatives()>1)
			    {
			        mStructure->NodeExtractDofValues(1,preJump_vel_j,preJump_vel_k);
			        mStructure->NodeExtractDofValues(2,preJump_acc_j,preJump_acc_k);
			    }
		        preJump_disp_Mean_j = disp_Mean_j; preJump_disp_Mean_k = disp_Mean_k;
		        preJump_disp_Ampl_j = disp_Ampl_j; preJump_disp_Ampl_k = disp_Ampl_k;

		        //go to the loop-beginning to start the next jump
			} else {
				//jump is too long, repeat again with a smaller Njump

				//update the value of the counter of iterations
				NumberOfIterations += 1;

				//restore the pre-jump time and statevs
				mStructure->SetPrevTime(SaveTime);
				mStructure->SetTime(SaveTime);
				mStructure->ElementFatigueRestoreStaticData();

				//restore displacements (vel and acc has not been changed and not need to be restored)
				mStructure->NodeMergeActiveDofValues(0,preJump_disp_Mean_j);
		    	mStructure->ElementTotalUpdateTmpStaticData();

		    	disp_Mean_j = preJump_disp_Mean_j; disp_Mean_k = preJump_disp_Mean_k;
		    	disp_Ampl_j = preJump_disp_Ampl_j; disp_Ampl_k = preJump_disp_Ampl_k;

		    	//decrease Njump
		    	if (Explicit || (( Heun || Theta) && NumberOfIterations >= 2)) {
		    		NjumpDouble = 0.75*Njump*std::sqrt(mHarmonicExtrapolationTolerance/JumpCriterionDouble);
		    		Njump = std::max(NjumpMin,1 + boost::numeric_cast<int>(NjumpDouble));
		    		NumberOfIterations = 0;
		    	}

	    		std::cout<< "REPEAT STEP, Njump = " << Njump << std::endl;

	    		//integrate single cycle for the next extrapolation of statevs, postprocessing=false
		        this->IntegrateSingleCycle(&disp_Mean_j,&disp_Mean_k,&disp_Ampl_j,&disp_Ampl_k,false);						//the statevs evolve through the cycle

		        //go to the loop-beginning to repeat the jump
			}
//			//test of statev update
//	        // straight forward integration after extrapolation
//	        mStructure->NodeMergeActiveDofValues(0,disp_Mean_j);
//	        mStructure->ElementTotalUpdateTmpStaticData();
//	    	std::cout << "mStress before USD" << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetPrevStress() << std::endl;
//	    	std::cout << "mOmegaCompr before USD" << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl;
//			// the state variables are already extrapolated, except of mPrevStrain and mPrevSigma
//			// call UpdateStaticData in order to update mPrevStrain and mPrevSigma too (this has to be done after the equilibrium has been found)
//	        mStructure->ElementTotalUpdateStaticData();
//		    std::cout << "mStress after USD " << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetPrevStress() << std::endl;
//		    std::cout << "mOmegaCompr after USD" << mStructure->ElementGetElementPtr(0)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl;
    	} while (SaveTime < rTimeDelta);
    }
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

    /*
    // ConstraintGetConstraintMatrixAfterGaussElimination(const NuTo::SparseMatrixCSRGeneral<double>& ... ) replaced by
    // const NuTo::SparseMatrixCSRGeneral<double> ConstraintGetConstraintMatrixAfterGaussElimination() const;
    // ---> No need for CmatTmp anymore! --- vhirtham

    NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
    mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
    */
    NuTo::SparseMatrixCSRVector2General<double> Cmat(mStructure->ConstraintGetConstraintMatrixAfterGaussElimination());
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
    mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);

	if (rFourierMode > 0 && mStructure->GetNumTimeDerivatives()>1)
		{
		double frequency(mHarmonicFactor(0,1));
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
//!@the structure is set to the updated mean displacements rDisp_Mean_j, and the RHS at the beginning of the cycle
void NuTo::JumpDirect::CalculateFourierCoefficients(NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Mean_j, NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Mean_k,
		NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Ampl_j, NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Ampl_k,
		const NuTo::FullVector<double,Eigen::Dynamic>* rIntForce_Mean_j, const NuTo::FullVector<double,Eigen::Dynamic>* rIntForce_Mean_k,
		const NuTo::FullVector<double,Eigen::Dynamic>* rIntForce_Max_j,  const NuTo::FullVector<double,Eigen::Dynamic>* rIntForce_Max_k)
{
    NuTo::FullVector<double,Eigen::Dynamic> extForce_j, extForce_k, extForceMax_j, extForceMax_k;
    NuTo::FullVector<double,Eigen::Dynamic> intForce_j(mStructure->GetNumActiveDofs()), intForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs()),
    		intForceMax_j(mStructure->GetNumActiveDofs()), intForceMax_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
    NuTo::FullVector<double,Eigen::Dynamic> residual_j, residual_k, residual_mod;
	NuTo::FullVector<double,Eigen::Dynamic> delta_disp_j, delta_disp_k, zero_displ_j(mStructure->GetNumActiveDofs());
    FullVector<double,Eigen::Dynamic> bRHS, bRHSmax;

    //calculate constraint matrix
    /*
    // ConstraintGetConstraintMatrixAfterGaussElimination(const NuTo::SparseMatrixCSRGeneral<double>& ... ) replaced by
    // const NuTo::SparseMatrixCSRGeneral<double> ConstraintGetConstraintMatrixAfterGaussElimination() const;
    // ---> No need for CmatTmp anymore! --- vhirtham

    NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
    mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
    */
    NuTo::SparseMatrixCSRVector2General<double> Cmat(mStructure->ConstraintGetConstraintMatrixAfterGaussElimination());
    SparseMatrixCSRVector2General<double> CmatT (Cmat.Transpose());

    // allocate space for stiffness matrix
    SparseMatrixCSRVector2General<double> stiffMatrix_jj(mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());


	if (mTimeDependentConstraint != -1 && mTimeDependentLoadCase != -1) {
		throw MechanicsException("[NuTo::JumpDirect::CalculateFourierCoefficients] the harmonic excitation is currently implemented for either time dependent constraint or load case.");
	}

	double BeginHarmonicExcitation;
	if (mTimeDependentConstraint != -1) {
		BeginHarmonicExcitation = mTimeDependentConstraintFactor(mTimeDependentConstraintFactor.GetNumRows()-1,0);
	} else if (mTimeDependentLoadCase != -1) {
		BeginHarmonicExcitation = mTimeDependentLoadFactor(mTimeDependentLoadFactor.GetNumRows()-1,0);
	} else {
		throw MechanicsException("[NuTo::JumpDirect::CalculateFourierCoefficients] the harmonic excitation can be only applied to the time dependent constraint or load case.");
	}

	//RemCom
	NuTo::FullVector<double,Eigen::Dynamic> rIntForceX_Max_j, rIntForceX_Max_k, rIntForceX_Mean_j, rIntForceX_Mean_k,
		rDispX_Mean_j, rDispX_Mean_k, rDispX_Ampl_j, rDispX_Ampl_k;
	rIntForceX_Max_j = *rIntForce_Max_j; rIntForceX_Max_k = *rIntForce_Max_k;
	rIntForceX_Mean_j = *rIntForce_Mean_j; rIntForceX_Mean_k = *rIntForce_Mean_k;
	rDispX_Mean_j = *rDisp_Mean_j; rDispX_Mean_k = *rDisp_Mean_k;
	rDispX_Ampl_j = *rDisp_Ampl_j; rDispX_Ampl_k = *rDisp_Ampl_k;

	// find the mean value (rFourier = 0)

	// set BC for the mean displacement, which is the displacement at the beginning of the harmonic excitation
	double timeDependentConstraintFactor(0);
    if (mTimeDependentConstraint!=-1)
    {
    	timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(BeginHarmonicExcitation);
    	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
    }

    bRHS = mStructure->ConstraintGetRHSAfterGaussElimination();
	// set mean displacement field
	mStructure->NodeMergeActiveDofValues(0,*rDisp_Mean_j);
	mStructure->ElementTotalUpdateTmpStaticData();
	// calculate internal forces and external forces
	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(intForce_j,intForce_k);
	CalculateExternalLoad(*mStructure, BeginHarmonicExcitation, extForce_j, extForce_k);

	residual_j = intForce_j - 0*extForce_j;
	residual_k = intForce_k - 0*extForce_k;

	residual_j -= *rIntForce_Mean_j;				//new place
	residual_k -= *rIntForce_Mean_k;				//new place


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
            mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);
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
	std::cout << "CFC: prior mean solver ===================================" << std::endl;

	// solve for the mean displacement
	NuTo::SparseMatrixCSRGeneral<double> hessianModSolver(stiffMatrix_jj);

	hessianModSolver.SetOneBasedIndexing();
    // allocate solver
    NuTo::SparseDirectSolverMUMPS mySolver;
	mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
	delta_disp_j*=-1;

	std::cout << "CFC: after mean solver ===================================" << std::endl;

	// calculate mean displacement
	*rDisp_Mean_j += delta_disp_j;
	*rDisp_Mean_k = bRHS - Cmat*(*rDisp_Mean_j);

	////////**************
	// the previous calculation of disp_j is based on the update scheme
	// the following calculation is based on direct calculation, that is:
	// Kmod * rDisp_mean_j = -( (K12 - CmatT*K22)*bRHS - Fext ), Fint = 0
	bool DirectCalculation(true);
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

        // contribution of Fint due to plastic strain, calculate Fint(u=0)	// new place
        if (mTimeDependentConstraint!=-1)
        {
        	mStructure->ConstraintSetRHS(mTimeDependentConstraint,0.);		// set zero rhs displacement
        }
      	zero_displ_j.Zero(mStructure->GetNumActiveDofs());
       	mStructure->NodeMergeActiveDofValues(0,zero_displ_j);				// set zero active displacement
       	mStructure->ElementTotalUpdateTmpStaticData();
       	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(intForce_j,intForce_k);	// calculate internal forces at zero displ.

       	// add plastic contribution to residual								// new place
        residual_j += intForce_j;											// new place
        residual_k += intForce_k;											// new place

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
    }
    	// update internal force
    	mStructure->NodeMergeActiveDofValues(0,*rDisp_Mean_j);
    	mStructure->ElementTotalUpdateTmpStaticData();
    	// calculate internal force
    	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(intForce_j,intForce_k);


	////////**************

	// find the amplitude (rFourier = 1)

	// we need to calculate the change of intForce(maximal displacement) - intForce(mean displacement) due to extrapolation of statevs
	// calculate time at which the harmonic displacement is maximal, this is a quarter of the period
	double MaximumHarmonicExcitation(BeginHarmonicExcitation + 0.25/mHarmonicFactor(0,1));
	// set BC for the maximal displacement
	double timeDependentConstraintFactorMax(0);
    if (mTimeDependentConstraint != -1)
    {
    	timeDependentConstraintFactorMax = this->CalculateTimeDependentConstraintFactor(MaximumHarmonicExcitation);
    	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactorMax);
    }

    //set load amplitude for a load controlled case
    if (mTimeDependentLoadCase != -1) {
    	CalculateExternalLoad(*mStructure, MaximumHarmonicExcitation, extForceMax_j, extForceMax_k);
    	extForce_j = extForceMax_j - extForce_j;
    	extForce_k = extForceMax_k - extForce_k;
	} else {
		extForce_j.Zero(mStructure->GetNumActiveDofs()); extForce_k.Zero(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
	}

    bRHSmax= mStructure->ConstraintGetRHSAfterGaussElimination();
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

	// calculate the change of residuals and stiffness
	residual_j.Zero(mStructure->GetNumActiveDofs()); residual_k.Zero(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());

	residual_j = intForce_j;
	residual_k = intForce_k;

	residual_j -= (*rIntForce_Max_j - *rIntForce_Mean_j);				//new place
	residual_k -= (*rIntForce_Max_k - *rIntForce_Mean_k);				//new place

	// the external forces do not experience any change via statevs and displ. evolutions, therefore excluded
//	if (mTimeDependentLoadCase != -1)									//new place
//	{
//		residual_j -= extForce_j;
//		residual_k -= extForce_k;
//	}

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
        mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);
		std::cout << "use full mass matrix " << std::endl;
	}

	double frequency(mHarmonicFactor(0,1));
	const double pi = boost::math::constants::pi<double>();
	double factor(2*pi*frequency);

	// the impact of mass seems to be redundant because it is already introduced in the stiffness matrix
//	if (mStructure->GetNumTimeDerivatives()>1) {
//		if (mUseLumpedMass)
//		{
//			residual_j -= lumped_massMatrix_j.asDiagonal()*(*rDisp_Ampl_j)*factor*factor;
//			residual_k -= lumped_massMatrix_k.asDiagonal()*(*rDisp_Ampl_k)*factor*factor;
//		}
//		else
//		{
//			//add mass terms
//			residual_j -= (massMatrix_jj*(*rDisp_Ampl_j)+massMatrix_jk*(*rDisp_Ampl_k))*factor*factor;
//			residual_k -= (massMatrix_kj*(*rDisp_Ampl_j)+massMatrix_kk*(*rDisp_Ampl_k))*factor*factor;
//		}
//	}

	if (Cmat.GetNumEntries()>0)
	{
		residual_mod = residual_j - CmatT*residual_k;
	} else {
		residual_mod = residual_j;
	}

	stiffMatrix_jj.SetZeroEntries();

	this->CalculateGlobalModifiedStiffness(&stiffMatrix_jj,1);
	std::cout << "CFC: prior ampl. solver ===================================" << std::endl;
	// solve for the ampl displacement
	hessianModSolver = stiffMatrix_jj;

	hessianModSolver.SetOneBasedIndexing();
	delta_disp_j.Zero(mStructure->GetNumActiveDofs());
	mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
	delta_disp_j*=-1;
	std::cout << "CFC: after ampl. solver ===================================" << std::endl;


	// calculate ampl displacement
	*rDisp_Ampl_j += delta_disp_j;
	*rDisp_Ampl_k = bRHS - Cmat*(*rDisp_Ampl_j);

	////////**************
	// the previous calculation of rDisp_Ampl_j is based on the update scheme
	// the following calculation is based on direct calculation, that is:
	// Kmod * rDisp_mean_j = -( (K12 - CmatT*K22)*delta_bRHS ), Fint = 0
	DirectCalculation = true;  //displ controlled set true!, load controlled set false (or true, which overestimates than the displ. amplitude)
    if (DirectCalculation)
    {
    	// calculate stiffness for the 1 mode
        K_jj.AddScal(massMatrix_jj,-factor*factor);
        K_jk.AddScal(massMatrix_jk,-factor*factor);
        K_kj.AddScal(massMatrix_kj,-factor*factor);
        K_kk.AddScal(massMatrix_kk,-factor*factor);

    	// contribution of bRHS to residual
        residual_j = K_jk * bRHS;
        residual_k = K_kk * bRHS;

    	// contribution of the amplitude of external forces to residual
    	if (mTimeDependentLoadCase != -1)
    	{
    		residual_j -= extForce_j;
    		residual_k -= extForce_k;
    	}

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

    //RemCom
//
//	double MaxDeltaFint_j(0.), MaxDeltaFint_k(0.);
//
//	std::cout << "mStructure->GetNumDofs() = " << mStructure->GetNumDofs() << " mStructure->GetNumActiveDofs() = " << mStructure->GetNumActiveDofs() << std::endl;
//	std::cout << "sizeof(intForce_j) = " << sizeof(intForce_j) << " sizeof(intForce_k) = " << sizeof(intForce_k) << " sizeof(*rintForce_Max_j)" <<
//			sizeof(*rIntForce_Max_j) << " sizeof(*rintForce_Max_k)" << sizeof(*rIntForce_Max_k) << std::endl;
//	std::cout << "*rIntForce_Max_j -*rIntForce_Mean_j" << std::endl;
//	std::cout << *rIntForce_Max_j - *rIntForce_Mean_j << std::endl;
//	std::cout << " rDispX_Mean_j | rDispX_Ampl_j | delta_displ_j | rFintPrior_Max_j  |  FintAfterJump_Max_j  |  Delta_rFintPrior_j  |  Delta_FintAfterJump_j" << std::endl;
//	for (int i = 0; i < mStructure->GetNumActiveDofs(); ++i) {
//		std::cout << "i = " << i << " | " << rDispX_Mean_j[i] << " | " << rDispX_Ampl_j[i] << " | " << delta_disp_j[i]
//			<< " | " <<	rIntForceX_Max_j[i] << " | " << intForceMax_j[i]
//			<< " | " << rIntForceX_Max_j[i] - rIntForceX_Mean_j[i] << " | " << intForce_j[i] << std::endl;
//		if (abs(rIntForceX_Max_j[i] - rIntForceX_Mean_j[i]-intForce_j[i]) > MaxDeltaFint_j) {
//			MaxDeltaFint_j = abs(rIntForceX_Max_j[i] - rIntForceX_Mean_j[i]-intForce_j[i]);
//		}
//
//	}
//
//	std::cout << " rDispX_Mean_k | rDispX_Ampl_k |rFintPrior_Max_k  |  FintAfterJump_Max_k  |  Delta_rFintPrior_k  |  Delta_FintAfterJump_k" << std::endl;
//	for (int i = 0; i < mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(); ++i) {
//		std::cout << "i = " << i << " | " << rDispX_Mean_k[i] << " | " << rDispX_Ampl_k[i]
//			<< " | " <<	rIntForceX_Max_k[i] << " | " << intForceMax_k[i]
//			<< " | " << rIntForceX_Max_k[i] - rIntForceX_Mean_k[i] << " | " << intForce_k[i] << std::endl;
//		if (abs(rIntForceX_Max_k[i] - rIntForceX_Mean_k[i] - intForce_k[i]) > MaxDeltaFint_k) {
//			MaxDeltaFint_k = abs(rIntForceX_Max_k[i] - rIntForceX_Mean_k[i] - intForce_k[i]);
//		}
//	}
//
//	std::cout << " MaxDeltaFint_j = " << MaxDeltaFint_j << ", MaxDeltaFint_k = " << MaxDeltaFint_k << std::endl;

	////////**************
    // set the new mean displacements and the respective displacement BC
    if (mTimeDependentConstraint!=-1)
    {
    	timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(BeginHarmonicExcitation);
    	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
    }
	// set mean displacement field
	mStructure->NodeMergeActiveDofValues(0,*rDisp_Mean_j);
	mStructure->ElementTotalUpdateTmpStaticData();
}

//!@brief updates DofTypes after the CalculateFourierCoefficients routine
// The CalculateFourierCoefficients routine determines the DISPLACEMENTS Dof only. For a coupled problem DISPLACEMENTS/DofType, the another DofType
// has to be updated too. This routine updates the following list of Dofs:
// NONLOCALEQSTRAIN
void NuTo::JumpDirect::CalculateFourierCoefficientsCoupledDofs(NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Mean_j, NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Mean_k,
		NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Ampl_j, NuTo::FullVector<double,Eigen::Dynamic>* rDisp_Ampl_k)
{
	std::map<NuTo::Node::eAttributes, bool> updatedDofs;
	bool updateCoupledDofs(false);

	// creating the map of those Dofs, which can be updated by this routine
	updatedDofs[NuTo::Node::NONLOCALEQSTRAIN] = false;
	// here add further Dofs

	// loop over those updatedDofs, which can be updated by this routine, in order to find out whether at least one Dof from updatedDofs is active in the structure
	for (auto it = updatedDofs.begin(); it != updatedDofs.end(); it++){
		if(mStructure->InterpolationTypeIsConstitutiveInput(it->first)){
			it->second = true;
			updateCoupledDofs = true;
			std::cout<< "[NuTo::JumpDirect::CalculateFourierCoefficientsCoupledDofs] the following Dof will be updated: " << NuTo::Node::AttributeToString(it->first) <<std::endl;
		}
	}

    // if any of Dofs has to be updated, do nothing and return
    if (updateCoupledDofs == false) {
		return;
	}

    if (mStructure->GetNumTimeDerivatives()>1)
    {
    	throw MechanicsException("[NuTo::JumpDirect::CalculateFourierCoefficientsCoupledDofs] the mass matrix is not implemented for the Coupled Dofs, only static integration is possible.");
    }

//    // print the Dofs to be updated
//	for (auto it = updatedDofs.begin(); it != updatedDofs.end(); it++){
//		if(it->second){
//			std::cout<< "[NuTo::JumpDirect::CalculateFourierCoefficientsCoupledDofs] the following Dof will be updated" << NuTo::Node::AttributeToString(it->first) <<std::endl;
//		}
//	}

	// starting update the Dofs

    NuTo::FullVector<double,Eigen::Dynamic> extForce_j, extForce_k, extForceMax_j, extForceMax_k;
    NuTo::FullVector<double,Eigen::Dynamic> intForce_j(mStructure->GetNumActiveDofs()), intForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs()),
    		intForceMax_j(mStructure->GetNumActiveDofs()), intForceMax_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
    NuTo::FullVector<double,Eigen::Dynamic> residual_j, residual_k, residual_mod;
	NuTo::FullVector<double,Eigen::Dynamic> delta_disp_j, delta_disp_k, zero_displ_j(mStructure->GetNumActiveDofs());
    FullVector<double,Eigen::Dynamic> bRHS, bRHSmax;

    //calculate constraint matrix
    NuTo::SparseMatrixCSRVector2General<double> Cmat(mStructure->ConstraintGetConstraintMatrixAfterGaussElimination());
    SparseMatrixCSRVector2General<double> CmatT (Cmat.Transpose());

    // allocate space for stiffness matrix
    SparseMatrixCSRVector2General<double> stiffMatrix_jj(mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());


	if (mTimeDependentConstraint != -1 && mTimeDependentLoadCase != -1) {
		throw MechanicsException("[NuTo::JumpDirect::CalculateFourierCoefficientsCoupledDofs] the harmonic excitation is currently implemented for either time dependent constraint or load case.");
	}

	double BeginHarmonicExcitation;
	if (mTimeDependentConstraint != -1) {
		BeginHarmonicExcitation = mTimeDependentConstraintFactor(mTimeDependentConstraintFactor.GetNumRows()-1,0);
	} else if (mTimeDependentLoadCase != -1) {
		BeginHarmonicExcitation = mTimeDependentLoadFactor(mTimeDependentLoadFactor.GetNumRows()-1,0);
	} else {
		throw MechanicsException("[NuTo::JumpDirect::CalculateFourierCoefficientsCoupledDofs] the harmonic excitation can be only applied to the time dependent constraint or load case.");
	}

	// update the mean value (rFourier = 0)

	// set BC for the mean displacement, which is the displacement at the beginning of the harmonic excitation
	double timeDependentConstraintFactor(0);
    if (mTimeDependentConstraint!=-1)
    {
    	timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(BeginHarmonicExcitation);
    	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
    }

    bRHS = mStructure->ConstraintGetRHSAfterGaussElimination();
	// set mean displacement field
	mStructure->NodeMergeActiveDofValues(0,*rDisp_Mean_j);
	mStructure->ElementTotalUpdateTmpStaticData();

	// calculate internal forces and external forces
	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(intForce_j,intForce_k);
	CalculateExternalLoad(*mStructure, BeginHarmonicExcitation, extForce_j, extForce_k);

	//	calculate residuals
	residual_j = intForce_j - extForce_j;
	residual_k = intForce_k - extForce_k;

	// an inertial term with the mass matrix should be added here for dynamic equilibrium

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

	// find the amplitude (rFourier = 1)

	// we need to calculate the change of intForce(maximal displacement) - intForce(mean displacement) due to extrapolation of statevs
	// calculate time at which the harmonic displacement is maximal, this is a quarter of the period
	double MaximumHarmonicExcitation(BeginHarmonicExcitation + 0.25/mHarmonicFactor(0,1));
	// set BC for the maximal displacement
	double timeDependentConstraintFactorMax(0);
    if (mTimeDependentConstraint != -1)
    {
    	timeDependentConstraintFactorMax = this->CalculateTimeDependentConstraintFactor(MaximumHarmonicExcitation);
    	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactorMax);
    }

    //set load amplitude for a load controlled case
    if (mTimeDependentLoadCase != -1) {
    	CalculateExternalLoad(*mStructure, MaximumHarmonicExcitation, extForceMax_j, extForceMax_k);
    	extForce_j = extForceMax_j - extForce_j;
    	extForce_k = extForceMax_k - extForce_k;
	} else {
		extForce_j.Zero(mStructure->GetNumActiveDofs()); extForce_k.Zero(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
	}

    bRHSmax= mStructure->ConstraintGetRHSAfterGaussElimination();
	// set maximal displacement
	mStructure->NodeMergeActiveDofValues(0,(*rDisp_Mean_j) + (*rDisp_Ampl_j));
	mStructure->ElementTotalUpdateTmpStaticData();

	// calculate intForce(maximal displacement)
	mStructure->BuildGlobalElasticGradientInternalPotentialSubVectors(intForceMax_j,intForceMax_k);
	intForce_j = intForceMax_j - intForce_j; // intForce is currently the intForce(meanDisplacement)
	intForce_k = intForceMax_k - intForce_k;

//	calculate bRHS
	bRHS = bRHSmax - bRHS;

	// calculate the residuals
	residual_j.Zero(mStructure->GetNumActiveDofs()); residual_k.Zero(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());

	residual_j = intForce_j - extForce_j;
	residual_k = intForce_k - extForce_k;

	// an inertial term with the mass matrix should be added here for dynamic equilibrium

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

    // set the new mean displacements and the respective displacement BC
    if (mTimeDependentConstraint!=-1)
    {
    	timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(BeginHarmonicExcitation);
    	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
    }
	// set mean displacement field
	mStructure->NodeMergeActiveDofValues(0,*rDisp_Mean_j);
	mStructure->ElementTotalUpdateTmpStaticData();

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
    /*
    // ConstraintGetConstraintMatrixAfterGaussElimination(const NuTo::SparseMatrixCSRGeneral<double>& ... ) replaced by
    // const NuTo::SparseMatrixCSRGeneral<double> ConstraintGetConstraintMatrixAfterGaussElimination() const;
    // ---> No need for CmatTmp anymore! --- vhirtham

    NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
    mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
    */
    NuTo::SparseMatrixCSRVector2General<double> Cmat(mStructure->ConstraintGetConstraintMatrixAfterGaussElimination());
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
            mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);
    		std::cout << "use full mass matrix " << std::endl;
    	}
    }

    //extract displacements, velocities and accelerations
    NuTo::FullVector<double,Eigen::Dynamic> disp_j,disp_k, vel_j, vel_k, acc_j, acc_k;
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

    bRHSprev = mStructure->ConstraintGetRHSAfterGaussElimination();
    mStructure->NodeMergeActiveDofValues(0,lastConverged_disp_j); //disp_k is internally calculated from the previous constraint matrix
    mStructure->ElementTotalUpdateTmpStaticData();

    //calculate internal force, update of history variables=false
    mStructure->BuildGlobalGradientInternalPotentialSubVectors(prevIntForce_j,prevIntForce_k,false);
    prevResidual_j = prevIntForce_j;
    prevResidual_k = prevIntForce_k;
//    std::cout << "     RES.: prevIntForce_j = " << prevIntForce_j.Norm() << std::endl;

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
//		    std::cout << "     RES.: prevIntForce_j + M*acc_j = " << prevResidual_j.Norm() << std::endl;
    	}
    }
    //add external force
    CalculateExternalLoad(*mStructure, curTime, prevExtForce_j, prevExtForce_k);
    prevResidual_j -= prevExtForce_j;
    prevResidual_k -= prevExtForce_k;
//    std::cout << "     RES.: prevExtForce_j = " << prevExtForce_j.Norm() << std::endl;


    if (Cmat.GetNumEntries()>0)
    {
        residual_mod=prevResidual_j - CmatT*prevResidual_k;
    }
    else
    {
        residual_mod=prevResidual_j;
    }

//    if (mStructure->GetNumTimeDerivatives()>1) {
//    	std::cout << "=== vel_j = " << lastConverged_vel_j[0] << std::endl;
//    }

//    if (rIncludePostProcess && residual_mod.Norm()>mToleranceForce){
//        std::cout << "residual in initial configuration " << residual_mod.Norm() << std::endl;
//        throw MechanicsException("[NuTo::JumpDirect::Solve] Configuration after the extrapolation jump is not in (dynamic) equilibrium.");
//    }

//    std::ofstream DamageFile; std::ofstream TotalInelasticEqStrainFile;
//    DamageFile.open("DamageJump.txt", std::ios::app);
//    TotalInelasticEqStrainFile.open("TotalInelasticEqStrainJump.txt", std::ios::app);

    // for 2d gradient model
    std::ofstream DamageFile;
    DamageFile.open("Damage.txt", std::ios::app);


    if (rIncludePostProcess){
    	mTime = curTime;
    	NuTo::NewmarkDirect::PostProcess(prevResidual_j, prevResidual_k);
//    	DamageFile << mTime << " " << mStructure->ElementGetElementPtr(14)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl; // for Brick8N was ElementGetElementPtr(9) and ElementGetElementPtr(141) dlya Brick8Nhole;
//    	TotalInelasticEqStrainFile << mTime << " " << mStructure->ElementGetElementPtr(14)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetKappaInelastic() << std::endl; // for Brick8N was ElementGetElementPtr(9) and ElementGetElementPtr(141) dlya Brick8Nhole;
        // for 2d gradient model
//    	DamageFile << mTime << " " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega()
//                       		<< " " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << " " << std::endl;
    	// 2D nonlocal test
//        DamageFile << mTime << " " << mStructure->ElementGetElementPtr(339)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega()
//        		<< " " << mStructure->ElementGetElementPtr(339)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << " "
//				<< mStructure->ElementGetElementPtr(339)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() <<
//				" " << mStructure->ElementGetElementPtr(406)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << " "
//				<< mStructure->ElementGetElementPtr(406)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << " "
//				<< mStructure->ElementGetElementPtr(406)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain()<< std::endl;
    }

    // iterate over the increments within the cycle
	const double pi = boost::math::constants::pi<double>();
    double frequency(mHarmonicFactor(0,1));
    double timeStep(1./(mHarmonicIncrementation*frequency));
    for (int incr = 1; incr <= mHarmonicIncrementation; ++incr) {
    	// set new time
    	curTime += timeStep;
    	mStructure->SetTime(curTime);

    	// set BC
        if (mTimeDependentConstraint!=-1)
        {
        	timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
        	mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
        }
    	bRHSend.Zero(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
        bRHSend = mStructure->ConstraintGetRHSAfterGaussElimination();

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
//	    std::cout << "     RES.: prevIntForce_j = " << prevIntForce_j.Norm() << std::endl;

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
                bRHShalf = mStructure->ConstraintGetRHSAfterGaussElimination();

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
//                acc_k = bRHSddot - (Cmat*acc_j);
//                vel_k = bRHSdot - (Cmat*vel_j);
                acc_k = (lastConverged_vel_k+lastConverged_acc_k*(timeStep*(0.5-mBeta)))*(-1./(timeStep*mBeta));
                acc_k += bRHSdot * (1./(timeStep*mBeta));
                vel_k = lastConverged_vel_k+lastConverged_acc_k*((1.-mGamma)*timeStep)+acc_k*(mGamma*timeStep);

               // direct calculation of vel and acc
//        		vel_j = 2*pi*frequency*(*rDisp_Ampl_j)*cos(2*pi*frequency*(curTime - startingCycleTime));
//        		acc_j = -pow(2*pi*frequency,2)*(*rDisp_Ampl_j)*sin(2*pi*frequency*(curTime - startingCycleTime));
//
//        		vel_k = 2*pi*frequency*(*rDisp_Ampl_k)*cos(2*pi*frequency*(curTime - startingCycleTime));
//        		acc_k = -pow(2*pi*frequency,2)*(*rDisp_Ampl_k)*sin(2*pi*frequency*(curTime - startingCycleTime));

//    			std::cout << " i | AfterJumpAccMax_k " << std::endl;
//    			for (int i = 0; i < mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(); ++i) {
//    				std::cout << "i = " << i << " | " << acc_k[i] << std::endl;
//    			}

                // calculate the dynamic contribution to the residuals
                if (mMuDampingMass>0)
                {
                	//add damping terms
                    if (mUseLumpedMass)
                    {
                    	prevResidual_j += lumped_massMatrix_j.asDiagonal()*vel_j * mMuDampingMass;
                    	prevResidual_k += lumped_massMatrix_k.asDiagonal()*vel_k * mMuDampingMass;
                    }
                    else
                    {
                		prevResidual_j += (massMatrix_jj*vel_j+massMatrix_jk*vel_k) *mMuDampingMass;
                		prevResidual_k += (massMatrix_kj*vel_j+massMatrix_kk*vel_k) *mMuDampingMass;
                    }
                }
                if (mUseLumpedMass)
                {
                	prevResidual_j += lumped_massMatrix_j.asDiagonal()*acc_j;
                	prevResidual_k += lumped_massMatrix_k.asDiagonal()*acc_k;
                }
                else
                {
                //add mass terms
                	prevResidual_j += (massMatrix_jj*acc_j+massMatrix_jk*acc_k);
                	prevResidual_k += (massMatrix_kj*acc_j+massMatrix_kk*acc_k);
//                  std::cout << "     RES.: prevIntForce_j + M*acc_j = " << prevResidual_j.Norm() << std::endl;
                }

            }

        	std::cout << " Cyclic increment " << incr << std::endl;
        	NuTo::NewmarkDirect::PostProcess(prevResidual_j, prevResidual_k);
//        	DamageFile << mTime << " " << mStructure->ElementGetElementPtr(14)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetOmegaCompr() << std::endl; // for Brick8N was ElementGetElementPtr(9) and ElementGetElementPtr(141) dlya Brick8Nhole;
//        	TotalInelasticEqStrainFile << mTime << " " << mStructure->ElementGetElementPtr(14)->GetStaticData(0)->AsDamageViscoPlasticity3D()->GetKappaInelastic() << std::endl; // for Brick8N was ElementGetElementPtr(9) and ElementGetElementPtr(141) dlya Brick8Nhole;
            // for 2d gradient model
//        	DamageFile << mTime << " " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega()
//                           		<< " " << mStructure->ElementGetElementPtr(189)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << " " << std::endl;
        	// 2D nonlocal test
//            DamageFile << mTime << " " << mStructure->ElementGetElementPtr(339)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega()
//            		<< " " << mStructure->ElementGetElementPtr(339)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << " "
//    				<< mStructure->ElementGetElementPtr(339)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain() <<
//    				" " << mStructure->ElementGetElementPtr(406)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetOmega() << " "
//    				<< mStructure->ElementGetElementPtr(406)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetKappa() << " "
//    				<< mStructure->ElementGetElementPtr(406)->GetStaticData(0)->AsGradientDamage2DFatigue()->GetPrevNonlocalEqStrain()<< std::endl;

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
//    DamageFile.close(); TotalInelasticEqStrainFile.close();
    // 2D nonlocal test
//    DamageFile.close();
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
