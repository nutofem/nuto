// $Id: NystroemBase.cpp 575 2011-09-20 18:05:35Z unger3 $

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

#include "nuto/math/SparseDirectSolverMUMPS.h"

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/NystroemBase.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::NystroemBase::NystroemBase (StructureBase* rStructure)  : TimeIntegrationBase (rStructure)
{
    mTimeStep = 0.;
    mUseDiagonalMassMatrix = true;
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::NystroemBase::Info()const
{
	TimeIntegrationBase::Info();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NystroemBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NystroemBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NystroemBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NystroemBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NystroemBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NystroemBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NystroemBase::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of NystroemBase" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase)
           & BOOST_SERIALIZATION_NVP(mTimeStep)
           & BOOST_SERIALIZATION_NVP(mUseDiagonalMassMatrix);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of NystroemBase" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
NuTo::Error::eError NuTo::NystroemBase::Solve(double rTimeDelta)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    try
    {
        if (mTimeStep==0.)
        {
        	if (this->HasCriticalTimeStep())
        	{
        		mTimeStep = this->CalculateCriticalTimeStep();
        	}
        	else
        	{
                throw MechanicsException("[NuTo::NystroemBase::Solve] time step not set for unconditionally stable algorithm.");
        	}
        }
        //allocate solver (usually not needed)
        NuTo::SparseDirectSolverMUMPS mySolver;

        std::cout << "modify computation of critical time step to include the dependence on the time integration scheme." <<std::endl;
        //calculate instead the smallest eigenfrequency, depending on the time integration this gives the critical time step

        std::cout << "time step " << mTimeStep << std::endl;
        std::cout << "number of time steps " << rTimeDelta/mTimeStep << std::endl;

        //renumber dofs and build constraint matrix
        mStructure->NodeBuildGlobalDofs();

        //calculate constraint matrix
        NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
        mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
        NuTo::SparseMatrixCSRVector2General<double> Cmat(CmatTmp);
        SparseMatrixCSRVector2General<double> CmatT (Cmat.Transpose());
        //FullVector<double,Eigen::Dynamic> bRHS;
        if (CmatT.GetNumEntries() > 0)
        {
            throw MechanicsException("[NuTo::NystroemBase::Solve] not implemented for constrained systems including multiple dofs.");
        }

        //calculate individual inverse mass matrix, use only lumped mass matrices - stored as fullvectors and then use asDiagonal()
        NuTo::FullVector<double,Eigen::Dynamic> invMassMatrix_j;
        NuTo::SparseMatrixCSRGeneral<double> fullMassMatrix_jj;

        //extract displacements, velocities and accelerations
        NuTo::FullVector<double,Eigen::Dynamic> disp_j, vel_j, tmp_k,
                                                disp_j_tmp,              //temporary argument for the evalution of f(y)
                                                disp_j_new, vel_j_new;   //new d and v at end of time step
        std::vector<NuTo::FullVector<double,Eigen::Dynamic> > acc_j_tmp; //intermediate values of the time derivatives (Z or l)

        NuTo::FullVector<double,Eigen::Dynamic> extForce_j, extForce_k;
        NuTo::FullVector<double,Eigen::Dynamic> outOfBalance_j(mStructure->GetNumActiveDofs()), outOfBalance_k;
        NuTo::FullVector<double,Eigen::Dynamic> intForce_j(mStructure->GetNumActiveDofs()),
        		                                intForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());

        //store last converged displacements, velocities and accelerations
        mStructure->NodeExtractDofValues(0,disp_j,tmp_k);
        mStructure->NodeExtractDofValues(1,vel_j,tmp_k);

        if (mUseDiagonalMassMatrix)
        {
            //calculate lumped mass matrix
        	invMassMatrix_j.resize(mStructure->GetNumActiveDofs());
            //intForce_j is just used as a tmp variable
			mStructure->BuildGlobalLumpedHession2(intForce_j,tmp_k);

			//check the sum of all entries
			std::cout << "the total mass is " << intForce_j.sum()/3. +  tmp_k.sum()/3. << std::endl;

			//invert the mass matrix
			invMassMatrix_j = intForce_j.cwiseInverse();
        }
        else
        {
        	NuTo::SparseMatrixCSRVector2General<double> full2MassMatrix_jj,full2MassMatrix_jk,full2MassMatrix_kj,full2MassMatrix_kk;
        	full2MassMatrix_jj.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
        	full2MassMatrix_jk.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
        	full2MassMatrix_kj.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
        	full2MassMatrix_kk.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
			mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::MASS,full2MassMatrix_jj,full2MassMatrix_jk,full2MassMatrix_kj,full2MassMatrix_kk);
			fullMassMatrix_jj=full2MassMatrix_jj;
			fullMassMatrix_jj.SetOneBasedIndexing();
	    }

        double curTime  = 0;
		CalculateExternalLoad(*mStructure, curTime, extForce_j, extForce_k);
		acc_j_tmp.resize(this->GetNumStages());
		std::vector<double> stageDerivativeFactor(this->GetNumStages()-1);
        while (curTime < rTimeDelta)
        {
         	//calculate for delta_t = 0
            std::cout << "curTime " << curTime <<   " (" << curTime/rTimeDelta << ") max Disp = "  <<  disp_j.maxCoeff() << std::endl;
        	disp_j_new = disp_j + vel_j*mTimeStep;
        	vel_j_new = vel_j;

			double prevTime(mTime);
			double prevCurTime(curTime);
        	for (int countStage=0; countStage<this->GetNumStages(); countStage++)
            {
        		//std::cout << "\n stage weight " << GetStageWeights(countStage) << std::endl;
        		double deltaTimeStage = this->GetStageTimeFactor(countStage)*mTimeStep;
				this->GetStageDerivativeFactor(stageDerivativeFactor, countStage);
				disp_j_tmp = disp_j+vel_j*deltaTimeStage;
				for (int countStage2=0; countStage2<countStage; countStage2++)
				{
					if (stageDerivativeFactor[countStage2]!=0.)
					{
						disp_j_tmp +=acc_j_tmp[countStage2]*(stageDerivativeFactor[countStage2]*mTimeStep*mTimeStep);
					}
				}

				if (this->HasTimeChanged(countStage)==true)
				{
					curTime=prevCurTime+deltaTimeStage;
					mTime=prevTime+deltaTimeStage;

                    //to be implemented for time dependent problems mStructure->SetCurrentTime(mTime);
					//an update of the external load factor and the time dependent constraint is only
					//necessary for a modified global time
					if (mTimeDependentConstraint!=-1)
    				{
				        if (CmatT.GetNumEntries() > 0)
				        {
				            throw MechanicsException("[NuTo::NystroemBase::Solve] solution with constraints not yet implemented.");
				        }
    					double timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
    					mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
    					//mStructure->ConstraintGetRHSAfterGaussElimination(bRHS);
    				}
    				//calculate external force
    				this->CalculateExternalLoad(*mStructure, curTime, extForce_j, extForce_k);
				}

				mStructure->NodeMergeActiveDofValues(0,disp_j_tmp);
				mStructure->ElementTotalUpdateTmpStaticData();

				//calculate internal force (with update of history variables = true)
				mStructure->BuildGlobalGradientInternalPotentialSubVectors(intForce_j,intForce_k,false);

				//std::cout << "norm of deltaForce " << (intForce_j-extForce_j).norm() << std::endl;

				//std::cout << "d_disp_j_tmp " << d_disp_j_tmp[countStage](0) << std::endl;
				if (mUseDiagonalMassMatrix)
				{
					//no system has to be solved
					acc_j_tmp[countStage]  = (invMassMatrix_j.asDiagonal()*(extForce_j-intForce_j));
				}
				else
				{
					//system of equations has to be solved
					mySolver.Solve(fullMassMatrix_jj, extForce_j-intForce_j, acc_j_tmp[countStage]);
				}
				//std::cout << "d_vel_j_tmp " << d_vel_j_tmp[countStage](0) << std::endl;
				//std::cout << "norm of acc " << (d_vel_j_tmp).norm() << std::endl;

				disp_j_new += acc_j_tmp[countStage]*(GetStageWeights1(countStage)*mTimeStep*mTimeStep);
				//std::cout << "disp_j_new " << disp_j_new(0) << std::endl;
				vel_j_new  += acc_j_tmp[countStage]*(GetStageWeights2(countStage)*mTimeStep);
				//std::cout << "vel_j_new " << vel_j_new(0) << std::endl;
            }

        	mTime = prevTime + mTimeStep;
        	curTime = prevCurTime + mTimeStep;

			//std::cout << "final disp_j_new " << disp_j_new(0) << std::endl;
			mStructure->NodeMergeActiveDofValues(0,disp_j_new);
            if (mMergeActiveDofValuesOrder1)
                mStructure->NodeMergeActiveDofValues(1,vel_j_new);
            if (mMergeActiveDofValuesOrder2)
            {
            	mStructure->NodeMergeActiveDofValues(2,(vel_j_new-vel_j)*(1./mTimeStep));
            }

			mStructure->ElementTotalUpdateTmpStaticData();
			mStructure->ElementTotalUpdateStaticData();
			//std::cout << "delta disp between time steps" <<  (disp_j-disp_j_new).norm() << std::endl;

            disp_j = disp_j_new;
            vel_j  = vel_j_new;

            //**********************************************
			//PostProcessing
			//**********************************************
            if (CmatT.GetNumEntries() > 0)
	        {
	            throw MechanicsException("[NuTo::NystroemBase::Solve] not implemented for constrained systems including multiple dofs.");
	        }
	        else
	        {
	        	// outOfBalance_j is automatically zero
			    //outOfBalance_j.Resize(intForce_j.GetNumRows());
	        	//the acceleration of the dofs k is given by the acceleration of the rhs of the constraint equation
	        	//this is calculated using finite differencs
	        	//make sure to recalculate the internal force and external force (if time factor is not 1)
				//if (mTimeDependentConstraint!=-1)
				//{
				//	throw MechanicsException("[NuTo::NystroemBase::Solve3] solution with constraints not yet implemented.");
				//}

				//acc_k = (bRHSprev-bRHShalf*2+bRHSend)*(4./(timeStep*timeStep))
	        	//outOfBalance_k = intForce_k - extForce_k + massMatrix_k.asDiagonal()*acc_k;
	        }

			//postprocess data for plotting
			//, NuTo::FullVector<double,Eigen::Dynamic>(intForce_k-extForce_k)
            this->PostProcess(outOfBalance_j, outOfBalance_k);
        }
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::NystroemBase::Solve] performing Nystroem iteration.");
        throw e;
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mStructure->GetLogger()<<"[NuTo::NystroemBase::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mStructure->GetLogger()<< "[NuTo::NystroemBase::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return NuTo::Error::SUCCESSFUL;

}


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NystroemBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
