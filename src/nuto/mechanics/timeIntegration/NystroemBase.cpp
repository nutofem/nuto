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

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/NystroemBase.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::NystroemBase::NystroemBase ()  : TimeIntegrationBase ()
{
    mTimeStep = 0.;
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
           & BOOST_SERIALIZATION_NVP(mTimeStep);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of NystroemBase" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
NuTo::Error::eError NuTo::NystroemBase::Solve(StructureBase& rStructure, double rTimeDelta)
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
            mTimeStep = rStructure.ElementTotalCalculateCriticalTimeStep();

        std::cout << "modify computation of critical time step to include the dependence on the time integration scheme." <<std::endl;
        //calculate instead the smallest eigenfrequency, depending on the time integration this gives the critical time step

        std::cout << "time step " << mTimeStep << std::endl;
        std::cout << "number of time steps " << rTimeDelta/mTimeStep << std::endl;

        //calculate the node ptrs for the postprocessing routines
        CalculateOutputDispNodesPtr(rStructure);

        //renumber dofs and build constraint matrix
        rStructure.NodeBuildGlobalDofs();

        //calculate constraint matrix
        NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
        rStructure.ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
        NuTo::SparseMatrixCSRVector2General<double> Cmat(CmatTmp);
        SparseMatrixCSRVector2General<double> CmatT (Cmat.Transpose());
        FullVector<double,Eigen::Dynamic> bRHS;
        if (CmatT.GetNumEntries() > 0)
        {
            throw MechanicsException("[NuTo::NystroemBase::Solve] not implemented for constrained systems including multiple dofs.");
        }

        //calculate individual inverse mass matrix, use only lumped mass matrices - stored as fullvectors and then use asDiagonal()
        NuTo::FullVector<double,Eigen::Dynamic> invMassMatrix_j(rStructure.GetNumActiveDofs());

        //extract displacements, velocities and accelerations
        NuTo::FullVector<double,Eigen::Dynamic> disp_j, vel_j, tmp_k,
                                                disp_j_tmp,              //temporary argument for the evalution of f(y)
                                                disp_j_new, vel_j_new;   //new d and v at end of time step
        std::vector<NuTo::FullVector<double,Eigen::Dynamic> > acc_j_tmp; //intermediate values of the time derivatives (Z or l)

        NuTo::FullVector<double,Eigen::Dynamic> extForce_j, extForce_k;
        NuTo::FullVector<double,Eigen::Dynamic> intForce_j(rStructure.GetNumActiveDofs()),
        		                                intForce_k(rStructure.GetNumDofs() - rStructure.GetNumActiveDofs());

        //store last converged displacements, velocities and accelerations
        rStructure.NodeExtractDofValues(0,disp_j,tmp_k);
        rStructure.NodeExtractDofValues(1,vel_j,tmp_k);

        //calculate lumped mass matrix
        //intForce_j is just used as a tmp variable
        rStructure.BuildGlobalLumpedHession2(intForce_j,tmp_k);

        //check the sum of all entries
        std::cout << "the total mass is " << intForce_j.sum()/3. +  tmp_k.sum()/3. << std::endl;

        //invert the mass matrix
        invMassMatrix_j = intForce_j.cwiseInverse();

        double curTime  = 0;
		CalculateExternalLoad(rStructure, curTime, extForce_j, extForce_k);
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

//to be implemented rStructure.SetCurrentTime(mTime);
					//an update of the external load factor and the time dependent constraint is only
					//necessary for a modified global time
					if (mTimeDependentConstraint!=-1)
    				{
    					throw MechanicsException("[NuTo::NystroemBase::Solve] solution with constraints not yet implemented.");
    					//double timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
    					//rStructure.ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
    					//rStructure.ConstraintGetRHSAfterGaussElimination(bRHS);
    				}
    				//calculate external force
    				this->CalculateExternalLoad(rStructure, curTime, extForce_j, extForce_k);
				}

				rStructure.NodeMergeActiveDofValues(0,disp_j_tmp);
				rStructure.ElementTotalUpdateTmpStaticData();

				//calculate internal force (with update of history variables = true)
				rStructure.BuildGlobalGradientInternalPotentialSubVectors(intForce_j,intForce_k,false);

				//std::cout << "norm of deltaForce " << (intForce_j-extForce_j).norm() << std::endl;

				//std::cout << "d_disp_j_tmp " << d_disp_j_tmp[countStage](0) << std::endl;
				acc_j_tmp[countStage]  = (invMassMatrix_j.asDiagonal()*(extForce_j-intForce_j));
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
			rStructure.NodeMergeActiveDofValues(0,disp_j_new);
			rStructure.ElementTotalUpdateTmpStaticData();
			rStructure.ElementTotalUpdateStaticData();
			//std::cout << "delta disp between time steps" <<  (disp_j-disp_j_new).norm() << std::endl;
            disp_j = disp_j_new;
            vel_j  = vel_j_new;

			//**********************************************
			//PostProcessing
			//**********************************************
            FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> plotVector(1,2);
			plotVector(0,0) = mTime;

			//temporarily change this
			//plotVector(0,1) = timeDependentConstraintFactor;

			for (unsigned int countNode=0; countNode<mVecOutputDispNodesPtr.size(); countNode++)
			{
				switch(mVecOutputDispNodesPtr[countNode] -> GetNumDisplacements())
				{
				case 1:
				{
					FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> displacements(1,1);
					mVecOutputDispNodesPtr[countNode]->GetDisplacements1D(displacements.data());
					plotVector.AppendColumns(displacements);
				}
				break;
				case 2:
				{
					FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> displacements(1,2);
					mVecOutputDispNodesPtr[countNode]->GetDisplacements2D(displacements.data());
					plotVector.AppendColumns(displacements);
				}
				break;
				case 3:
				{
					FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> displacements(1,3);
					mVecOutputDispNodesPtr[countNode]->GetDisplacements3D(displacements.data());
					plotVector.AppendColumns(displacements);
				}
				break;
				default:
					break;
				}
			}

			//postprocess data for plotting
            this->PostProcess(rStructure, plotVector);
        }
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::NystroemBase::Solve] performing Newton-Raphson iteration.");
        throw e;
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        rStructure.GetLogger()<<"[NuTo::NystroemBase::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        rStructure.GetLogger()<< "[NuTo::NystroemBase::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return NuTo::Error::SUCCESSFUL;

}


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NystroemBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
