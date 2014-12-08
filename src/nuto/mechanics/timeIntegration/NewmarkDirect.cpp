// $Id: NewmarkDirect.cpp 575 2011-09-20 18:05:35Z unger3 $

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
#include "nuto/math/SparseDirectSolverMKLPardiso.h"

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/timeIntegration/ResultGroupNodeForce.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::NewmarkDirect::NewmarkDirect (StructureBase* rStructure)  : NewmarkBase (rStructure)
{
    mMinLineSearchStep = 0.01;
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::NewmarkDirect::Info()const
{
    NewmarkBase::Info();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NewmarkDirect::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NewmarkDirect::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NewmarkDirect::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NewmarkDirect::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NewmarkDirect::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NewmarkDirect::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NewmarkDirect::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of NewmarkDirect" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NewmarkBase)
           & BOOST_SERIALIZATION_NVP(mMinLineSearchStep);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of NewmarkDirect" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


NuTo::Error::eError NuTo::NewmarkDirect::Solve(double rTimeDelta)
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
        if (mMaxTimeStep==0)
            throw MechanicsException("[NuTo::NewmarkDirect::Solve] max time step is set to zero.");

        //renumber dofs and build constraint matrix
        mStructure->NodeBuildGlobalDofs();

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
        if (this->IsDynamic())
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
            //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> m11Full(massMatrix_jj);
            //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenValues;
            //m11Full.EigenValuesSymmetric(eigenValues);
            //std::cout << "eigenvalues mass" << "\n";
            //std::cout << eigenValues.Trans() << "\n";
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
        mStructure->NodeExtractDofValues(0,lastConverged_disp_j, lastConverged_disp_k);
//        mInternalEnergy = mStructure->ElementTotalGetInternalEnergy();
//        FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> plotHistory(1,7);
//        plotHistory(0,0) = mTime;
//        plotHistory(0,2) = mInternalEnergy;
//        plotHistory(0,6) = mInternalEnergy;
        if (this->IsDynamic())
        {
            mStructure->NodeExtractDofValues(1,lastConverged_vel_j,lastConverged_vel_k);
            mStructure->NodeExtractDofValues(2,lastConverged_acc_j,lastConverged_acc_k);

//            mKineticEnergy = 0.5*(lastConverged_vel_j.Dot(massMatrix_jj*lastConverged_vel_j+massMatrix_jk*lastConverged_vel_k)+lastConverged_vel_k.Dot(massMatrix_kj*lastConverged_vel_j+massMatrix_kk*lastConverged_vel_k));
//            mExternalEnergy = 0;
//            mDampedEnergy = 0.;

//            plotHistory(0,3) = mKineticEnergy;
//            plotHistory(0,4) = mDampedEnergy;
//            plotHistory(0,5) = mExternalEnergy;
//            plotHistory(0,6)+= mKineticEnergy ;
        }

        double curTime  = 0;

        // initialize the structure times
        mStructure->SetPrevTime(curTime);
        mStructure->SetTime(curTime);

        //apply constraints for last converged time step
        double timeDependentConstraintFactor(0);
        if (mTimeDependentConstraint!=-1)
        {
            timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
            mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
        }
//        plotHistory(0,1) = timeDependentConstraintFactor;
        mStructure->ConstraintGetRHSAfterGaussElimination(bRHSprev);
        mStructure->NodeMergeActiveDofValues(0,lastConverged_disp_j); //disp_k is internally calculated from the previous constraint matrix
        mStructure->ElementTotalUpdateTmpStaticData();

        //calculate internal force, update of history variables=false
        mStructure->BuildGlobalGradientInternalPotentialSubVectors(prevIntForce_j,prevIntForce_k,false);
        prevResidual_j = prevIntForce_j;
        prevResidual_k = prevIntForce_k;

        if (this->IsDynamic())
        {
            if (mMuDampingMass>0)
            {
                //add damping terme
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
				//add mass terme
				prevResidual_j += (massMatrix_jj*lastConverged_acc_j+massMatrix_jk*lastConverged_acc_k);
				prevResidual_k += (massMatrix_kj*lastConverged_acc_j+massMatrix_kk*lastConverged_acc_k);
        	}
        }
        //add external force
        CalculateExternalLoad(*mStructure, curTime, prevExtForce_j, prevExtForce_k);
        prevResidual_j -= prevExtForce_j;
        prevResidual_k -= prevExtForce_k;

        if (Cmat.GetNumEntries()>0)
        {
            residual_mod=prevResidual_j - CmatT*prevResidual_k;
        }
        else
        {
            residual_mod=prevResidual_j;
        }

        if (residual_mod.Norm()>mToleranceForce)
            throw MechanicsException("[NuTo::NewmarkDirect::Solve] Initial configuration is not in (dynamic) equilibrium.");
        std::cout << "residual in initial configuration " << residual_mod.Norm() << std::endl;

/*        FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> plotVector0(plotHistory);
        for (int countGroup=0; countGroup<mVecGroupNodesReactionForces.GetNumRows(); countGroup++)
        {
            FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> reactionForce(1,mStructure->GetDimension());

            boost::ptr_map<int,GroupBase>::iterator itGroup = mStructure->mGroupMap.find(mVecGroupNodesReactionForces(countGroup,0));
            if (itGroup==mStructure->mGroupMap.end())
                throw MechanicsException("[NuTo::NewmarkDirect::Solve] node group with the given identifier for the reaction forces does not exist.");
            if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
                throw MechanicsException("[NuTo::NewmarkDirect::Solve] Group is not a node group (reaction forces).");
            Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();
            assert(nodeGroup!=0);

            //all nodes have to have the same dimension
            if(nodeGroup->GetNumMembers()<1)
                throw MechanicsException("[NuTo::StructureBase::NodeGroupGetCoordinates] Group has no members.");

            //std::cout << "residual k " << residual_k.transpose() << std::endl;
            for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++)
            {
                for (int countDimension=0; countDimension<mStructure->GetDimension(); countDimension++)
                {
                    int theDof = itNode->second->GetDofDisplacement(countDimension);
                    if (theDof<mStructure->GetNumActiveDofs())
                    {
                        reactionForce(0,countDimension)+=prevResidual_j(theDof);
                    }
                    else
                    {
                        reactionForce(0,countDimension)+=prevResidual_k(theDof-mStructure->GetNumActiveDofs());
                    }
                }
            }
            //std::cout << "reaction force for node group " <<  countGroup << " : " << reactionForce.Trans() << std::endl;
            plotVector0.AppendColumns(reactionForce);
        }
*/
        PostProcess(prevResidual_j, prevResidual_k);

        double timeStep = mTimeStep;
        while (curTime < rTimeDelta)
        {
            //apply constraints for last converged time step
            if (mTimeDependentConstraint!=-1)
            {
				timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
				mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
            }
            mStructure->ConstraintGetRHSAfterGaussElimination(bRHSprev);
            mStructure->NodeMergeActiveDofValues(0,lastConverged_disp_j); //disp_k is internally calculated from the previous constraint matrix
            mStructure->ElementTotalUpdateTmpStaticData();

            //add external force
            CalculateExternalLoad(*mStructure, curTime, prevExtForce_j, prevExtForce_k);

            //increase time step
            curTime += timeStep;

            if (timeStep<mMinTimeStep)
                throw MechanicsException("[NuTo::NewmarkDirect::Solve] time step is smaller than minimum - no convergence is obtained.");


            //if the time difference towards the end is very small, increase the time step, if it exceeds the total time, decrease the time step
            if (mAutomaticTimeStepping)
            {
                if (rTimeDelta-curTime<0.5*timeStep)
                {
                    timeStep += rTimeDelta-curTime;
                    curTime = rTimeDelta;
                }
            }

            // set new structure time at the end of the time increment
            mStructure->SetTime(curTime);

            if (this->IsDynamic())
            {
                if (1.-mGamma*mMuDampingMass*timeStep<=0)
                    throw MechanicsException("[NuTo::NewmarkDirect::Solve] Factor for the mass matrix is negative - reduce your time step.");
            }

            mStructure->NodeMergeActiveDofValues(0,lastConverged_disp_j);
            mStructure->ElementTotalUpdateTmpStaticData();

            //calculate the initial out-of-balance force
            disp_j = lastConverged_disp_j;
            disp_k = lastConverged_disp_k;

            //calculate internal force, update of history variables=false
            mStructure->BuildGlobalGradientInternalPotentialSubVectors(prevIntForce_j,prevIntForce_k,false);
            residual_j = prevIntForce_j;
            residual_k = prevIntForce_k;

            //increase time step by half and calculate the constraint matrix (used to approximate the velocity and acceleration of the rhs)
            if (this->IsDynamic())
            {
                if (mTimeDependentConstraint!=-1)
                {
					timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
					mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
                }
                mStructure->ConstraintGetRHSAfterGaussElimination(bRHShalf);
            }

            //apply constraints for the new time step (modified bRHS)
            if (mTimeDependentConstraint!=-1)
            {
				timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
				mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
            }
            mStructure->ConstraintGetRHSAfterGaussElimination(bRHSend);
            FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deltaBRHS(bRHSend-bRHSprev);

            //std::cout << "bRHSprev " << bRHSprev.Trans() << "\n";
            //std::cout << "bRHShalf " << bRHShalf.Trans() << "\n";
            //std::cout << "bRHSend " << bRHSend.Trans() << "\n";
            //std::cout << "bRHSdot " << bRHSdot.Trans() << "\n";
            //std::cout << "bRHSddot " << bRHSddot.Trans() << "\n";

            if (this->IsDynamic())
            {
                //calculate approximations to the time derivates of the rhs of the constraint matrix
                bRHSdot = (bRHSprev*5.-bRHShalf*8.+bRHSend*3.)*(-1./(timeStep));
                bRHSddot = (bRHSprev-bRHShalf*2+bRHSend)*(4./(timeStep*timeStep));

                //calculate new accelerations and velocities of independent dofs
                acc_j = (lastConverged_vel_j+lastConverged_acc_j*(timeStep*(0.5-mBeta)))*(-1./(timeStep*mBeta));
                vel_j = lastConverged_vel_j+lastConverged_acc_j*((1.-mGamma)*timeStep)+acc_j*(mGamma*timeStep);

                //calculate new accelerations and velocities of dependent dofs
                acc_k = bRHSddot - (Cmat*acc_j);
                vel_k = bRHSdot - (Cmat*vel_j);

				if (mMuDampingMass>0)
				{
					//add damping terme
					if (mUseLumpedMass)
					{
						residual_j += lumped_massMatrix_j.asDiagonal()*vel_j * mMuDampingMass;
						residual_k += lumped_massMatrix_k.asDiagonal()*vel_k * mMuDampingMass;
					}
					else
					{
	                    //add damping terme
	                    residual_j += (massMatrix_jj*vel_j+massMatrix_jk*vel_k) *mMuDampingMass;
	                    residual_k += (massMatrix_kj*vel_j+massMatrix_kk*vel_k) *mMuDampingMass;
					}
				}
				if (mUseLumpedMass)
				{
					residual_j += lumped_massMatrix_j.asDiagonal()*acc_j;
					residual_k += lumped_massMatrix_k.asDiagonal()*acc_k;
				}
				else
				{
	                //add mass terme
	                residual_j += (massMatrix_jj*acc_j+massMatrix_jk*acc_k);
	                residual_k += (massMatrix_kj*acc_j+massMatrix_kk*acc_k);
				}
            }
            //add external force
            CalculateExternalLoad(*mStructure, curTime, extForce_j, extForce_k);
            residual_j -= extForce_j;
            residual_k -= extForce_k;

            //calculate stiffness and disp force vector (the displacements of the dependent dofs do not take into account the modified rhs)
            if (Cmat.GetNumEntries()>0)
            {
                stiffMatrix_jj.SetZeroEntries();stiffMatrix_jk.SetZeroEntries();stiffMatrix_kj.SetZeroEntries();stiffMatrix_kk.SetZeroEntries();
                mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::STIFFNESS, stiffMatrix_jj, stiffMatrix_jk, stiffMatrix_kj, stiffMatrix_kk);
                if (this->IsDynamic())
                {
                    double factor((1.+mGamma*mMuDampingMass*timeStep)/(timeStep*timeStep*mBeta));
    				if (mUseLumpedMass)
    				{
                        stiffMatrix_jj.AddScalDiag(lumped_massMatrix_j,factor);
                        stiffMatrix_kk.AddScalDiag(lumped_massMatrix_k,factor);
    				}
    				else
    				{
                        stiffMatrix_jj.AddScal(massMatrix_jj,factor);
                        stiffMatrix_jk.AddScal(massMatrix_jk,factor);
                        stiffMatrix_kj.AddScal(massMatrix_kj,factor);
                        stiffMatrix_kk.AddScal(massMatrix_kk,factor);
    				}
                }
                //calculate the residual for the zero time step
                residual_j += stiffMatrix_jk * deltaBRHS;
                residual_k += stiffMatrix_kk * deltaBRHS;
                residual_mod=residual_j - CmatT*residual_k;

                //add damping and mass to full hessian
                stiffMatrix_jj -= CmatT * stiffMatrix_kj + stiffMatrix_jk * Cmat;
                stiffMatrix_jj += CmatT * stiffMatrix_kk * Cmat;
            }
            else
            {
                stiffMatrix_jj.SetZeroEntries();stiffMatrix_jk.SetZeroEntries();
                mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::STIFFNESS, stiffMatrix_jj, stiffMatrix_jk);

                //calculate the residual for the zero time step
                residual_j += stiffMatrix_jk * deltaBRHS;
                residual_mod=residual_j;

                //add damping and mass to full hessian
                if (this->IsDynamic())
                {
                    double factor((1.+mGamma*mMuDampingMass*timeStep)/(timeStep*timeStep*mBeta));
    				if (mUseLumpedMass)
    				{
    					stiffMatrix_jj.AddScalDiag(lumped_massMatrix_j,factor);
    				}
    				else
    				{
                        stiffMatrix_jj.AddScal(massMatrix_jj,factor);
                        stiffMatrix_jk.AddScal(massMatrix_jk,factor);
    				}
                }
            }


            if (mCheckCoefficientMatrix)
            {
                mStructure->CheckCoefficientMatrix_0(1.e-6, false);
                mStructure->ElementCheckCoefficientMatrix_0(1.e-6);

//                    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> tmp;
//                    mStructure->ElementCheckCoefficientMatrix_0(1.e-6, 2,tmp ,true);
            }

			//solve for trial state
			NuTo::SparseMatrixCSRGeneral<double> hessianModSolver(stiffMatrix_jj);
//FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> hessianModFull(stiffMatrix_jj);
//std::cout << "stiffness\n" << hessianModFull << std::endl;

			hessianModSolver.SetOneBasedIndexing();
			mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
			delta_disp_j*=-1;

            //calculate trial state
            disp_j += delta_disp_j;
            disp_k = bRHSend - Cmat*disp_j;
            if (this->IsDynamic())
            {
                acc_j  += delta_disp_j * (1./(timeStep*timeStep*mBeta));
                vel_j  += delta_disp_j * (mGamma/(timeStep*mBeta));

                acc_k = bRHSddot - (Cmat*acc_j);
                vel_k = bRHSdot - (Cmat*vel_j);
            }

            //apply displacements
            //std::cout << "norm of delta disp " << delta_disp_j.Norm() << std::endl;
            mStructure->NodeMergeActiveDofValues(0,disp_j);
            mStructure->ElementTotalUpdateTmpStaticData();

            //calculate internal force, update of history variables=false
            mStructure->BuildGlobalGradientInternalPotentialSubVectors(intForce_j, intForce_k,false);
            residual_j = intForce_j;
            residual_k = intForce_k;
            if (this->IsDynamic())
            {
                if (mMuDampingMass>0)
                {
                    if (mUseLumpedMass)
                    {
                    	residual_j += (lumped_massMatrix_j.asDiagonal()*vel_j) *mMuDampingMass;
                    }
                    else
                    {
						//add damping terme
						residual_j += (massMatrix_jj*vel_j+massMatrix_jk*vel_k) *mMuDampingMass;
						residual_k += (massMatrix_kj*vel_j+massMatrix_kk*vel_k) *mMuDampingMass;
                    }
                }

                if (mUseLumpedMass)
                {
					residual_j += lumped_massMatrix_j.asDiagonal()*acc_j;

                }
                else
                {
					//add mass terme
					residual_j += massMatrix_jj*acc_j+massMatrix_jk*acc_k;
					residual_k += massMatrix_kj*acc_j+massMatrix_kk*acc_k;
                }
            }

            //add external force
            residual_j -= extForce_j;
            residual_k -= extForce_k;

            //calculate the initial residual
            if (Cmat.GetNumEntries()>0)
            {
                residual_mod=residual_j - CmatT*residual_k;
            }
            else
            {
                residual_mod=residual_j;
            }

            //calculate norm of initial residual (first time step)
            double normResidual(residual_mod.Norm());
            //std::cout << "max Residual " << residual_mod.Max() << std::endl;

            int iteration(0);
            while(normResidual>mToleranceForce && iteration<mMaxNumIterations)
            {
                //calculate stiffness
                if (Cmat.GetNumEntries()>0)
                {
                    stiffMatrix_jj.SetZeroEntries();stiffMatrix_jk.SetZeroEntries();stiffMatrix_kj.SetZeroEntries();stiffMatrix_kk.SetZeroEntries();
                    mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::STIFFNESS, stiffMatrix_jj, stiffMatrix_jk, stiffMatrix_kj, stiffMatrix_kk);
                    if (this->IsDynamic())
                    {
                        //add damping and mass to full hessian
                        double factor((1.+mGamma*mMuDampingMass*timeStep)/(timeStep*timeStep*mBeta));
                        if (mUseLumpedMass)
                        {
                        	stiffMatrix_jj.AddScalDiag(lumped_massMatrix_j,factor);
                        }
                        else
                        {
                            stiffMatrix_jj.AddScal(massMatrix_jj,factor);
                            stiffMatrix_jk.AddScal(massMatrix_jk,factor);
                        }
                    }

                    stiffMatrix_jj -= CmatT * stiffMatrix_kj + stiffMatrix_jk * Cmat;
                    stiffMatrix_jj += CmatT * stiffMatrix_kk * Cmat;
                }
                else
                {
                    stiffMatrix_jj.SetZeroEntries();stiffMatrix_jk.SetZeroEntries();
                    mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::STIFFNESS, stiffMatrix_jj, stiffMatrix_jk);

                    if (this->IsDynamic())
                    {
                        //add damping and mass to full hessian
                        double factor((1.+mGamma*mMuDampingMass*timeStep)/(timeStep*timeStep*mBeta));
                        if (mUseLumpedMass)
                        {
                        	stiffMatrix_jj.AddScalDiag(lumped_massMatrix_j,factor);
                        }
                        else
                        {
							stiffMatrix_jj.AddScal(massMatrix_jj,factor);
							stiffMatrix_jk.AddScal(massMatrix_jk,factor);
                        }
                    }
                }

                if (mCheckCoefficientMatrix)
                {
                    mStructure->CheckCoefficientMatrix_0(1.e-6, false);
                    mStructure->ElementCheckCoefficientMatrix_0(1.e-6);

//                    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> tmp;
//                    mStructure->ElementCheckCoefficientMatrix_0(1.e-6, 2,tmp ,true);
                }

/*
//check stiffness matrix
stiffMatrix_jj.SetZeroEntries();stiffMatrix_jk.SetZeroEntries();
mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::STIFFNESS, stiffMatrix_jj, stiffMatrix_jk);
//std::cout << "stiffness exact1\n" <<  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>(stiffMatrix_jj) << std::endl;
FullVector<double, Eigen::Dynamic> check_disp_j1,check_disp_j2,check_disp_k1,check_disp_k2;
mStructure->NodeExtractDofValues(0,check_disp_j1,check_disp_k1);
//std::cout << "active dof values " << check_disp_j1 << std::endl;
FullVector<double, Eigen::Dynamic> check_intforce_j1(check_disp_j1),check_intforce_k1(check_disp_k1),
		check_intforce_j2(check_disp_j1),check_intforce_k2(check_disp_k1);
FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> stiffMatrix_jj_cdf(stiffMatrix_jj);
stiffMatrix_jj_cdf.setZero();

mStructure->BuildGlobalGradientInternalPotentialSubVectors(check_intforce_j1,check_intforce_k1);
//std::cout << "check_intforce_j1 \n" << check_intforce_j1 << std::endl;
double delta(1e-8);
for (int count=0; count<check_disp_j1.GetNumRows(); count++)
{
	check_disp_j2 = check_disp_j1;
	check_disp_j2(count)+=delta;
	//mStructure->NodeMergeActiveDofValues(0,check_disp_j2);
	mStructure->NodeMergeDofValues(0,check_disp_j2,check_disp_k1);
	mStructure->ElementTotalUpdateTmpStaticData();
	mStructure->BuildGlobalGradientInternalPotentialSubVectors(check_intforce_j2,check_intforce_k2);
	//std::cout << "check_intforce_j2 \n" << check_intforce_j2 << std::endl;
	stiffMatrix_jj_cdf.col(count) = ((check_intforce_j2-check_intforce_j1).transpose())*(1./delta);
}
if ((stiffMatrix_jj_cdf-FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>(stiffMatrix_jj)).cwiseAbs().maxCoeff()>1e-3)
{
	std::cout << "stiffness exact2\n" <<  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>(stiffMatrix_jj) << std::endl;
	std::cout << "stiffness cdf\n" <<  stiffMatrix_jj_cdf << std::endl;
	std::cout << "delta \n" <<  stiffMatrix_jj_cdf-FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>(stiffMatrix_jj) << std::endl;
}
else
{
	std::cout << "error stiffness is " << (stiffMatrix_jj_cdf-FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>(stiffMatrix_jj)).cwiseAbs().maxCoeff() << std::endl;
}
mStructure->NodeMergeDofValues(0,check_disp_j1,check_disp_k1);
*/

                //solve for new state
                NuTo::SparseMatrixCSRGeneral<double> hessianModSolver(stiffMatrix_jj);
                hessianModSolver.SetOneBasedIndexing();
                mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
                delta_disp_j*=-1;

                //std::cout << "norm of delta disp " << delta_disp_j.Norm() <<  std::endl;
                //std::cout << "delta disp " << std::endl;
                //std::cout << delta_disp_j.Trans() << std::endl;

                //perform a line search
                double alpha(1);
                double trialNormResidual(0);
                do
                {
                     //calculate trial state
                    trial_disp_j = disp_j + delta_disp_j * alpha;
                    delta_disp_k = (Cmat*delta_disp_j)*(-1.);
                    trial_disp_k = disp_k + delta_disp_k * alpha;

                    if (this->IsDynamic())
                    {
                        trial_acc_j  = acc_j+delta_disp_j * (alpha/(timeStep*timeStep*mBeta));
                        trial_acc_k  = acc_k+delta_disp_k * (alpha/(timeStep*timeStep*mBeta));

                        trial_vel_j  = vel_j+delta_disp_j * (alpha*mGamma/(timeStep*mBeta));
                        trial_vel_k  = vel_k+delta_disp_k * (alpha*mGamma/(timeStep*mBeta));
                    }

                    //apply displacements
                    mStructure->NodeMergeActiveDofValues(0,trial_disp_j);
                    mStructure->ElementTotalUpdateTmpStaticData();

                    //calculate internal force, update of history variables=false
                    mStructure->BuildGlobalGradientInternalPotentialSubVectors(intForce_j, intForce_k,false);
                    residual_j = intForce_j;
                    residual_k = intForce_k;

                    if (this->IsDynamic())
                    {
                        if (mMuDampingMass>0)
                        {
                            if (mUseLumpedMass)
                            {
								residual_j += (lumped_massMatrix_j.asDiagonal()*trial_vel_j) *mMuDampingMass;
								residual_k += (lumped_massMatrix_k.asDiagonal()*trial_vel_k) *mMuDampingMass;
                            }
                            else
                            {
								//add damping terme
								residual_j += (massMatrix_jj*trial_vel_j+massMatrix_jk*trial_vel_k) *mMuDampingMass;
								residual_k += (massMatrix_jk.TransMult(trial_vel_j)+massMatrix_kk*trial_vel_k) *mMuDampingMass;
                            }
                        }

                        //add mass terme
                        if (mUseLumpedMass)
                        {
							residual_j += lumped_massMatrix_j.asDiagonal()*trial_acc_j;
							residual_k += lumped_massMatrix_k.asDiagonal()*trial_acc_k;
                        }
                        else
                        {
                            residual_j += (massMatrix_jj*trial_acc_j+massMatrix_jk*trial_acc_k);
                            residual_k += (massMatrix_jk.TransMult(trial_acc_j)+massMatrix_kk*trial_acc_k);
                        }
                    }

                    //add external force (assumed to be independent of the deformation)
                    residual_j -=extForce_j;
                    residual_k -=extForce_k;

                    //std::cout << "residual_j\n " << residual_j << std::endl;
                    //std::cout << "residual_k\n " << residual_k << std::endl;

                    //calculate the initial residual
                    if (Cmat.GetNumEntries()>0)
                    {
                        residual_mod = residual_j - CmatT*residual_k;
                    }
                    else
                    {
                        residual_mod = residual_j;
                    }
                    //std::cout << "residual_mod\n " << residual_mod << std::endl;

                    trialNormResidual=residual_mod.Norm();

                    mStructure->GetLogger() << "  linesearch alpha " << alpha <<" previous residual " << normResidual << " trial residual " << trialNormResidual <<  "\n";

                    alpha*=0.5;
                }
                while(alpha>mMinLineSearchStep && trialNormResidual>(1.-alpha)*normResidual);
                   //in the first iteration, there is no line search

                if (alpha>mMinLineSearchStep)
                {
                    //improvement is achieved, go to next Newton step
                    disp_j = trial_disp_j;
                    disp_k = trial_disp_k;

                    if (this->IsDynamic())
                    {
                        acc_j  = trial_acc_j;
                        acc_k  = trial_acc_k;

                        vel_j  = trial_vel_j;
                        vel_k  = trial_vel_k;
                    }

                    normResidual = trialNormResidual;

                    iteration++;
                }
                else
                {
                    //and leave
                    iteration = mMaxNumIterations;
                }
            } // end of while(normResidual<mToleranceForce && iteration<mMaxNumIterations)

            if (normResidual<=mToleranceForce)
            {
                //converged solution
                mStructure->GetLogger() << " *** UpdateStaticData *** from NewMarkDirect \n";
                mStructure->ElementTotalUpdateStaticData();

                //calculate energies (at this point, use the trapezoidal rule for all energies, for some energies, you can compare that to the direct calculation
                // e.g. internal energy for elastic materials or the kinetic energy
                //double internalEnergyExact = mStructure->ElementTotalGetInternalEnergy();
                //FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> delta_disp_j(disp_j-lastConverged_disp_j);
                //FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> delta_disp_k(disp_k-lastConverged_disp_k);

                //mInternalEnergy+= 0.5*(delta_disp_j.Dot(intForce_j+prevIntForce_j)+
                //                       delta_disp_k.Dot(intForce_k+prevIntForce_k));

                //if (this->IsDynamic())
                //{
                    //double kineticEnergyExact = 0.5*(vel_j.Dot(massMatrix_jj*vel_j+massMatrix_jk*vel_k)+vel_k.Dot(massMatrix_kj*vel_j+massMatrix_kk*vel_k));
                //    mKineticEnergy+= 0.5*(delta_disp_j.Dot(massMatrix_jj*(acc_j+lastConverged_acc_j)+massMatrix_jk*(acc_k+lastConverged_acc_k))+
                //                      delta_disp_k.Dot(massMatrix_kj*(acc_j+lastConverged_acc_j)+massMatrix_kk*(acc_k+lastConverged_acc_k)));

                //    if (mMuDampingMass>0)
                //    {
                //        mDampedEnergy+=mMuDampingMass*0.5*(delta_disp_j.Dot(massMatrix_jj*(vel_j+lastConverged_vel_j)+massMatrix_jk*(vel_k+lastConverged_vel_k))+
                //                                           delta_disp_k.Dot(massMatrix_kj*(vel_j+lastConverged_vel_j)+massMatrix_kk*(vel_k +lastConverged_vel_k)));
                //    }
                //}

                //mExternalEnergy += 0.5*(delta_disp_j.Dot(extForce_j+prevExtForce_j) + delta_disp_k.Dot(extForce_k+prevExtForce_k)+
                //                        delta_disp_j.Dot(residual_j+prevResidual_j) + delta_disp_k.Dot(residual_k+prevResidual_k));

                //this is should be moved to the general postprocessing routine
                //FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> plotVector(1,7);
                //plotVector(0,0) = mTime+timeStep;
                //plotVector(0,1) = timeDependentConstraintFactor;
                //plotVector(0,2) = mInternalEnergy;
                //plotVector(0,3) = mKineticEnergy;
                //plotVector(0,4) = mDampedEnergy;
                //plotVector(0,5) = mExternalEnergy;
                //plotVector(0,6) = mInternalEnergy + mKineticEnergy + mDampedEnergy - mExternalEnergy;

                //plotHistory.AppendRows(plotVector);

                //std::cout << "mInternalEnergy " << mInternalEnergy << "(exact " << internalEnergyExact << "), kinetic energy " << mKineticEnergy << "(exact " << kineticEnergyExact << "), damped energy " << mDampedEnergy ;
                //std::cout << ", external energy " << mExternalEnergy << "\n";
                //std::cout << "energy stored in the system (internal + kinetic) " << mInternalEnergy+mKineticEnergy << "\n";
                //std::cout << "energy removed from the system (damping - external) " << mDampedEnergy-mExternalEnergy << "\n";
                //std::cout << "total energy " << mInternalEnergy+mKineticEnergy + mDampedEnergy-mExternalEnergy << "\n";


                //store converged step
                lastConverged_disp_j = disp_j;
                lastConverged_disp_k = disp_k;
                if (this->IsDynamic())
                {
                    lastConverged_vel_j = vel_j;
                    lastConverged_vel_k = vel_k;
                    lastConverged_acc_j = acc_j;
                    lastConverged_acc_k = acc_k;
                }

                prevResidual_j = residual_j;
                prevResidual_k = residual_k;

                //update nodal data
                mStructure->NodeMergeActiveDofValues(0,disp_j);
                if (this->IsDynamic())
                {
                    mStructure->NodeMergeDofValues(1,vel_j,vel_k);
                    mStructure->NodeMergeDofValues(2,acc_j,acc_k);
                }
                mStructure->ElementTotalUpdateTmpStaticData();

                //update structure time
                mStructure->SetPrevTime(curTime);

                //postprocess
                mTime+=timeStep;

                mStructure->GetLogger() << "Convergence after " << iteration << " iterations at time " << mTime << "(timestep " << timeStep << ").\n";
                //std::cout << "plot Vector " << plotVector << std::endl;

				//perform Postprocessing
                mStructure->GetLogger() << " *** PostProcess *** from NewMarkDirect \n";
				PostProcess(prevResidual_j, prevResidual_k);

                //eventually increase next time step
                if (mAutomaticTimeStepping && iteration<0.25*mMaxNumIterations)
                {
                    timeStep*=1.5;
                    if (timeStep>mMaxTimeStep)
                        timeStep = mMaxTimeStep;
                }
            }
            else
            {
                mStructure->GetLogger() << "No convergence with timestep " << timeStep << "\n";
                //no convergence
                if (mAutomaticTimeStepping)
                {
                    //no convergence, reduce the time step and start from scratch
                    curTime -= timeStep;
                    timeStep*=0.5;
                }
                else
                {
                    throw MechanicsException("[NuTo::NewmarkDirect::Solve] no convergence with the current maximum number of iterations, either use automatic time stepping, reduce the time step or the minimal line search cut back factor.");
                }
            }
        }
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::NewmarkDirect::Solve] performing Newton-Raphson iteration.");
        throw e;
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mStructure->GetLogger()<<"[NuTo::NewmarkDirect::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mStructure->GetLogger()<< "[NuTo::NewmarkDirect::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return NuTo::Error::SUCCESSFUL;

}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::NewmarkDirect::GetTypeId()const
{
    return std::string("NewmarkDirect");
}


#ifdef ENABLE_SERIALIZATION
//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
//! @brief ... save the object to a file
void NuTo::NewmarkDirect::Restore (const std::string &filename, std::string rType )
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
                throw MechanicsException ( "[NewmarkDirect::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[NewmarkDirect::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[NewmarkDirect::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
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
        throw MechanicsException ( "[NewmarkDirect::Restore]Unhandled exception." );
    }
}

//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::NewmarkDirect::Save (const std::string &filename, std::string rType )const
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
            throw MechanicsException ( "[NewmarkDirect::Save]File type not implemented." );
        }
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NewmarkDirect::Save]File save exception in boost - " ) +std::string ( e.what() ) );
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
        throw MechanicsException ( "[NewmarkDirect::Save]Unhandled exception." );
    }
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NewmarkDirect)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
