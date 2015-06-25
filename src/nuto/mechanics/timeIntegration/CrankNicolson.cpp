// $Id: CrankNicolson.cpp 575 2011-09-20 18:05:35Z unger3 $

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
#include "nuto/mechanics/timeIntegration/CrankNicolson.h"
#include "nuto/mechanics/timeIntegration/ResultGroupNodeForce.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"

#include "nuto/mechanics/elements/ElementBase.h"  // delete me
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h" // delete me
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataDamageViscoPlasticity3D.h" // delete me


//! @brief constructor
//! @param mDimension number of nodes
NuTo::CrankNicolson::CrankNicolson (StructureBase* rStructure)  : TimeIntegrationBase (rStructure)
{
    mMinLineSearchStep = 0.01;
    mVisualizeResidualTimeStep = 0;
    mPerformLineSearch = true;
	mToleranceForce = 1e-6;
	mMaxNumIterations = 20;
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::CrankNicolson::Info()const
{
    TimeIntegrationBase::Info();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::CrankNicolson::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::CrankNicolson::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::CrankNicolson::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::CrankNicolson::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::CrankNicolson::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::CrankNicolson::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::CrankNicolson::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of CrankNicolson" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase)
           & BOOST_SERIALIZATION_NVP(mMinLineSearchStep);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of CrankNicolson" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


NuTo::Error::eError NuTo::CrankNicolson::Solve(double rFinalTime)
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
        std::cout << "start Crank Nicolson solver " << std::endl;
    	if (mMaxTimeStep==0)
            throw MechanicsException("[NuTo::CrankNicolson::Solve] max time step is set to zero.");

        //renumber dofs and build constraint matrix
        mStructure->NodeBuildGlobalDofs();

        //calculate constraint matrix
        NuTo::SparseMatrixCSRGeneral<double> CmatTmp;
        mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
        NuTo::SparseMatrixCSRVector2General<double> Cmat(CmatTmp);
        SparseMatrixCSRVector2General<double> CmatT (Cmat.Transpose());
        int numEntriesCMat(Cmat.GetNumEntries());
        FullVector<double,Eigen::Dynamic> bRHSprev, bRHS;


        //Hessian matrix
        SparseMatrixCSRVector2General<double> hessian_jj(mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());

        // allocate space for stiffness matrix
        SparseMatrixCSRVector2General<double> stiffMatrix_jj(mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
        SparseMatrixCSRVector2General<double> stiffMatrix_jk(mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
        SparseMatrixCSRVector2General<double> stiffMatrix_kj(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumActiveDofs());
        SparseMatrixCSRVector2General<double> stiffMatrix_kk(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(), mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());


        //calculate individual mass matrix (I was lazy here and have just used general matrices, but symmetric is certainly better here)
        NuTo::SparseMatrixCSRVector2General<double> dampingMatrix_jj,dampingMatrix_jk,dampingMatrix_kj,dampingMatrix_kk;
        NuTo::SparseMatrixCSRVector2General<double> massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk;

        //allocate matrices of higher order (1 = damping, 2 = mass for orders larger than 1)
        if (mStructure->GetNumTimeDerivatives()>0)
        {
            dampingMatrix_jj.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
            dampingMatrix_jk.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
            dampingMatrix_kj.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
            dampingMatrix_kk.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
            mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::DAMPING,dampingMatrix_jj,dampingMatrix_jk,dampingMatrix_kj,dampingMatrix_kk);
            if (mStructure->GetNumTimeDerivatives())
            {
                massMatrix_jj.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
                massMatrix_jk.Resize(mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
                massMatrix_kj.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumActiveDofs());
                massMatrix_kk.Resize(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
                mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::MASS,massMatrix_jj,massMatrix_jk,massMatrix_kj,massMatrix_kk);
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
        mStructure->NodeExtractDofValues(0,lastConverged_disp_j, lastConverged_disp_k);

        if (mStructure->GetNumTimeDerivatives())
        {
            mStructure->NodeExtractDofValues(1,lastConverged_vel_j,lastConverged_vel_k);
            if (mStructure->GetNumTimeDerivatives()>1)
                mStructure->NodeExtractDofValues(2,lastConverged_acc_j,lastConverged_acc_k);
        }
        //this is the current starting time of the simulation (think about several load cycles calculated before that)
        double curTime(mTime);

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

        if (mStructure->GetNumTimeDerivatives()>0)
        {
            prevResidual_j += dampingMatrix_jj*lastConverged_vel_j+dampingMatrix_jk*lastConverged_vel_k;
            prevResidual_k += dampingMatrix_kj*lastConverged_vel_j+dampingMatrix_kk*lastConverged_vel_k;
            if (mStructure->GetNumTimeDerivatives()>1)
            {
                prevResidual_j += (massMatrix_jj*lastConverged_acc_j+massMatrix_jk*lastConverged_acc_k);
                prevResidual_k += (massMatrix_kj*lastConverged_acc_j+massMatrix_kk*lastConverged_acc_k);
            }
        }

        //add external force
        CalculateExternalLoad(*mStructure, curTime, prevExtForce_j, prevExtForce_k);
        prevResidual_j -= prevExtForce_j;
        prevResidual_k -= prevExtForce_k;

        residual_mod=prevResidual_j - CmatT*prevResidual_k;

//@COMMENT check, if this could be avoided, only required if check for initial equilibrium required
        //std::cout << "residual in initial configuration " << residual_mod.Norm() << std::endl;
        if (residual_mod.Norm()>mToleranceForce)
            throw MechanicsException("[NuTo::CrankNicolson::Solve] Initial configuration is not in (dynamic) equilibrium.");

        PostProcess(prevResidual_j, prevResidual_k);

        double timeStep = mTimeStep;
        while (curTime < rFinalTime)
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

            //add previous external force (because delta is required)
            CalculateExternalLoad(*mStructure, curTime, prevExtForce_j, prevExtForce_k);

            //increase time step
            curTime += timeStep;

            if (timeStep<mMinTimeStep)
                throw MechanicsException("[NuTo::CrankNicolson::Solve] time step is smaller than minimum - no convergence is obtained.");

            //check if curTime is close to final time
            if (curTime>rFinalTime || (rFinalTime-curTime)<0.2*timeStep)
            {
                timeStep+=rFinalTime-curTime;
                curTime = rFinalTime;
            }

            // set new structure time at the end of the new time increment
            mStructure->SetTime(curTime);

            //add residual contribution from external force
            CalculateExternalLoad(*mStructure, curTime, extForce_j, extForce_k);
            residual_j = prevExtForce_j-extForce_j;
            residual_k = prevExtForce_k-extForce_k;

            //apply constraints for the new time step (modified bRHS)
            //remember that this does not change the nodal values (at the nodes)
            if (mTimeDependentConstraint!=-1)
            {
                timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
                mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
            }
            mStructure->ConstraintGetRHSAfterGaussElimination(bRHS);
            FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deltaBRHS(bRHS-bRHSprev);

            //calculate the residual contribution from stiffness
            stiffMatrix_jj.SetZeroEntries();stiffMatrix_jk.SetZeroEntries();stiffMatrix_kj.SetZeroEntries();stiffMatrix_kk.SetZeroEntries();
            mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::STIFFNESS, stiffMatrix_jj, stiffMatrix_jk, stiffMatrix_kj, stiffMatrix_kk);
            residual_j += stiffMatrix_jk * deltaBRHS;
            if (mStructure->GetNumTimeDerivatives()>0)
            {
                dampingMatrix_jj.SetZeroEntries();dampingMatrix_jk.SetZeroEntries();dampingMatrix_kj.SetZeroEntries();dampingMatrix_kk.SetZeroEntries();
                mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::DAMPING, dampingMatrix_jj, dampingMatrix_jk, dampingMatrix_kj, dampingMatrix_kk);
                residual_j -= (dampingMatrix_jj*lastConverged_vel_j + (dampingMatrix_jk * (deltaBRHS *(-1./timeStep)+lastConverged_vel_k)))*2.;
                if (mStructure->GetNumTimeDerivatives()>1)
                {
                    massMatrix_jj.SetZeroEntries();massMatrix_jk.SetZeroEntries();massMatrix_kj.SetZeroEntries();massMatrix_kk.SetZeroEntries();
                    mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::MASS, massMatrix_jj, massMatrix_jk, massMatrix_kj, massMatrix_kk);
                    residual_j -= massMatrix_jj * (lastConverged_vel_j*(4./timeStep)+lastConverged_acc_j*2.)+
                                  massMatrix_jk * (deltaBRHS *(-4./(timeStep*timeStep))+lastConverged_vel_k*(4./timeStep)+lastConverged_acc_k*2.);
                }
            }

            if (numEntriesCMat==0)
            {
                residual_mod = residual_j;
            }
            else
            {
                //add residual contribution from external force
                residual_k = prevExtForce_k-extForce_k;

                //calculate the residual contribution from stiffness
                residual_k += stiffMatrix_kk * deltaBRHS;

                //add residual contribution from damping
                residual_k -= (dampingMatrix_kj*lastConverged_vel_j + (dampingMatrix_kk * (deltaBRHS) *(-1./timeStep)+lastConverged_vel_k))*2.;

                //add residual contribution from mass
                residual_k -= massMatrix_kj * (lastConverged_vel_j*(4./timeStep)+lastConverged_acc_j*2.)+
                          massMatrix_kk * (deltaBRHS *(-4./(timeStep*timeStep))+lastConverged_vel_k*(4./timeStep)+lastConverged_acc_k*2.);

                //add residual contribution from dependent dofs
                residual_mod = residual_j - CmatT*residual_k;
            }
//std::cout << "expected rhs for trial state\n" << residual_mod.transpose() << std::endl<< std::endl;
            std::cout << "norm of predicted residual in trial state:" << residual_mod.Norm() << std::endl<< std::endl;

            //build hessian for trial state
            hessian_jj = stiffMatrix_jj;
            if (numEntriesCMat>0)
            {
                hessian_jj += CmatT * stiffMatrix_kk * Cmat - CmatT * stiffMatrix_kj + stiffMatrix_jk * Cmat;
            }
            if (mStructure->GetNumTimeDerivatives()>0)
            {
                hessian_jj+=dampingMatrix_jj*(2./(timeStep));
                if (numEntriesCMat>0)
                {
                    hessian_jj+=(CmatT * dampingMatrix_kk * Cmat - CmatT * dampingMatrix_kj + dampingMatrix_jk * Cmat )*(2./(timeStep));
                }
                if (mStructure->GetNumTimeDerivatives()>1)
                {
                    hessian_jj+=massMatrix_jj*(4./(timeStep*timeStep));
                    if (numEntriesCMat>0)
                    {
                        hessian_jj+=(CmatT * massMatrix_kk * Cmat - CmatT * massMatrix_kj + massMatrix_jk * Cmat )*(4./(timeStep*timeStep));
                    }
                }
            }

            //solve for trial state
            NuTo::SparseMatrixCSRGeneral<double> hessianModSolver(hessian_jj);
            hessianModSolver.SetOneBasedIndexing();
            mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
            delta_disp_j*=-1;

//std::cout << "delta disp prediction\n" << delta_disp_j.transpose() << std::endl<< std::endl;

            //calculate trial state
            disp_j = lastConverged_disp_j + delta_disp_j;
//std::cout << "disp prediction\n" << disp_j.transpose() << std::endl<< std::endl;
            disp_k = bRHS - Cmat*disp_j;
            if (mStructure->GetNumTimeDerivatives()>0)
            {
                delta_disp_k = disp_k - lastConverged_disp_k;
                vel_j  = delta_disp_j * (2./timeStep) - lastConverged_vel_j;
                vel_k  = delta_disp_k * (2./timeStep) - lastConverged_vel_k;

                if (mStructure->GetNumTimeDerivatives()>1)
                {
                    acc_j = delta_disp_j * (4./(timeStep*timeStep)) - lastConverged_vel_j * (4./(timeStep)) - lastConverged_acc_j;
                    acc_k = delta_disp_k * (4./(timeStep*timeStep)) - lastConverged_vel_k * (4./(timeStep)) - lastConverged_acc_k;
                }
//std::cout << "acc_j prediction\n" << acc_j.transpose() << std::endl<< std::endl;
                //std::cout << "lastConverged_acc_j prediction\n" << lastConverged_acc_j.transpose() << std::endl<< std::endl;
                //std::cout << "4./(timeStep*timeStep)\n" << 4./(timeStep*timeStep) << std::endl<< std::endl;
            }

            //apply displacements
            //std::cout << "norm of delta disp " << delta_disp_j.Norm() << std::endl;
            mStructure->NodeMergeActiveDofValues(0,disp_j);
            mStructure->ElementTotalUpdateTmpStaticData();

            //calculate internal force, update of history variables=false
            mStructure->BuildGlobalGradientInternalPotentialSubVectors(intForce_j, intForce_k,false);
            residual_j = intForce_j;
            residual_k = intForce_k;
//std::cout << "intForce_j prediction\n" << intForce_j.transpose() << std::endl<< std::endl;

            if (mStructure->GetNumTimeDerivatives()>0)
            {
                residual_j += dampingMatrix_jj*vel_j+dampingMatrix_jk*vel_k;
                residual_k += dampingMatrix_kj*vel_j+dampingMatrix_kk*vel_k;
//std::cout << "damping j prediction\n" << (dampingMatrix_jj*vel_j+dampingMatrix_jk*vel_k).transpose() << std::endl<< std::endl;
                if (mStructure->GetNumTimeDerivatives()>1)
                {
                    residual_j += (massMatrix_jj*acc_j+massMatrix_jk*acc_k);
                    residual_k += (massMatrix_kj*acc_j+massMatrix_kk*acc_k);
//std::cout << "mass j prediction\n" << (massMatrix_jj*acc_j+massMatrix_jk*acc_k).transpose() << std::endl<< std::endl;
                }
            }

            //add external force
            residual_j -= extForce_j;
            residual_k -= extForce_k;
//std::cout << "extForce_j prediction\n" << extForce_j.transpose() << std::endl<< std::endl;

            //calculate the initial residual
            if (numEntriesCMat>0)
            {
                residual_mod=residual_j - CmatT*residual_k;
            }
            else
            {
                residual_mod=residual_j;
            }

            //calculate norm of initial residual (first time step)
            double normResidual(residual_mod.Norm());
//std::cout << "residual after prediction state\n" << residual_mod.transpose() << std::endl<< std::endl;
            std::cout << "norm of residual after prediction state:" << normResidual << std::endl<< std::endl;

            int iteration(0);
            while(normResidual>mToleranceForce && iteration<mMaxNumIterations)
            {
                //calculate hessian
                stiffMatrix_jj.SetZeroEntries();stiffMatrix_jk.SetZeroEntries();stiffMatrix_kj.SetZeroEntries();stiffMatrix_kk.SetZeroEntries();
                mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::STIFFNESS, stiffMatrix_jj, stiffMatrix_jk, stiffMatrix_kj, stiffMatrix_kk);
                if (mStructure->GetNumTimeDerivatives()>0)
                {
                    dampingMatrix_jj.SetZeroEntries();dampingMatrix_jk.SetZeroEntries();dampingMatrix_kj.SetZeroEntries();dampingMatrix_kk.SetZeroEntries();
                    mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::DAMPING, dampingMatrix_jj, dampingMatrix_jk, dampingMatrix_kj, dampingMatrix_kk);
                    if (mStructure->GetNumTimeDerivatives()>1)
                    {
                        massMatrix_jj.SetZeroEntries();massMatrix_jk.SetZeroEntries();massMatrix_kj.SetZeroEntries();massMatrix_kk.SetZeroEntries();
                        mStructure->BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType::MASS, massMatrix_jj, massMatrix_jk, massMatrix_kj, massMatrix_kk);
                    }
                }
                //build hessian
                hessian_jj = stiffMatrix_jj;
                if (numEntriesCMat>0)
                {
                    hessian_jj += CmatT * stiffMatrix_kk * Cmat - CmatT * stiffMatrix_kj + stiffMatrix_jk * Cmat;
                }
                if (mStructure->GetNumTimeDerivatives()>0)
                {
                    hessian_jj+=dampingMatrix_jj*(2./(timeStep));
                    if (numEntriesCMat>0)
                    {
                        hessian_jj+=(CmatT * dampingMatrix_kk * Cmat - CmatT * dampingMatrix_kj + dampingMatrix_jk * Cmat )*(2./(timeStep));
                    }
                    if (mStructure->GetNumTimeDerivatives()>1)
                    {
                        hessian_jj+=massMatrix_jj*(4./(timeStep*timeStep));
                        if (numEntriesCMat>0)
                        {
                            hessian_jj+=(CmatT * massMatrix_kk * Cmat - CmatT * massMatrix_kj + massMatrix_jk * Cmat )*(4./(timeStep*timeStep));
                        }
                    }
                }

                //solve for new state
                NuTo::SparseMatrixCSRGeneral<double> hessianModSolver(hessian_jj);
                hessianModSolver.SetOneBasedIndexing();
                mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
                delta_disp_j*=-1;

//std::cout << "delta disp \n" << delta_disp_j.transpose() << std::endl<< std::endl;
                //perform a line search
                double alpha(1);
                double trialNormResidual(0);
                do
                {
                     //calculate trial state
                    trial_disp_j = disp_j + delta_disp_j * alpha;
//std::cout << "trial_disp_j \n" << trial_disp_j.transpose() << std::endl<< std::endl;
                    delta_disp_k = (Cmat*delta_disp_j)*(-1.);
                    trial_disp_k = disp_k + delta_disp_k * alpha;

                    if (mStructure->GetNumTimeDerivatives()>0)
                    {
                        trial_vel_j  = vel_j+delta_disp_j * (2./(timeStep));
                        trial_vel_k  = vel_k+delta_disp_k * (2./(timeStep));
                        if (mStructure->GetNumTimeDerivatives()>1)
                        {
                            trial_acc_j  = acc_j+delta_disp_j * (4./(timeStep*timeStep));
                            trial_acc_k  = acc_k+delta_disp_k * (4./(timeStep*timeStep));
                        }
//std::cout << "line search trial_acc_j\n" << trial_acc_j.transpose() << std::endl<< std::endl;
                    }

                    //apply displacements
                    mStructure->NodeMergeActiveDofValues(0,trial_disp_j);
                    if (mStructure->GetNumTimeDerivatives()>0)
                    {
                        mStructure->NodeMergeDofValues(1,trial_vel_j,trial_vel_k);
                        if (mStructure->GetNumTimeDerivatives()>1)
                            mStructure->NodeMergeDofValues(2,trial_acc_j,trial_acc_k);
                    }
                    mStructure->ElementTotalUpdateTmpStaticData();

                    //calculate internal force, update of history variables=false
                    mStructure->BuildGlobalGradientInternalPotentialSubVectors(intForce_j, intForce_k,false);
                    residual_j = intForce_j;
                    residual_k = intForce_k;
//std::cout << "intForce_j\n" << intForce_j.transpose() << std::endl<< std::endl;
//std::cout << "residual_j\n" << residual_j.transpose() << std::endl<< std::endl;

                    if (mStructure->GetNumTimeDerivatives()>0)
                    {
                        residual_j += dampingMatrix_jj*trial_vel_j+dampingMatrix_jk*trial_vel_k;
                        residual_k += dampingMatrix_kj*trial_vel_j+dampingMatrix_kk*trial_vel_k;
                        if (mStructure->GetNumTimeDerivatives()>1)
                        {
                            residual_j += (massMatrix_jj*trial_acc_j+massMatrix_jk*trial_acc_k);
                            residual_k += (massMatrix_kj*trial_acc_j+massMatrix_kk*trial_acc_k);
//std::cout << "mass j\n" << (massMatrix_jj*trial_acc_j+massMatrix_jk*trial_acc_k).transpose() << std::endl<< std::endl;
//std::cout << "residual_j\n" << residual_j.transpose() << std::endl<< std::endl;
                        }
                    }

                    //add external force (assumed to be independent of the deformation)
                    residual_j -=extForce_j;
                    residual_k -=extForce_k;
//std::cout << "ext force j\n" << extForce_j.transpose() << std::endl<< std::endl;
//std::cout << "residual_j\n" << residual_j.transpose() << std::endl<< std::endl;

                    //std::cout << "residual_j\n " << residual_j << std::endl;
                    //std::cout << "residual_k\n " << residual_k << std::endl;

                    //calculate the initial residual
                    if (numEntriesCMat>0)
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

                    if (mVisualizeResidual)
                    {
                        VisualizeResidual(residual_mod,trial_disp_j, mVisualizeResidualTimeStep);
                        mVisualizeResidualTimeStep++;
                    }
                }
                while(mPerformLineSearch && alpha>mMinLineSearchStep && trialNormResidual>(1.-alpha)*normResidual);

//std::cout << "disp_j after convergence\n" << trial_disp_j.transpose() << std::endl<< std::endl;
                if (alpha>mMinLineSearchStep || !mPerformLineSearch)
                {
                    //improvement is achieved, go to next Newton step
                    disp_j = trial_disp_j;
                    disp_k = trial_disp_k;

                    if (mStructure->GetNumTimeDerivatives()>0)
                    {
                        vel_j  = trial_vel_j;
                        vel_k  = trial_vel_k;
                        if (mStructure->GetNumTimeDerivatives()>1)
                        {
                            acc_j  = trial_acc_j;
                        }
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
                mStructure->ElementTotalUpdateStaticData();

                //store converged step
                lastConverged_disp_j = disp_j;
                lastConverged_disp_k = disp_k;
                mStructure->NodeMergeActiveDofValues(0,disp_j);
                if (mStructure->GetNumTimeDerivatives()>0)
                {
                    lastConverged_vel_j = vel_j;
                    lastConverged_vel_k = vel_k;
                    mStructure->NodeMergeDofValues(1,vel_j,vel_k);
                    if (mStructure->GetNumTimeDerivatives()>1)
                    {
						lastConverged_acc_j = acc_j;
						lastConverged_acc_k = acc_k;
	                    mStructure->NodeMergeDofValues(2,acc_j,acc_k);
                    }
                }
                mStructure->ElementTotalUpdateTmpStaticData();

                //update global time
                mTime+=timeStep;

                mStructure->GetLogger() << "Convergence after " << iteration << " iterations at time " << mTime << "(timestep " << timeStep << ").\n";
                //std::cout << "plot Vector " << plotVector << std::endl;

                //perform Postprocessing
                //mStructure->GetLogger() << " *** PostProcess *** from Crank Nicolson \n";
                PostProcess(residual_j, residual_k);

                //update structure time
                mStructure->SetPrevTime(curTime);

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
                    throw MechanicsException("[NuTo::CrankNicolson::Solve] no convergence with the current maximum number of iterations, either use automatic time stepping, reduce the time step or the minimal line search cut back factor.");
                }
            }
        }
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::CrankNicolson::Solve] performing Newton-Raphson iteration.");
        throw e;
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mStructure->GetLogger()<<"[NuTo::CrankNicolson::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mStructure->GetLogger()<< "[NuTo::CrankNicolson::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return NuTo::Error::SUCCESSFUL;

}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::CrankNicolson::GetTypeId()const
{
    return std::string("CrankNicolson");
}


#ifdef ENABLE_SERIALIZATION
//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
//! @brief ... save the object to a file
void NuTo::CrankNicolson::Restore (const std::string &filename, std::string rType )
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
                throw MechanicsException ( "[CrankNicolson::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[CrankNicolson::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[CrankNicolson::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
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
        throw MechanicsException ( "[CrankNicolson::Restore]Unhandled exception." );
    }
}

//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::CrankNicolson::Save (const std::string &filename, std::string rType )const
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
            throw MechanicsException ( "[CrankNicolson::Save]File type not implemented." );
        }
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[CrankNicolson::Save]File save exception in boost - " ) +std::string ( e.what() ) );
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
        throw MechanicsException ( "[CrankNicolson::Save]Unhandled exception." );
    }
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::CrankNicolson)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
