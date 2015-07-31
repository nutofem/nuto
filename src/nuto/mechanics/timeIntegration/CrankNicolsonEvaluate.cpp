// $Id: CrankNicolsonEvaluate.cpp 575 2011-09-20 18:05:35Z unger3 $

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

#if defined HAVE_PARDISO
    #include <nuto/math/SparseDirectSolverPardiso.h>
#elif defined HAVE_MUMPS
    #include <nuto/math/SparseDirectSolverMUMPS.h>
#else
    std::cout << "Solver not available - can't solve system of equations " << std::endl;
#endif

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/CrankNicolsonEvaluate.h"
#include "nuto/mechanics/timeIntegration/ResultGroupNodeForce.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/mechanics/structures/StructureOutputFullSubvectorsDouble2.h"
#include "nuto/mechanics/structures/StructureOutputSparseSubmatrices4.h"

#include "nuto/mechanics/elements/ElementBase.h"  // delete me
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h" // delete me
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataDamageViscoPlasticity3D.h" // delete me


//! @brief constructor
//! @param mDimension number of nodes
NuTo::CrankNicolsonEvaluate::CrankNicolsonEvaluate (StructureBase* rStructure)  : TimeIntegrationBase (rStructure)
{
    mMinLineSearchStep = 0.01;
    mVisualizeResidualTimeStep = 0;
    mPerformLineSearch = true;
	mToleranceForce = 1e-6;
	mMaxNumIterations = 20;
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::CrankNicolsonEvaluate::Info()const
{
    TimeIntegrationBase::Info();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::CrankNicolsonEvaluate::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::CrankNicolsonEvaluate::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::CrankNicolsonEvaluate::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::CrankNicolsonEvaluate::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::CrankNicolsonEvaluate::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::CrankNicolsonEvaluate::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::CrankNicolsonEvaluate::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of CrankNicolsonEvaluate" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase)
           & BOOST_SERIALIZATION_NVP(mMinLineSearchStep);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of CrankNicolsonEvaluate" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


NuTo::Error::eError NuTo::CrankNicolsonEvaluate::Solve(double rFinalTime)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
    Eigen::initParallel();
    omp_set_num_threads(mStructure->GetNumProcessors());
    Eigen::setNbThreads(mStructure->GetNumProcessors());
#endif
    start=clock();
#endif
    //allocate solver
#if defined HAVE_PARDISO
        NuTo::SparseDirectSolverPardiso mySolver(mStructure->GetNumProcessors(),0);
#elif defined HAVE_MUMPS
        NuTo::SparseDirectSolverMUMPS mySolver;
#else
        std::cout << "Solver not available - can't solve system of equations " << std::endl;
#endif
#ifdef SHOW_TIME
        mySolver.SetShowTime(mStructure->GetShowTime());
#endif
    try
    {

        /*---------------------------------*\
        |***********************************|
        |   Setup time integration scheme   |
        |***********************************|
        \*---------------------------------*/





        /*---------------------------------*\
        |            Basic setup            |
        \*---------------------------------*/

        std::cout << "start Crank Nicolson solver " << std::endl;
#if defined HAVE_PARDISO
        std::cout << "Employing PARDISO-Solver with " << mStructure->GetNumProcessors() << " processors." << std::endl;
#elif defined HAVE_MUMPS
        std::cout << "Employing MUMPS-Solver." << std::endl;
#endif

    	if (mMaxTimeStep==0)
        {
            throw MechanicsException("[NuTo::CrankNicolsonEvaluate::Solve] max time step is set to zero.");
        }

        // renumber dofs and build constraint matrix
        mStructure->NodeBuildGlobalDofs();

        // this is the current starting time of the simulation (think about several load cycles calculated before that)
        double curTime  = mTime;
        double timeStep = mTimeStep;

        // initialize the structure times
        mStructure->SetPrevTime(curTime);
        mStructure->SetTime(curTime);





        /*---------------------------------*\
        |        Allocate Variables         |
        \*---------------------------------*/



        // Matrices
        // --------

        // allocate space for hessian matrix
        SparseMatrixCSRVector2General<double>   hessian_jj(mStructure->GetNumActiveDofs(),
                                                           mStructure->GetNumActiveDofs());


        // allocate space for stiffness matrix
        SparseMatrixCSRVector2General<double>   stiffMatrix_jj(mStructure->GetNumActiveDofs(),
                                                               mStructure->GetNumActiveDofs());
        SparseMatrixCSRVector2General<double>   stiffMatrix_jk(mStructure->GetNumActiveDofs(),
                                                               mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());
        SparseMatrixCSRVector2General<double>   stiffMatrix_kj(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),
                                                               mStructure->GetNumActiveDofs());
        SparseMatrixCSRVector2General<double>   stiffMatrix_kk(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs(),
                                                               mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());


        //calculate individual mass matrix (I was lazy here and have just used general matrices, but symmetric is certainly better here)
        // Declare damping matrix
        SparseMatrixCSRVector2General<double>   dampingMatrix_jj,
                                                dampingMatrix_jk,
                                                dampingMatrix_kj,
                                                dampingMatrix_kk;


        //calculate individual mass matrix (I was lazy here and have just used general matrices, but symmetric is certainly better here)
        // Declare mass matrix
        SparseMatrixCSRVector2General<double>   massMatrix_jj,
                                                massMatrix_jk,
                                                massMatrix_kj,
                                                massMatrix_kk;

        // Constraint matrices
        SparseMatrixCSRVector2General<double>   Cmat    (mStructure->ConstraintGetConstraintMatrixAfterGaussElimination());
        SparseMatrixCSRVector2General<double>   CmatT   (Cmat.Transpose());




        // Vectors
        // -------

        NuTo::FullVector<double,Eigen::Dynamic> disp_j,
                                                disp_k,
                                                vel_j,
                                                vel_k,
                                                acc_j,
                                                acc_k,
                                                delta_disp_j,
                                                delta_disp_k;

        NuTo::FullVector<double,Eigen::Dynamic> trial_disp_j,
                                                trial_disp_k,
                                                trial_vel_j,
                                                trial_vel_k,
                                                trial_acc_j,
                                                trial_acc_k;

        NuTo::FullVector<double,Eigen::Dynamic> lastConverged_disp_j,
                                                lastConverged_disp_k,
                                                lastConverged_vel_j,
                                                lastConverged_vel_k,
                                                lastConverged_acc_j,
                                                lastConverged_acc_k;

        NuTo::FullVector<double,Eigen::Dynamic> extForce_j,
                                                extForce_k,
                                                prevExtForce_j,
                                                prevExtForce_k;

        NuTo::FullVector<double,Eigen::Dynamic> intForce_j(mStructure->GetNumActiveDofs()),
                                                intForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs()),
                                                prevIntForce_j(mStructure->GetNumActiveDofs()),
                                                prevIntForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());

        NuTo::FullVector<double,Eigen::Dynamic> residual_j,
                                                residual_k,
                                                residual_mod,
                                                prevResidual_j,
                                                prevResidual_k;

        // constraints
        FullVector<double,Eigen::Dynamic>       bRHSprev,
                                                bRHS;



        // Skalars
        // -------

        int         numEntriesCMat                  (Cmat.GetNumEntries());





        /*---------------------------------*\
        |     Declare Structure Outputs     |
        \*---------------------------------*/



        // Matrices
        // --------

        NuTo::StructureOutputSparseSubmatrices4     outputStiffness     (stiffMatrix_jj,
                                                                         stiffMatrix_jk,
                                                                         stiffMatrix_kj,
                                                                         stiffMatrix_kk);
        outputStiffness.SetConstant(true);


        NuTo::StructureOutputSparseSubmatrices4     outputDamping       (dampingMatrix_jj,
                                                                         dampingMatrix_jk,
                                                                         dampingMatrix_kj,
                                                                         dampingMatrix_kk);
        outputDamping.SetConstant(true);


        NuTo::StructureOutputSparseSubmatrices4     outputMass          (massMatrix_jj,
                                                                         massMatrix_jk,
                                                                         massMatrix_kj,
                                                                         massMatrix_kk);
        outputMass.SetConstant(true);



        // Vectors
        // -------

        NuTo::StructureOutputFullSubvectorsDouble2  outputIntForce      (intForce_j,
                                                                         intForce_k);

        NuTo::StructureOutputFullSubvectorsDouble2  outputPrevIntForce  (prevIntForce_j,
                                                                         prevIntForce_k);




        /*---------------------------------*\
        |    Declare and fill Output Maps   |
        \*---------------------------------*/

        std::map<NuTo::StructureEnum::eOutput, NuTo::StructureOutputBase*> OutputMapHessianComponents;
        std::map<NuTo::StructureEnum::eOutput, NuTo::StructureOutputBase*> OutputMapInternalForce;
        std::map<NuTo::StructureEnum::eOutput, NuTo::StructureOutputBase*> OutputMapFirstEvaluate;

        OutputMapFirstEvaluate      [StructureEnum::eOutput::INTERNAL_GRADIENT] = &outputPrevIntForce;
        OutputMapInternalForce      [StructureEnum::eOutput::INTERNAL_GRADIENT] = &outputIntForce;

        OutputMapHessianComponents  [StructureEnum::eOutput::STIFFNESS]           = &outputStiffness;
        OutputMapFirstEvaluate      [StructureEnum::eOutput::STIFFNESS]           = &outputStiffness;

        //allocate matrices of higher order (1 = damping, 2 = mass for orders larger than 1)
        if (mStructure->GetNumTimeDerivatives()>0)
        {
            OutputMapHessianComponents[StructureEnum::eOutput::DAMPING] = &outputDamping;
            OutputMapFirstEvaluate    [StructureEnum::eOutput::DAMPING] = &outputDamping;

            if (mStructure->GetNumTimeDerivatives()>1)
            {
                OutputMapHessianComponents[StructureEnum::eOutput::MASS] = &outputMass;
                OutputMapFirstEvaluate    [StructureEnum::eOutput::MASS] = &outputMass;
            }
        }







        /*---------------------------------*\
        |***********************************|
        |         Start Calculation         |
        |***********************************|
        \*---------------------------------*/





        /*---------------------------------*\
        |   Extract last conv. dof values   |
        \*---------------------------------*/

        //store last converged displacements, velocities and accelerations
        mStructure->NodeExtractDofValues(0, lastConverged_disp_j, lastConverged_disp_k);

        if (mStructure->GetNumTimeDerivatives())
        {
            mStructure->NodeExtractDofValues(1, lastConverged_vel_j, lastConverged_vel_k);
            if (mStructure->GetNumTimeDerivatives()>1)
            {
                mStructure->NodeExtractDofValues(2, lastConverged_acc_j, lastConverged_acc_k);
            }
        }





        /*---------------------------------*\
        |         Update Constraints        |
        \*---------------------------------*/

        for(auto itTDC : mMapTimeDependentConstraint)
        {
            mStructure->ConstraintSetRHS(itTDC.first,
                                         itTDC.second->GetTimeDependentFactor(curTime));
        }

        bRHSprev = mStructure->ConstraintGetRHSAfterGaussElimination();
        mStructure->NodeMergeActiveDofValues(0,lastConverged_disp_j); //disp_k is internally calculated from the previous constraint matrix
        mStructure->ElementTotalUpdateTmpStaticData();





        /*---------------------------------*\
        |          Evaluate Outputs         |
        \*---------------------------------*/

        mStructure->Evaluate(OutputMapFirstEvaluate);

        // remove all constant outputs from output maps to avoid unneccessary calculations
        for(auto itOutputMap : OutputMapFirstEvaluate)
        {
            if (itOutputMap.second->GetConstant())
            {
                OutputMapHessianComponents.erase(itOutputMap.first);
            }
        }





        /*---------------------------------*\
        |    Calculate previous residual    |
        \*---------------------------------*/

        prevResidual_j = prevIntForce_j;
        prevResidual_k = prevIntForce_k;

        if (mStructure->GetNumTimeDerivatives()>0)
        {
            if (outputDamping.GetConstant())
            {
                prevResidual_j += dampingMatrix_jj*lastConverged_vel_j+dampingMatrix_jk*lastConverged_vel_k;
                prevResidual_k += dampingMatrix_kj*lastConverged_vel_j+dampingMatrix_kk*lastConverged_vel_k;
            }
            if (mStructure->GetNumTimeDerivatives()>1)
            {
                if(outputMass.GetConstant())
                {
                    prevResidual_j += (massMatrix_jj*lastConverged_acc_j+massMatrix_jk*lastConverged_acc_k);
                    prevResidual_k += (massMatrix_kj*lastConverged_acc_j+massMatrix_kk*lastConverged_acc_k);
                }
            }
        }

        //add external force
        CalculateExternalLoad(*mStructure, curTime, prevExtForce_j, prevExtForce_k);
        prevResidual_j -= prevExtForce_j;
        prevResidual_k -= prevExtForce_k;

        residual_mod=prevResidual_j - CmatT*prevResidual_k;

        //@COMMENT check, if this could be avoided, only required if check for initial equilibrium required
        //std::cout << "residual in initial configuration " << residual_mod.Norm() << std::endl;

        if (residual_mod.Norm()>mToleranceForce && mCheckEquilibriumOnStart)
        {
            throw MechanicsException("[NuTo::CrankNicolsonEvaluate::Solve] Initial configuration is not in (dynamic) equilibrium.");
        }
        PostProcess(prevResidual_j, prevResidual_k);







        /*---------------------------------*\
        |***********************************|
        |          Start Main Loop          |
        |***********************************|
        \*---------------------------------*/

        while (curTime < rFinalTime)
        {





            /*---------------------------------*\
            |    Prev. Constraints and Loads    |
            \*---------------------------------*/

            for(auto itTDC : mMapTimeDependentConstraint)
            {
                mStructure->ConstraintSetRHS(itTDC.first,
                                             itTDC.second->GetTimeDependentFactor(curTime));
            }
            bRHSprev = mStructure->ConstraintGetRHSAfterGaussElimination();
            mStructure->NodeMergeActiveDofValues(0,lastConverged_disp_j); //disp_k is internally calculated from the previous constraint matrix
            mStructure->ElementTotalUpdateTmpStaticData();


            //add previous external force (because delta is required)
            CalculateExternalLoad(*mStructure, curTime, prevExtForce_j, prevExtForce_k);





            /*---------------------------------*\
            |             Adjust Time           |
            \*---------------------------------*/

            //increase time step
            curTime += timeStep;

            if (timeStep<mMinTimeStep)
                throw MechanicsException("[NuTo::CrankNicolsonEvaluate::Solve] time step is smaller than minimum - no convergence is obtained.");

            //check if curTime is close to final time
            if (curTime>rFinalTime || (rFinalTime-curTime)<0.2*timeStep)
            {
                timeStep+=rFinalTime-curTime;
                curTime = rFinalTime;
            }

            // set new structure time at the end of the new time increment
            mStructure->SetTime(curTime);





            /*---------------------------------*\
            |    Update Constraints and Loads   |
            \*---------------------------------*/

            CalculateExternalLoad(*mStructure, curTime, extForce_j, extForce_k);
            residual_j = prevExtForce_j-extForce_j;
            residual_k = prevExtForce_k-extForce_k;

            for(auto itTDC : mMapTimeDependentConstraint)
            {
                mStructure->ConstraintSetRHS(itTDC.first,
                                             itTDC.second->GetTimeDependentFactor(curTime));
            }
            bRHS = mStructure->ConstraintGetRHSAfterGaussElimination();
            FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deltaBRHS(bRHS-bRHSprev);

             NuTo::FullMatrix<double, Eigen::Dynamic,Eigen::Dynamic> test(Cmat);






            /*---------------------------------*\
            | Evaluate Stiffness, Damping, Mass |
            \*---------------------------------*/

            mStructure->Evaluate(OutputMapHessianComponents);





            /*---------------------------------*\
            |         Calculate Residual        |
            \*---------------------------------*/

            residual_j += stiffMatrix_jk * deltaBRHS;
            if (mStructure->GetNumTimeDerivatives()>0)
            {
                if(outputDamping.GetConstant())
                {
                    residual_j -= (dampingMatrix_jj*lastConverged_vel_j + (dampingMatrix_jk * (deltaBRHS *(-1./timeStep)+lastConverged_vel_k)))*2.;
                }
                if (mStructure->GetNumTimeDerivatives()>1)
                {
                    if(outputMass.GetConstant())
                    {
                    residual_j -= massMatrix_jj * (lastConverged_vel_j*(4./timeStep)+lastConverged_acc_j*2.)+
                                  massMatrix_jk * (deltaBRHS *(-4./(timeStep*timeStep))+lastConverged_vel_k*(4./timeStep)+lastConverged_acc_k*2.);
                    }
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





            /*---------------------------------*\
            |   Build hessian for trial state   |
            \*---------------------------------*/

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





            /*---------------------------------*\
            |         Solve trial state         |
            \*---------------------------------*/

            NuTo::SparseMatrixCSRGeneral<double> hessianModSolver(hessian_jj);
            hessianModSolver.SetOneBasedIndexing();
            mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
            delta_disp_j*=-1;





            /*---------------------------------*\
            |  Calculate new disp, vel and acc  |
            \*---------------------------------*/

            //calculate trial state
            disp_j = lastConverged_disp_j + delta_disp_j;
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
            }

            //apply displacements, velocities and accelerations
            mStructure->NodeMergeActiveDofValues(0,disp_j);
            if (mStructure->GetNumTimeDerivatives()>0)
            {
                mStructure->NodeMergeActiveDofValues(1,vel_j);
                if (mStructure->GetNumTimeDerivatives()>1)
                {
                    mStructure->NodeMergeActiveDofValues(2,acc_j);
                }
            }

            mStructure->ElementTotalUpdateTmpStaticData();





            /*---------------------------------*\
            |         Calculate Residual        |
            \*---------------------------------*/

            mStructure->Evaluate(OutputMapInternalForce);
            residual_j = intForce_j;
            residual_k = intForce_k;

            if (mStructure->GetNumTimeDerivatives()>0)
            {
                if (outputDamping.GetConstant())
                {
                    residual_j += dampingMatrix_jj*vel_j+dampingMatrix_jk*vel_k;
                    residual_k += dampingMatrix_kj*vel_j+dampingMatrix_kk*vel_k;
                }
                if (mStructure->GetNumTimeDerivatives()>1)
                {
                    if(outputMass.GetConstant())
                    {
                        residual_j += (massMatrix_jj*acc_j+massMatrix_jk*acc_k);
                        residual_k += (massMatrix_kj*acc_j+massMatrix_kk*acc_k);
                    }
                }
            }

            //add external force
            residual_j -= extForce_j;
            residual_k -= extForce_k;

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
            std::cout << "norm of residual after prediction state:" << normResidual << std::endl<< std::endl;







            /*---------------------------------*\
            |***********************************|
            |          Start Iteration          |
            |***********************************|
            \*---------------------------------*/

            int iteration(0);
            while(normResidual>mToleranceForce && iteration<mMaxNumIterations)
            {





                /*---------------------------------*\
                |           Build hessian           |
                \*---------------------------------*/

                mStructure->Evaluate(OutputMapHessianComponents);

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





                /*---------------------------------*\
                |               Solve               |
                \*---------------------------------*/

                NuTo::SparseMatrixCSRGeneral<double> hessianModSolver(hessian_jj);
                hessianModSolver.SetOneBasedIndexing();
                mySolver.Solve(hessianModSolver, residual_mod, delta_disp_j);
                delta_disp_j*=-1;







                /*---------------------------------*\
                |***********************************|
                |        Perform Line Search        |
                |***********************************|
                \*---------------------------------*/

                // Because of the do-loop, the following code is executed (on purpose) even if linesearch is disabled
                double alpha(1.);
                double trialNormResidual(0.);
                do
                {





                    /*---------------------------------*\
                    | Calculate trial disp, vel and acc |
                    \*---------------------------------*/

                     //calculate trial state
                    trial_disp_j = disp_j + delta_disp_j * alpha;;
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
                    }

                    //apply displacements
                    mStructure->NodeMergeActiveDofValues(0,trial_disp_j);
                    if (mStructure->GetNumTimeDerivatives()>0)
                    {
                        mStructure->NodeMergeDofValues(1,trial_vel_j,trial_vel_k);
                        if (mStructure->GetNumTimeDerivatives()>1)
                        {
                            mStructure->NodeMergeDofValues(2,trial_acc_j,trial_acc_k);
                        }
                    }

                    mStructure->ElementTotalUpdateTmpStaticData();





                    /*---------------------------------*\
                    |         Calculate Residual        |
                    \*---------------------------------*/

                    mStructure->Evaluate(OutputMapInternalForce);

                    residual_j = intForce_j;
                    residual_k = intForce_k;

                    if (mStructure->GetNumTimeDerivatives()>0)
                    {
                        if(outputDamping.GetConstant())
                        {
                            //residual_j += dampingMatrix_jj*trial_vel_j+dampingMatrix_jk*trial_vel_k;
                            //residual_k += dampingMatrix_kj*trial_vel_j+dampingMatrix_kk*trial_vel_k;
                        }
                        if (mStructure->GetNumTimeDerivatives()>1)
                        {
                            if(outputMass.GetConstant())
                            {
                                residual_j += (massMatrix_jj*trial_acc_j+massMatrix_jk*trial_acc_k);
                                residual_k += (massMatrix_kj*trial_acc_j+massMatrix_kk*trial_acc_k);
                            }
                        }
                    }

                    //add external force (assumed to be independent of the deformation)
                    residual_j -=extForce_j;
                    residual_k -=extForce_k;

                    //calculate the initial residual
                    if (numEntriesCMat>0)
                    {
                        residual_mod = residual_j - CmatT*residual_k;
                    }
                    else
                    {
                        residual_mod = residual_j;
                    }

                    trialNormResidual=residual_mod.Norm();

                    if(mPerformLineSearch)
                    {
                        mStructure->GetLogger() << "  linesearch alpha " << alpha <<" previous residual " << normResidual << " trial residual " << trialNormResidual <<  "\n";
                    }
                    alpha*=0.5;

                    if (mVisualizeResidual)
                    {
                        VisualizeResidual(residual_mod,trial_disp_j, mVisualizeResidualTimeStep);
                        mVisualizeResidualTimeStep++;
                    }
                }
                while(mPerformLineSearch && alpha>mMinLineSearchStep && trialNormResidual>(1.-alpha)*normResidual);
                // END of linesearch





                /*---------------------------------*\
                |      Update disp, vel and acc     |
                \*---------------------------------*/

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

            } // END of iteration





            /*---------------------------------*\
            |  Update last conv disp, vel & acc |
            \*---------------------------------*/

            if (normResidual<=mToleranceForce)
            {
                //converged solution
                mStructure->ElementTotalUpdateStaticData();

                //store converged step
                lastConverged_disp_j = disp_j;
                lastConverged_disp_k = disp_k;
                mStructure->NodeMergeActiveDofValues(0,disp_j);         // not necessary???
                if (mStructure->GetNumTimeDerivatives()>0)
                {
                    lastConverged_vel_j = vel_j;
                    lastConverged_vel_k = vel_k;
                    mStructure->NodeMergeDofValues(1,vel_j,vel_k);      // not necessary???
                    if (mStructure->GetNumTimeDerivatives()>1)
                    {
						lastConverged_acc_j = acc_j;
						lastConverged_acc_k = acc_k;
                        mStructure->NodeMergeDofValues(2,acc_j,acc_k);  // not necessary???
                    }
                }
                mStructure->ElementTotalUpdateTmpStaticData();





                /*---------------------------------*\
                |             Adjust Time           |
                \*---------------------------------*/

                mTime+=timeStep;





                /*---------------------------------*\
                |         Post Processing           |
                \*---------------------------------*/

                mStructure->GetLogger() << "Convergence after " << iteration << " iterations at time " << mTime << "(timestep " << timeStep << ").\n";

                //perform Postprocessing
                //mStructure->GetLogger() << " *** PostProcess *** from Crank Nicolson \n";
                PostProcess(residual_j, residual_k);

                //update structure time
                mStructure->SetPrevTime(curTime);





                /*---------------------------------*\
                |      Automatic Timestepping       |
                \*---------------------------------*/

                if (mAutomaticTimeStepping && iteration<0.25*mMaxNumIterations)
                {
                    timeStep*=1.5;
                    if (timeStep>mMaxTimeStep)
                        timeStep = mMaxTimeStep;
                }
            }






            /*---------------------------------*\
            |***********************************|
            |   Handle not converged solutions  |
            |***********************************|
            \*---------------------------------*/

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
                    throw MechanicsException("[NuTo::CrankNicolsonEvaluate::Solve] no convergence with the current maximum number of iterations, either use automatic time stepping, reduce the time step or the minimal line search cut back factor.");
                }
            }
        } // END of main loop
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::CrankNicolsonEvaluate::Solve] performing Newton-Raphson iteration.");
        throw e;
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mStructure->GetLogger()<<"[NuTo::CrankNicolsonEvaluate::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mStructure->GetLogger()<< "[NuTo::CrankNicolsonEvaluate::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return NuTo::Error::SUCCESSFUL;

}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::CrankNicolsonEvaluate::GetTypeId()const
{
    return std::string("CrankNicolsonEvaluate");
}


#ifdef ENABLE_SERIALIZATION
//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
//! @brief ... save the object to a file
void NuTo::CrankNicolsonEvaluate::Restore (const std::string &filename, std::string rType )
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
                throw MechanicsException ( "[CrankNicolsonEvaluate::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[CrankNicolsonEvaluate::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[CrankNicolsonEvaluate::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
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
        throw MechanicsException ( "[CrankNicolsonEvaluate::Restore]Unhandled exception." );
    }
}

//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::CrankNicolsonEvaluate::Save (const std::string &filename, std::string rType )const
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
            throw MechanicsException ( "[CrankNicolsonEvaluate::Save]File type not implemented." );
        }
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[CrankNicolsonEvaluate::Save]File save exception in boost - " ) +std::string ( e.what() ) );
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
        throw MechanicsException ( "[CrankNicolsonEvaluate::Save]Unhandled exception." );
    }
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::CrankNicolsonEvaluate)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
