#include "mechanics/timeIntegration/TimeIntegrationBase.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "base/ErrorEnum.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/MechanicsException.h"
#include "base/Timer.h"


#include "mechanics/structures/StructureBase.h"
#include "mechanics/structures/unstructured/StructureFETI.h"
#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputDummy.h"
#include <boost/mpi.hpp>

#include "base/CallbackInterface.h"
#include "math/SparseMatrixCSRGeneral.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "mechanics/timeIntegration/NewmarkFeti.h"

#include <cmath>


Eigen::MatrixXd NuTo::NewmarkFeti::GatherInterfaceRigidBodyModes(const Eigen::MatrixXd& interfaceRigidBodyModes, const int numRigidBodyModesGlobal)
{

    std::vector<int> recvCount;
    std::vector<int> displ;
    MpiGatherRecvCountAndDispls(recvCount, displ, interfaceRigidBodyModes.size());

    const int numInterfaceEqs               = interfaceRigidBodyModes.rows();
    MatrixXd interfaceRigidBodyModesGlobal   = MatrixXd::Zero(numInterfaceEqs,numRigidBodyModesGlobal);

    MPI_Allgatherv(interfaceRigidBodyModes.data(),
                   interfaceRigidBodyModes.size(),
                   MPI_DOUBLE,
                   interfaceRigidBodyModesGlobal.data(),
                   recvCount.data(),
                   displ.data(),
                   MPI_DOUBLE,
                   MPI_COMM_WORLD);

    return interfaceRigidBodyModesGlobal;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Eigen::VectorXd NuTo::NewmarkFeti::GatherRigidBodyForceVector(Eigen::VectorXd &rigidBodyForceVectorLocal, const int numRigidBodyModesGlobal)
{

    const int numRigidBodyModesLocal        = rigidBodyForceVectorLocal.rows();
    std::vector<int> recvCount;
    std::vector<int> displ;
    MpiGatherRecvCountAndDispls(recvCount, displ, numRigidBodyModesLocal);

    VectorXd rigidBodyForceVectorGlobal      =VectorXd::Zero(numRigidBodyModesGlobal);
    MPI_Allgatherv(rigidBodyForceVectorLocal.data(),
                   rigidBodyForceVectorLocal.size(),
                   MPI_DOUBLE,
                   rigidBodyForceVectorGlobal.data(),
                   recvCount.data(),
                   displ.data(),
                   MPI_DOUBLE,
                   MPI_COMM_WORLD);

    return rigidBodyForceVectorGlobal;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void NuTo::NewmarkFeti::MpiGatherRecvCountAndDispls(std::vector<int> &recvCount, std::vector<int> &displs, const int numValues)
{
    boost::mpi::communicator world;
    const int numProcesses = world.size();
    // recvCount:
    // Contais the number of elements that are received from each process.
    recvCount.clear();
    recvCount.resize(numProcesses, 0);

    boost::mpi::all_gather<int>(world,numValues,recvCount);

    // displs:
    // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i.
    displs.clear();
    displs.resize(numProcesses, 0);
    for (int i = 1; i < numProcesses; ++i)
        displs[i] = displs[i-1] + recvCount[i-1];

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int NuTo::NewmarkFeti::BiCgStab(const MatrixXd &projection, VectorXd &x, const VectorXd &rhs)
{

    StructureFETI* structure = static_cast<StructureFETI*>(mStructure);
    const SparseMatrix& B      = structure->GetConnectivityMatrix();
    const SparseMatrix Btrans = B.transpose();


    // solve the linear system Ax = b
    VectorXd Ax = B * mSolver.solve(Btrans*x);
    MPI_Allreduce(MPI_IN_PLACE, Ax.data(), Ax.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    VectorXd& Ay = Ax;
    VectorXd& Az = Ax;

    const int n = rhs.rows();
    VectorXd r  = rhs - Ax;
//    VectorXd r0 = r;

    VectorXd rProj   = projection * r;
    VectorXd rProj0  = rProj;

    double rProj0_sqnorm = rProj0.squaredNorm();
    double rhs_sqnorm = rhs.squaredNorm();

    double rho    = 1.;
    double alpha  = 1.;
    double w      = 1.;

    VectorXd v = VectorXd::Zero(n);
    VectorXd p = VectorXd::Zero(n);
    VectorXd y(n);
    VectorXd z(n);
    VectorXd s(n);
    VectorXd t(n);

//    const double tol    = mCpgTolerance;
    //        const double tol2   = tol*tol*rhs_sqnorm;
    const double tol2   = 1.0e-8*rhs_sqnorm; // the phase-field model is worse conditioned it seems. Lower tolerance required?

    //const double eps    = std::numeric_limits<double>::epsilon();
    //const double eps2   = eps*eps;

    int i = 0;
    int restarts = 0;



    const int maxIters = mCpgMaxIterations;


    while (rProj.squaredNorm() > tol2 and i<maxIters )
    {
        double rho_old = rho;

        rho = rProj0.dot(rProj);



        if (std::fabs(rho) < 1.0e-20)
        {

            if (structure->mRank == 0)
            {
                std::cout << "The new residual vector became too orthogonal to the arbitrarily chosen direction r0.\n"
                             "Let's restart with a new r0!" << std::endl;
            }

            // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
            // Let's restart with a new r0:

            Ax.noalias() = B * mSolver.solve(Btrans*x);
            MPI_Allreduce(MPI_IN_PLACE, Ax.data(), Ax.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


            rProj  = projection*(rhs - Ax);
            rProj0 = rProj;
            rho = rProj0_sqnorm = rProj.squaredNorm();
            if(restarts++ == 0)
                i = 0;
        }


        double beta = (rho/rho_old) * (alpha / w);

        p = rProj + beta * (p - w * v);


        //      y = precond.solve(p);
        y = p;

        Ay.noalias() = B * mSolver.solve(Btrans*y);
        MPI_Allreduce(MPI_IN_PLACE, Ay.data(), Ay.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        v.noalias() = projection * Ay;

        alpha = rho / rProj0.dot(v);

        s = rProj - alpha * v;

        //      z = precond.solve(s);
        z = s;

        Az.noalias() = B * mSolver.solve(Btrans*z);
        MPI_Allreduce(MPI_IN_PLACE, Az.data(), Az.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        t.noalias() = projection * Az;

        double tmp = t.dot(t);

        if(tmp>0.0)
            w = t.dot(s) / tmp;
        else
            w = 0.0;

        x += alpha * y + w * z;
        rProj = z - w * t;

        if (structure->mRank == 0)
            std::cout << "BiCGStab rel. error = " << rProj.squaredNorm()/rhs_sqnorm << "\t at iteration = " << i << "/" << maxIters << std::endl;

        ++i;



    }
    //    tol_error = sqrt(r.squaredNorm()/rhs_sqnorm);
    return i;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int NuTo::NewmarkFeti::CPG(const Eigen::MatrixXd &projection, Eigen::VectorXd &x, const Eigen::VectorXd &rhs)
{

    boost::mpi::communicator world;
    VectorXd        Ap;

    StructureFETI* structure = static_cast<StructureFETI*>(mStructure);

    const SparseMatrix& B       = structure->GetConnectivityMatrix();
    const SparseMatrix Btrans   = B.transpose();

    world.barrier();

    const int numActiveDofs     = mStructure->GetNumActiveDofs(Node::eDof::DISPLACEMENTS);
    const int numRigidBodyModes = structure->GetNumRigidBodyModes();

    VectorXd tmp;
    tmp.setZero(numActiveDofs + numRigidBodyModes);
    tmp.head(numActiveDofs) = mSolver.solve( (Btrans*x).head(numActiveDofs));
    VectorXd Ax = B * tmp;


//    std::cout << "Rank: " << structure->mRank << " tmp " << tmp << std::endl;
    world.barrier();


    boost::mpi::all_reduce(world,boost::mpi::inplace(Ax.data()), Ax.size(),std::plus<double>());

    world.barrier();
    // initial residual
    VectorXd r = rhs - Ax;
    VectorXd z;

    // initial projected search direction
    VectorXd w = projection * r;
    VectorXd p = w;

    // precondition
//    p =  projection * mLocalPreconditioner * w;
//    boost::mpi::all_reduce(world,boost::mpi::inplace(p.data()), p.size(),std::plus<double>());
    p = w;

    const double rhs_sqnorm = rhs.squaredNorm();
    const double threshold = mCpgTolerance * mCpgTolerance * rhs_sqnorm;

    double absNew = w.dot(p);
    int iteration = 0;
    while(iteration < mCpgMaxIterations)
    {
        // at every iteration i the r has to be recomputd which is quite expensive
        world.barrier();

        tmp.head(numActiveDofs) = mSolver.solve( (Btrans*p).head(numActiveDofs));
        Ap.noalias() = B * tmp;

        world.barrier();
        boost::mpi::all_reduce(world,boost::mpi::inplace(Ap.data()), Ap.size(),std::plus<double>());

        // step size
        double alpha = absNew / p.dot(Ap);

        // update solution
        x += alpha * p;

        // update projected residual
        w -= alpha * projection * Ap;

        // precondition
//        z =  projection * mLocalPreconditioner * w;
//        boost::mpi::all_reduce(world,boost::mpi::inplace(z.data()), z.size(),std::plus<double>());
        z = w;

        if (w.squaredNorm() < threshold)
            break;

        double absOld = absNew;
        absNew = w.dot(z);                        // update the absolute value of r
        double beta = absNew / absOld;            // calculate the Gram-Schmidt value used to create the new search direction
        p = z + beta * p;                         // update search direction

        structure->GetLogger() << "CPG rel. error = "   << w.squaredNorm()/rhs_sqnorm
                               << "\t at iteration = "  << iteration << "/" << mCpgMaxIterations <<"\n";


        ++iteration;


    }

    structure->GetLogger() << "\nConverged! \n"
                           << "CPG rel. error = "   << w.squaredNorm()/rhs_sqnorm
                           << "\t at iteration = "  << iteration << "/" << mCpgMaxIterations <<"\n";

    return iteration;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


NuTo::StructureOutputBlockVector NuTo::NewmarkFeti::FetiSolve(const VectorXd& residual_mod, const std::set<Node::eDof> &activeDofSet, VectorXd &deltaLambda)
{


    boost::mpi::communicator world;

    StructureFETI* structure = static_cast<StructureFETI*>(mStructure);



    const SparseMatrix& B                    = structure->GetConnectivityMatrix();
    const SparseMatrix  Btrans               = B.transpose();

    const int numRigidBodyModesLocal        = structure->GetNumRigidBodyModes();
    const int numRigidBodyModesGlobal       = structure->GetNumRigidBodyModesTotal();

    const int numTotalActiveDofs            = mStructure->GetNumTotalActiveDofs();
    const int numTotalDofs                  = mStructure->GetNumTotalDofs();

    const MatrixXd& rigidBodyModes              = structure->GetRigidBodyModes();

    world.barrier();

    // G.size() = (number of Lagrange multipliers) x (total number of rigid body modes)
    const MatrixXd& G          =   structure->GetG();

    // GtransGinv.size() = (total number of rigid body modes) x (total number of rigid body modes)
    const MatrixXd GtransGinv =   (G.transpose() * G).inverse();


    // R_s^T f
    VectorXd rigidBodyForceVectorLocal = rigidBodyModes.transpose() * residual_mod;

    //
    //     | R_1^T f        |
    // e = |     :          |
    //     | R_{N_s}^T f    |
    //
    VectorXd rigidBodyForceVectorGlobal = GatherRigidBodyForceVector(rigidBodyForceVectorLocal, numRigidBodyModesGlobal);



    // initial guess for lambda
    deltaLambda = G * GtransGinv * rigidBodyForceVectorGlobal;

    world.barrier();

    //*****************************************
    //
    //       | K_11^-1 * r |
    // d = B |  0          |
    //
    //*****************************************
    VectorXd pseudoInvVec;
    pseudoInvVec.setZero(numTotalDofs);
    pseudoInvVec.head(numTotalActiveDofs)   = mSolver.solve(residual_mod.head(numTotalActiveDofs));
    VectorXd displacementGap            = B * pseudoInvVec;
    boost::mpi::all_reduce(world,boost::mpi::inplace(displacementGap.data()), displacementGap.size(),std::plus<double>());


    MPI_Barrier(MPI_COMM_WORLD);

    VectorXd zeroVec = G.transpose() * deltaLambda - rigidBodyForceVectorGlobal;

    mStructure->GetLogger() << "Rank: \t" << structure->mRank << "\n Gtrans * lambda - e = 0? \n" << zeroVec.norm() << "\n \n";
    mStructure->GetLogger() << "Rank: \t" << structure->mRank << "\n residual_mod.norm() \n" << residual_mod.norm() << "\n \n";

    MPI_Barrier(MPI_COMM_WORLD);

    structure->GetLogger() << "\n*************************************\n";
    structure->GetLogger() << "       Start Interface Problem           ";
    structure->GetLogger() << "\n*************************************\n\n";

    MPI_Barrier(MPI_COMM_WORLD);

    int iteration = CPG(structure->GetProjectionMatrix(),deltaLambda,displacementGap);
    //int iteration = BiCgStab(projection,x,rhs);

    if (iteration >= mCpgMaxIterations)
        throw MechanicsException(__PRETTY_FUNCTION__,"Maximum number of iterations exceeded.");

    structure->GetLogger() << "\n*************************************\n";
    structure->GetLogger() << "       End Interface Problem             ";
    structure->GetLogger() << "\n*************************************\n\n";

    MPI_Barrier(MPI_COMM_WORLD);
    zeroVec = G.transpose() * deltaLambda - rigidBodyForceVectorGlobal;


    structure->GetLogger() << "Gtrans * lambda - e = 0? \n" << zeroVec.norm() << "\n \n";

    structure->GetLogger() << "Rtrans ( f - Btrans*lambda ) = 0? \n" << rigidBodyModes.transpose() * ( residual_mod - B.transpose() *deltaLambda ) << "\n \n";


    MPI_Barrier(MPI_COMM_WORLD);

    //*****************************************
    //
    //         | K_11^-1 * Btrans*lambda    |
    // tmp = B |  0                         |
    //
    //*****************************************
    pseudoInvVec.setZero(mStructure->GetNumTotalDofs());
    pseudoInvVec.head(mStructure->GetNumTotalActiveDofs()) = mSolver.solve((Btrans*deltaLambda).head(mStructure->GetNumTotalActiveDofs()));

    VectorXd tmp = B * pseudoInvVec;
    boost::mpi::all_reduce(world,boost::mpi::inplace(tmp.data()), tmp.size(),std::plus<double>());

    VectorXd alphaGlobal =VectorXd::Zero(numRigidBodyModesGlobal);
    alphaGlobal = GtransGinv * G.transpose() * (displacementGap - tmp);


    MPI_Barrier(MPI_COMM_WORLD);
    VectorXd zeroVecTmp = G * alphaGlobal - (displacementGap - tmp);
    structure->GetLogger() << "Rank: \t" << structure->mRank << "\n G*alpha - (d- F*lambda) = 0? \n" << zeroVecTmp.norm() << "\n \n";

    VectorXd zeroVecTmp2 = structure->GetProjectionMatrix().transpose() * (displacementGap - tmp);
    structure->GetLogger() << "Rank: \t" << structure->mRank << "\n Ptrans * (d- F*lambda) = 0? \n" << zeroVecTmp2.norm() << "\n \n";

    MPI_Barrier(MPI_COMM_WORLD);



    std::vector<int> recvCount;
    std::vector<int> displs;
    MpiGatherRecvCountAndDispls(recvCount, displs, numRigidBodyModesLocal);
    VectorXd alphaLocal = alphaGlobal.segment(displs[structure->mRank], numRigidBodyModesLocal);


    VectorXd delta_dof_active;
    delta_dof_active.setZero(mStructure->GetNumTotalDofs());
    delta_dof_active.head(mStructure->GetNumTotalActiveDofs())    = mSolver.solve( (residual_mod - Btrans * deltaLambda).head(structure->GetNumTotalActiveDofs()));

    delta_dof_active   = delta_dof_active - rigidBodyModes * alphaLocal;



    StructureOutputBlockVector  delta_dof_dt0(structure->GetDofStatus());
    int offset = 0;
    for (const auto& dof : activeDofSet)
    {
        int numActiveDofs       = structure->GetNumActiveDofs(dof);
        delta_dof_dt0.J[dof]    = delta_dof_active.segment(offset, structure->GetNumActiveDofs(dof));
        /// \todo This will not work for multi field problems, please implement it more generally
        delta_dof_dt0.K[dof]    = delta_dof_active.tail(numRigidBodyModesLocal);
        offset += numActiveDofs;
    }




    return delta_dof_dt0;


}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


NuTo::eError NuTo::NewmarkFeti::Solve(double rTimeDelta)
{

    try
    {
        boost::mpi::communicator world;
        StructureFETI* structure = static_cast<StructureFETI*>(mStructure);
        mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

        mStructure->GetLogger() << "********************************************" << "\n";
        mStructure->GetLogger() << "**        Calculate rigid body modes      **" << "\n";
        mStructure->GetLogger() << "********************************************" << "\n\n";

        structure->CalculateRigidBodyModesTotalFETI();

        const int numRigidBodyModes = structure->GetNumRigidBodyModes();

        mStructure->GetLogger() << "************************************************" << "\n";
        mStructure->GetLogger() << "**        Assemble connectivity matrix B      **" << "\n";
        mStructure->GetLogger() << "************************************************" << "\n\n";

        structure->AssembleConnectivityMatrix();
        const SparseMatrix& B       = structure->GetConnectivityMatrix();
        const SparseMatrix Btrans   = B.transpose();

        mStructure->GetLogger() << "******************************************************" << "\n";
        mStructure->GetLogger() << "**        Calculate interface rigid body modes      **" << "\n";
        mStructure->GetLogger() << "******************************************************" << "\n\n";

        structure->CalculateInterfaceRigidBodyModes();

        mStructure->GetLogger() << "********************************************" << "\n";
        mStructure->GetLogger() << "**        Calculate projection            **" << "\n";
        mStructure->GetLogger() << "********************************************" << "\n\n";

        structure->CalculateG();
        structure->CalculateProjectionMatrix();

//        structure->CheckProjectionMatrix();
//        structure->CheckProjectionOfCoarseGrid();

        const int numActiveDofs     = mStructure->GetNumTotalActiveDofs();
        const int numTotalDofs      = mStructure->GetNumTotalDofs();


        VectorXd deltaLambda            = VectorXd::Zero(B.rows());
        VectorXd lambda                 = VectorXd::Zero(B.rows());
        VectorXd lastConverged_lambda   = VectorXd::Zero(B.rows());

        double curTime  = mTime;
        double timeStep = mTimeStep;
        mStructure->SetPrevTime(curTime);
        mStructure->SetTime(curTime);


        CalculateStaticAndTimeDependentExternalLoad();

        // sets hasInteractingConstraints to true to also assemble K_KJ and K_KK
        // only necessary for CheckRigidBodyModes
        mStructure->DofStatusSetHasInteractingConstraints(true);
        const DofStatus& dofStatus = mStructure->GetDofStatus();

        assert(dofStatus.HasInteractingConstraints()==true);

        if(mStepActiveDofs.empty())
        {
            mStepActiveDofs.push_back(mStructure->DofTypesGetActive());
        }
        else
        {
            for(unsigned int i=0; i<mStepActiveDofs.size();++i)
            {
                if(mStepActiveDofs[i].empty())
                {
                    throw MechanicsException(__PRETTY_FUNCTION__, "Calculation step " +std::to_string(i)+ " has no active DOFs.");
                }
            }
        }


        //********************************************
        //        Allocate Variables
        //********************************************

        StructureOutputBlockMatrix  hessian0(dofStatus, true);

        StructureOutputBlockVector  delta_dof_dt0(dofStatus, true);

        StructureOutputBlockVector  dof_dt0(dofStatus, true); // e.g. disp

        StructureOutputBlockVector  lastConverged_dof_dt0(dofStatus, true); // e.g. disp
        StructureOutputBlockVector  lastConverged_dof_dt1(dofStatus, true); // e.g. velocity
        StructureOutputBlockVector  lastConverged_dof_dt2(dofStatus, true); // e.g. accelerations

        StructureOutputBlockVector  extForce(dofStatus, true);
        StructureOutputBlockVector  intForce(dofStatus, true);
        StructureOutputBlockVector  residual(dofStatus, true);

        StructureOutputBlockVector  prevResidual(dofStatus, true);


        // for constraints
        // ---------------
//        const auto& cmat = mStructure->GetConstraintMatrix();

        BlockScalar normResidual(dofStatus);

        //********************************************
        //        Declare and fill Output maps
        //********************************************

        // Declare output maps

        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradient;
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradientAndHessian0;
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalHessian0;

        evalInternalGradient                [eStructureOutput::INTERNAL_GRADIENT]   = &intForce;

        evalInternalGradientAndHessian0     [eStructureOutput::INTERNAL_GRADIENT]   = &intForce;
        evalInternalGradientAndHessian0     [eStructureOutput::HESSIAN0]            = &hessian0;

        evalHessian0                        [eStructureOutput::HESSIAN0]            = &hessian0;

        //********************************************
        //        Declare and fill Input map
        //********************************************

        ConstitutiveInputMap inputMap;
        inputMap[Constitutive::eInput::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(
                    eCalculateStaticData::EULER_BACKWARD);

        ExtractDofValues(lastConverged_dof_dt0, lastConverged_dof_dt1, lastConverged_dof_dt2);

        UpdateAndGetAndMergeConstraintRHS(curTime, lastConverged_dof_dt0);

        PostProcess(residual);





        mStructure->GetLogger() << "********************************************" << "\n";
        mStructure->GetLogger() << "**        Start time stepping             **" << "\n";
        mStructure->GetLogger() << "********************************************" << "\n\n";

        while (curTime < rTimeDelta)
        {


            mStructure->DofTypeActivateAll();

            // calculate Delta_BRhs and Delta_ExtForce
            auto bRHS           = UpdateAndGetAndMergeConstraintRHS(mTime, lastConverged_dof_dt0);

            curTime += timeStep;
            SetTimeAndTimeStep(curTime, timeStep, rTimeDelta);     //check whether harmonic excitation, check whether curTime is too close to the time data
            mStructure->SetTime(curTime);

            auto deltaBRHS = UpdateAndGetConstraintRHS(curTime) - bRHS;
            extForce = CalculateCurrentExternalLoad(curTime);

            for (const auto& activeDofSet : mStepActiveDofs)
            {
                mStructure->DofTypeSetIsActive(activeDofSet);

                // ******************************************************
                mStructure->Evaluate(inputMap, evalInternalGradientAndHessian0);
                // ******************************************************



                // Calculate the right hand side
                residual    = extForce - intForce;
                residual.J -= hessian0.JK * deltaBRHS;

                VectorXd rhs;
                rhs.setZero(numTotalDofs);
                rhs.head(numActiveDofs)     = residual.J.Export();
                rhs.tail(numRigidBodyModes) = residual.K.Export();
                rhs -= Btrans * lambda;

                // K_{JJ}^{-1}
                mSolver.compute(hessian0.JJ.ExportToEigenSparseMatrix());

                delta_dof_dt0   = FetiSolve(rhs, activeDofSet, deltaLambda);


                /// \todo This is necessary for displacement controlled loading but I am not sure if it is needed in TFETI
//                delta_dof_dt0.K = deltaBRHS;

//                StructureOutputBlockVector temporary(delta_dof_dt0);

//                delta_dof_dt0.J[Node::eDof::DISPLACEMENTS].setZero(numActiveDofs);
//                delta_dof_dt0.J[Node::eDof::DISPLACEMENTS] = temporary.J.Export().head(numActiveDofs);
//                delta_dof_dt0.K[Node::eDof::DISPLACEMENTS] = temporary.J.Export().tail(numRigidBodyModes);




                dof_dt0 = lastConverged_dof_dt0 + delta_dof_dt0;

                lambda  = lastConverged_lambda + deltaLambda;

                mStructure->NodeMergeDofValues(dof_dt0);

                // ******************************************************
                mStructure->Evaluate(inputMap, evalInternalGradient);
                // ******************************************************


                residual =  extForce - intForce;
                residual.J = BlockFullVector<double>(residual.J.Export() - (Btrans * lambda).head(numActiveDofs), dofStatus);

                normResidual = residual.J.CalculateInfNorm();
                for (const auto& dof : activeDofSet)
                    boost::mpi::all_reduce(world,boost::mpi::inplace(normResidual[dof]),std::plus<double>());


                mStructure->GetLogger() << "Residual: \t" << normResidual << "\n\n";


                int iteration = 0;
                while( not(normResidual < mToleranceResidual) and iteration < mMaxNumIterations)
                {
                    // this block should not be executed for a linear elastic problem and needs to be reworked
                    std::exit(EXIT_FAILURE);

//                    // ******************************************************
//                    mStructure->Evaluate(inputMap, evalHessian0);
//                    // ******************************************************

//                    mSolver.compute(hessian0.JJ.ExportToEigenSparseMatrix());

//                    delta_dof_dt0 = FetiSolve(residual.J.Export(), activeDofSet, deltaLambda);
//                    delta_dof_dt0.K = cmat*delta_dof_dt0.J*(-1.);


//                    dof_dt0 += delta_dof_dt0;
//                    lambda += deltaLambda;
//                    mStructure->NodeMergeDofValues(dof_dt0);

//                    // ******************************************************
//                    mStructure->Evaluate(inputMap, evalInternalGradient);
//                    // ******************************************************

//                    residual =  extForce - intForce;
//                    residual.J = BlockFullVector<double>(residual.J.Export() - Btrans * lambda, dofStatus);

//                    normResidual = residual.J.CalculateInfNorm();
//                    for (const auto& dof : activeDofSet)
//                        boost::mpi::all_reduce(world,boost::mpi::inplace(normResidual[dof]),std::plus<double>());


//                    PrintInfoIteration(normResidual,iteration);
//                    iteration++;

                } // end of while(normResidual<mToleranceForce && iteration<mMaxNumIterations)


                world.barrier();

                if (normResidual < mToleranceResidual)
                {
                    //converged solution
                    mStructure->ElementTotalUpdateStaticData();

                    //store converged step
                    lastConverged_dof_dt0 = dof_dt0;
                    lastConverged_lambda  = lambda;

                    prevResidual = residual;

                    mStructure->NodeMergeDofValues(dof_dt0);

                    //update structure time
                    mStructure->SetPrevTime(curTime);

                    mTime+=timeStep;


                    mStructure->GetLogger() << "Convergence after " << iteration << " iterations at time " << mTime << " (timestep " << timeStep << ").\n";
                    mStructure->GetLogger() << "Residual: \t" << normResidual << "\n\n";


                    //eventually increase next time step
                    if (mAutomaticTimeStepping && iteration<0.25*mMaxNumIterations)
                    {
                        timeStep*=1.5;
                        if (timeStep>mMaxTimeStep)
                            timeStep = mMaxTimeStep;
                    }

                    //perform Postprocessing
                    PostProcess(prevResidual);

                    if (mCallback && mCallback->Exit(*mStructure))
                        return NuTo::eError::SUCCESSFUL;

                }
                else
                {
                    if (structure->mRank == 0)
                        mStructure->GetLogger() << "No convergence with timestep " << timeStep << "\n\n";

                    //no convergence
                    if (mAutomaticTimeStepping)
                    {
                        //no convergence, reduce the time step and start from scratch
                        curTime -= timeStep;
                        timeStep *= 0.5;
                        if (timeStep < mMinTimeStep) {
                            mStructure->GetLogger() << "The minimal time step achieved, the actual time step is " << timeStep << "\n";
                            throw MechanicsException(__PRETTY_FUNCTION__, "No convergence, the current time step is too short.");
                        }
                    }
                    else
                    {
                        throw MechanicsException(__PRETTY_FUNCTION__, "No convergence with the current maximum number of iterations, either use automatic time stepping, reduce the time step or the minimal line search cut back factor.");
                    }

                }

            } // end for mStepActiveDofs


        } // end while



    }
    catch (MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, " ERROR performing FETI.");
        throw;
    }


    return NuTo::eError::SUCCESSFUL;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

