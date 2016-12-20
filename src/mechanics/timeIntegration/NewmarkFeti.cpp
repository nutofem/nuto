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
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "mechanics/timeIntegration/NewmarkFeti.h"

#include <cmath>


Eigen::MatrixXd NuTo::NewmarkFeti::GatherInterfaceRigidBodyModes(Eigen::MatrixXd &interfaceRigidBodyModes, const int numRigidBodyModesGlobal)
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


NuTo::BlockScalar NuTo::NewmarkFeti::CalculateNormResidual(BlockFullVector<double> &residual_mod, const std::set<Node::eDof> &activeDofSet)
{
    StructureFETI* structure = static_cast<StructureFETI*>(mStructure);
    for (const auto& interface : structure->mInterfaces)
        for (const auto& nodePair: interface.mNodeIdsMap)
        {
            for (const auto& dof : activeDofSet)
            {
                const std::vector<int> dofIds = structure->NodeGetDofIds(nodePair.second, dof);
                for (const auto& id : dofIds)
                    residual_mod[dof][id] = 0;

            }

        }

    BlockScalar normResidual = residual_mod.CalculateInfNorm();

    for (const auto& dof : activeDofSet)
        MPI_Allreduce(MPI_IN_PLACE,  &normResidual[dof], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);





    return normResidual;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int NuTo::NewmarkFeti::BiCgStab(const MatrixXd &projection, VectorXd &x, const VectorXd &rhs)
{

    StructureFETI* structure = static_cast<StructureFETI*>(mStructure);
    const SparseMatrix B      = structure->GetConnectivityMatrix();
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

    const SparseMatrix B      = structure->GetConnectivityMatrix();
    const SparseMatrix Btrans = B.transpose();


    std::cout << "Rank: " << structure->mRank << " => mSolver " << mSolver.rows() << "\t B.rows() \t" << B.rows() << "\t B.cols() \t" << B.cols() << std::endl;
    world.barrier();

    const int numActiveDofs     = mStructure->GetNumActiveDofs(Node::eDof::DISPLACEMENTS);
    const int numRigidBodyModes = structure->mNumRigidBodyModes;

    VectorXd tmp;
    tmp.setZero(numActiveDofs + numRigidBodyModes);
    tmp.head(numActiveDofs) = mSolver.solve( (Btrans*x).head(numActiveDofs));
    VectorXd Ax = B * tmp;


//    std::cout << "Rank: " << structure->mRank << " tmp " << tmp << std::endl;
    world.barrier();


    boost::mpi::all_reduce(world,boost::mpi::inplace(Ax.data()), Ax.size(),std::plus<double>());

    world.barrier();
    //initial residual
    VectorXd r = rhs - Ax;
    VectorXd z;

    //initial projected search direction
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

        if (structure->mRank == 0)
            std::cout << "CPG rel. error = " << w.squaredNorm()/rhs_sqnorm << "\t at iteration = " << iteration << "/" << mCpgMaxIterations << std::endl;


        ++iteration;


    }

    return iteration;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


NuTo::StructureOutputBlockVector NuTo::NewmarkFeti::FetiSolve(VectorXd residual_mod, const std::set<Node::eDof> &activeDofSet, VectorXd &deltaLambda)
{

    boost::mpi::communicator world;

    StructureFETI* structure = static_cast<StructureFETI*>(mStructure);

    const SparseMatrix& B                    = structure->GetConnectivityMatrix();
    const SparseMatrix& Btrans               = B.transpose();

    const int& numRigidBodyModesLocal       = structure->mNumRigidBodyModes;
    const int numRigidBodyModesGlobal       = boost::mpi::all_reduce(world,numRigidBodyModesLocal, std::plus<int>());

    const auto& rigidBodyModes              = structure->GetRigidBodyModes();

    world.barrier();

    // G.size() = (number of Lagrange multipliers) x (total number of rigid body modes)
    const MatrixXd& G          =   structure->mG;

    // GtransGinv.size() = (total number of rigid body modes) x (total number of rigid body modes)
    const MatrixXd GtransGinv =   (G.transpose() * G).inverse();


    if (structure->mRank == 0)
        std::cout << "residual_mod \n" << residual_mod << std::endl;

    world.barrier();
    if (structure->mRank == 1)
        std::cout << "residual_mod \n" << residual_mod << std::endl;

    world.barrier();
    std::cin.get();

    // R_s^T f
    VectorXd rigidBodyForceVectorLocal = rigidBodyModes.transpose() * residual_mod;

    // e = [ R_1^T f; R_2^T f; ...; R_{N_s}^T f ]
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
    VectorXd tempVec;
    tempVec.setZero(mStructure->GetNumTotalDofs());
    tempVec.head(mStructure->GetNumTotalActiveDofs()) = mSolver.solve(residual_mod.head(mStructure->GetNumTotalActiveDofs()));
    VectorXd displacementGap    = B * tempVec;
    boost::mpi::all_reduce(world,boost::mpi::inplace(displacementGap.data()), displacementGap.size(),std::plus<double>());


    VectorXd rhs    = displacementGap;
    VectorXd x      = deltaLambda;




    if (structure->mRank == 0)
        std::cout << "rhs: " << rhs << std::endl;

    if (structure->mRank == 0)
        std::cout << "x: " << x << std::endl;


    if (structure->mRank == 0)
    {
        structure->GetLogger() << "\n*************************************\n";
        structure->GetLogger() << "       Start Interface Problem           ";
        structure->GetLogger() << "\n*************************************\n\n";
    }

    world.barrier();

    int iteration = CPG(structure->mProjectionMatrix,x,rhs);
    //int iteration = BiCgStab(projection,x,rhs);

    if (iteration >= mCpgMaxIterations)
        throw MechanicsException(__PRETTY_FUNCTION__,"Maximum number of iterations exceeded.");

    if (structure->mRank == 0)
    {
        structure->GetLogger() << "\n*************************************\n";
        structure->GetLogger() << "       End Interface Problem             ";
        structure->GetLogger() << "\n*************************************\n\n";
    }



    world.barrier();
    deltaLambda = x;



    //*****************************************
    //
    //         | K_11^-1 * lambda |
    // tmp = B |  0               |
    //
    //*****************************************
    tempVec.head(structure->GetNumActiveDofs(Node::eDof::DISPLACEMENTS)) = mSolver.solve((Btrans*x).head(structure->GetNumActiveDofs(Node::eDof::DISPLACEMENTS)));

    VectorXd tmp = B * tempVec;
    boost::mpi::all_reduce(world,boost::mpi::inplace(tmp.data()), tmp.size(),std::plus<double>());

    VectorXd alphaGlobal =VectorXd::Zero(numRigidBodyModesGlobal);
    alphaGlobal = GtransGinv * G.transpose() * (displacementGap - tmp);


    world.barrier();

    std::vector<int> recvCount;
    std::vector<int> displs;
    MpiGatherRecvCountAndDispls(recvCount, displs, numRigidBodyModesLocal);
    VectorXd alphaLocal = alphaGlobal.segment(displs[structure->mRank], numRigidBodyModesLocal);


    std::cout << "alphaLocal \t" << alphaLocal << std::endl;

    VectorXd delta_dof_active;
    delta_dof_active.setZero(mStructure->GetNumTotalDofs());
    delta_dof_active    = mSolver.solve( (residual_mod - Btrans * deltaLambda).head(structure->GetNumTotalActiveDofs()));
    delta_dof_active   -= rigidBodyModes * alphaLocal;

    StructureOutputBlockVector  delta_dof_dt0(structure->GetDofStatus());
    int offset = 0;
    for (const auto& dof : activeDofSet)
    {
        int numActiveDofs       = structure->GetNumActiveDofs(dof);
        delta_dof_dt0.J[dof]    = delta_dof_active.segment(offset, structure->GetNumActiveDofs(dof)+numRigidBodyModesLocal);
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
        structure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

        structure->AssembleConnectivityMatrix();
        const SparseMatrix& B       = structure->GetConnectivityMatrix();
        const SparseMatrix& Btrans  = B.transpose();

        VectorXd deltaLambda            = VectorXd::Zero(B.rows());
        VectorXd lambda                 = VectorXd::Zero(B.rows());
        VectorXd lastConverged_lambda   = VectorXd::Zero(B.rows());

        double curTime  = mTime;
        double timeStep = mTimeStep;
        mStructure->SetPrevTime(curTime);
        mStructure->SetTime(curTime);

        mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

        CalculateStaticAndTimeDependentExternalLoad();

        const DofStatus& dofStatus = mStructure->GetDofStatus();


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
        const auto& cmat = mStructure->GetConstraintMatrix();

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

        std::cout << "after first post process" << std::endl;

        //********************************************
        //        Calculate rigid body modes
        //********************************************

        const int numActiveDofs     = mStructure->GetNumTotalActiveDofs();
        const int numRigidBodyModes = structure->mNumRigidBodyModes;
        const int numTotalDofs      = mStructure->GetNumTotalDofs();

        mStructure->Evaluate(inputMap, evalHessian0);

        SparseMatrix hessian0_JJ = hessian0.JJ.ExportToEigenSparseMatrix();
        SparseMatrix hessian0_JK = hessian0.JK.ExportToEigenSparseMatrix();

        mSolver.compute(hessian0_JJ);

        // The dimension of the rigid body modes depends on the number of types
        structure->mRigidBodyModes.setZero(numTotalDofs, numRigidBodyModes);

        std::cout << "structure->mRank \t" << structure->mRank << "\n structure->mRigidBodyModes.rows() \t" << structure->mRigidBodyModes.rows() << "\n structure->mRigidBodyModes.cols() \t" << structure->mRigidBodyModes.cols() << std::endl;

        if (structure->IsFloating())
        {
            structure->mRigidBodyModes.topRows(numActiveDofs) = mSolver.solve(hessian0_JK);
            structure->mRigidBodyModes *= -1; // this is stupid please fix asap
        }

        structure->mRigidBodyModes.bottomRows(numRigidBodyModes) = Eigen::MatrixXd::Identity(numRigidBodyModes,numRigidBodyModes);

        //********************************************
        //        Calculate interface rigid body modes
        //********************************************

        const int& numLagrangeMultipliers        = B.rows();

        structure->mInterfaceRigidBodyModes.resize(numLagrangeMultipliers, numRigidBodyModes);

        world.barrier();

        structure->mInterfaceRigidBodyModes = B * structure->mRigidBodyModes;


        //********************************************
        //        Calculate projection
        //********************************************

        const int numRigidBodyModesGlobal = boost::mpi::all_reduce(world,numRigidBodyModes, std::plus<int>());

        std::cout << "numRigidBodyModesGlobal \t" << numRigidBodyModesGlobal << std::endl;

        const MatrixXd G       =   GatherInterfaceRigidBodyModes(structure->mInterfaceRigidBodyModes, numRigidBodyModesGlobal);

        std::cout << "structure->mRank \t" << structure->mRank << "\n G.rows() \t" << G.rows() << "\n G.cols() \t" << G.cols() << std::endl;

        structure->mProjectionMatrix = MatrixXd::Identity(G.rows(), G.rows()) - G * (G.transpose() * G).inverse() * G.transpose();
        structure->mG = G;

        //********************************************
        //        Start time stepping
        //********************************************
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

                std::cout << "Rank: \t" << structure->mRank << "\n extForce \n" << extForce.J.Export().Norm() << std::endl;
                std::cout << "Rank: \t" << structure->mRank << "\n residual \n" << residual.J.Export().Norm() << std::endl;

                VectorXd rhs;
                rhs.setZero(numTotalDofs);

                rhs.head(numActiveDofs)     = residual.J.Export();
                rhs.tail(numRigidBodyModes) = residual.K.Export();
                rhs -= Btrans * lambda;

                // K_{JJ}^{-1}
                hessian0_JJ = hessian0.JJ.ExportToEigenSparseMatrix();
                mSolver.compute(hessian0_JJ);


                std::cout << "Rank: \t"         << structure->mRank                     << std::endl
                          << "Active dofs: \t"  << mStructure->GetNumTotalActiveDofs()  << std::endl
                          << "Total dofs: \t"   << structure->GetNumTotalDofs()         << std::endl;

                world.barrier();

                std::cout << "Rank: \t" << structure->mRank << " Before FetiSolve" << std::endl;

                delta_dof_dt0   = FetiSolve(rhs, activeDofSet, deltaLambda);

                std::cout << "Rank: \t" << structure->mRank << " After FetiSolve" << std::endl;
                world.barrier();

                std::cout << "Rank: \t"         << structure->mRank                     << std::endl
                          << "\t delta_dof_dt0.J.Export().rows(): \t"  << delta_dof_dt0.J.Export().rows()      << std::endl
                          << "\t delta_dof_dt0.K.Export().rows(): \t"  << delta_dof_dt0.K.Export().rows()      << std::endl;

                delta_dof_dt0.K = deltaBRHS;

                StructureOutputBlockVector temporary(delta_dof_dt0);

                delta_dof_dt0.J[Node::eDof::DISPLACEMENTS].setZero(numActiveDofs);
                delta_dof_dt0.J[Node::eDof::DISPLACEMENTS] = temporary.J.Export().head(numActiveDofs);
                delta_dof_dt0.K[Node::eDof::DISPLACEMENTS] = temporary.J.Export().tail(numRigidBodyModes);


                std::cout << "Rank: \t"         << structure->mRank                     << std::endl
                          << "\t delta_dof_dt0.J.Export().rows(): \t"  << delta_dof_dt0.J.Export().rows()      << std::endl
                          << "\t delta_dof_dt0.K.Export().rows(): \t"  << delta_dof_dt0.K.Export().rows()      << std::endl;

                dof_dt0 = lastConverged_dof_dt0 + delta_dof_dt0;

                std::cout << "Rank: \t"         << structure->mRank                     << std::endl
                          << "\t dof_dt0: \t"   << dof_dt0                              << std::endl;

                lambda  = lastConverged_lambda + deltaLambda;


                mStructure->NodeMergeDofValues(dof_dt0);

                // ******************************************************
                mStructure->Evaluate(inputMap, evalInternalGradient);
                // ******************************************************

                mTime+=timeStep;
                PostProcess(prevResidual);

                residual                        =   extForce - intForce;

                rhs.head(numActiveDofs)         =   residual.J.Export();
                rhs.tail(numRigidBodyModes)     =   residual.K.Export();
                rhs                            -=   Btrans * lambda;

                std::cout << "Rank: \t" << structure->mRank << " norm "  << rhs.norm() << std::endl;

                normResidual = residual.J.CalculateInfNorm();

                mStructure->GetLogger() << "Rank " << structure->mRank <<" normResidual= \t" << normResidual << "\n" <<"tolerance= \t" << mToleranceResidual <<  "\n";


                for (const auto& dof : activeDofSet)
                    boost::mpi::all_reduce(world,boost::mpi::inplace(normResidual[dof]),std::plus<double>());

                if (structure->mRank == 0)
                    std::cout << "Rank: \t" << structure->mRank << "\n normResidual "  << normResidual << std::endl;

                world.barrier();
                std::cin.get();

                int iteration = 0;
                while( not(normResidual < mToleranceResidual) and iteration < mMaxNumIterations)
                {

                    // ******************************************************
                    mStructure->Evaluate(inputMap, evalHessian0);
                    // ******************************************************

                    mSolver.compute(hessian0.JJ.ExportToEigenSparseMatrix());

                    delta_dof_dt0 = FetiSolve(residual.J.Export(), activeDofSet, deltaLambda);
                    delta_dof_dt0.K = cmat*delta_dof_dt0.J*(-1.);


                    dof_dt0 += delta_dof_dt0;
                    lambda += deltaLambda;
                    mStructure->NodeMergeDofValues(dof_dt0);

                    // ******************************************************
                    mStructure->Evaluate(inputMap, evalInternalGradient);
                    // ******************************************************

                    residual =  extForce - intForce;
                    residual.J = BlockFullVector<double>(residual.J.Export() - Btrans * lambda, dofStatus);

                    normResidual = residual.J.CalculateInfNorm();
                    for (const auto& dof : activeDofSet)
                        boost::mpi::all_reduce(world,boost::mpi::inplace(normResidual[dof]),std::plus<double>());


                    PrintInfoIteration(normResidual,iteration);
                    iteration++;

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

                    if (structure->mRank == 0)
                    {
                        mStructure->GetLogger() << "Convergence after " << iteration << " iterations at time " << mTime << " (timestep " << timeStep << ").\n";
                        mStructure->GetLogger() << "Residual: \t" << normResidual << "\n\n";
                    }

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
