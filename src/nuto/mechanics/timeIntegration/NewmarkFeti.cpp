#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/base/Timer.h"


#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/StructureFETI.h"
#include "nuto/mechanics/structures/StructureBaseEnum.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/mechanics/structures/StructureOutputDummy.h"
#include <boost/mpi.hpp>

#include "nuto/base/CallbackInterface.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include <eigen3/Eigen/Dense>
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "nuto/mechanics/timeIntegration/NewmarkFeti.h"

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



int NuTo::NewmarkFeti::BiCgStab(const NuTo::NewmarkFeti::MatrixXd &projection, NuTo::NewmarkFeti::VectorXd &x, const NuTo::NewmarkFeti::VectorXd &rhs)
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
    VectorXd r0 = r;

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

    const double tol    = mCpgTolerance;
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

    VectorXd Ax = B * mSolver.solve(Btrans*x);
    boost::mpi::all_reduce(world,boost::mpi::inplace(Ax.data()), Ax.size(),std::plus<double>());

    //initial residual
    VectorXd r = rhs - Ax;
    VectorXd z;

    //initial projected search direction
    VectorXd w = projection * r;
    VectorXd p = w;

    // precondition
    p =  projection * mLocalPreconditioner * w;
    boost::mpi::all_reduce(world,boost::mpi::inplace(p.data()), p.size(),std::plus<double>());
//    p = w;

    const double rhs_sqnorm = rhs.squaredNorm();
    const double threshold = mCpgTolerance * mCpgTolerance * rhs_sqnorm;

    double absNew = w.dot(p);
    int iteration = 0;
    while(iteration < mCpgMaxIterations)
    {
        // at every iteration i the r has to be recomputd which is quite expensive
        Ap.noalias() = B * mSolver.solve(Btrans*p);
        boost::mpi::all_reduce(world,boost::mpi::inplace(Ap.data()), Ap.size(),std::plus<double>());

        // step size
        double alpha = absNew / p.dot(Ap);

        // update solution
        x += alpha * p;

        // update projected residual
        w -= alpha * projection * Ap;

        // precondition
        z =  projection * mLocalPreconditioner * w;
        boost::mpi::all_reduce(world,boost::mpi::inplace(z.data()), z.size(),std::plus<double>());
//        z = w;

        // project the r to guarantee the constraints
//        w = projection * r;

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



NuTo::StructureOutputBlockVector NuTo::NewmarkFeti::FetiSolve(BlockFullVector<double> residual_mod, const std::set<Node::eDof> &activeDofSet, VectorXd &deltaLambda)
{

    boost::mpi::communicator world;

    StructureFETI* structure = static_cast<StructureFETI*>(mStructure);

    const SparseMatrix& B                    = structure->GetConnectivityMatrix();
    const SparseMatrix& Btrans               = B.transpose();
    const int& numLagrangeMultipliers        = B.rows();

    structure->mNumRigidBodyModes           = mSolver.cols() -  mSolver.rank();
    const int& numRigidBodyModesLocal       = structure->mNumRigidBodyModes;

    const int numRigidBodyModesGlobal = boost::mpi::all_reduce(world,numRigidBodyModesLocal, std::plus<int>());

    world.barrier();

    structure->mInterfaceRigidBodyModes.resize(numLagrangeMultipliers, numRigidBodyModesLocal);

    std::cout << "Rank: " << structure->mRank << " => number of rigid body modes " << numRigidBodyModesLocal << std::endl;
    world.barrier();

    if (structure->mNumRigidBodyModes > 0)
    {
        // extract the last columns of the Q matrix which represents the null space of K
        Eigen::MatrixXd rbHelper;
        rbHelper.setZero(mSolver.rows(), numRigidBodyModesLocal);
        rbHelper.bottomLeftCorner(numRigidBodyModesLocal,numRigidBodyModesLocal) = Eigen::MatrixXd::Identity(numRigidBodyModesLocal, numRigidBodyModesLocal);
        structure->mRigidBodyModes = mSolver.matrixQ() * rbHelper;
        std::cout << "Rank: " << structure->mRank << " => K*R.maxCoeff() " << mSolver.solve(structure->mRigidBodyModes).maxCoeff() << std::endl;

        structure->mInterfaceRigidBodyModes = B * structure->mRigidBodyModes;

    }



    const auto& rigidBodyModes              = structure->GetRigidBodyModes();
    auto& interfaceRigidBodyModes           = structure->GetInterfaceRigidBodyModes();

    world.barrier();


    // G.size() = (number of Lagrange multipliers) x (total number of rigid body modes)
    MatrixXd G          =   GatherInterfaceRigidBodyModes(interfaceRigidBodyModes, numRigidBodyModesGlobal);

    // GtransGinv.size() = (total number of rigid body modes) x (total number of rigid body modes)
    MatrixXd GtransGinv =   (G.transpose() * G).inverse();

    // The projection matrix guarantees that the solution for lambda satisfies the constraint at every iteration
    MatrixXd projection =   MatrixXd::Identity(G.rows(), G.rows()) - G * GtransGinv * G.transpose();

    VectorXd rigidBodyForceVectorLocal(numRigidBodyModesLocal);
    if (structure->mNumRigidBodyModes > 0)
        rigidBodyForceVectorLocal = rigidBodyModes.transpose() * residual_mod.Export();

    VectorXd rigidBodyForceVectorGlobal = GatherRigidBodyForceVector(rigidBodyForceVectorLocal, numRigidBodyModesGlobal);

    world.barrier();

    VectorXd displacementGap    = B * mSolver.solve(residual_mod.Export());

    boost::mpi::all_reduce(world,boost::mpi::inplace(displacementGap.data()), displacementGap.size(),std::plus<double>());

    deltaLambda = G * GtransGinv * rigidBodyForceVectorGlobal;






    VectorXd rhs    = displacementGap;
    VectorXd x      = deltaLambda;

    VectorXd alphaGlobal =VectorXd::Zero(numRigidBodyModesGlobal);


    if (structure->mRank == 0)
    {
        structure->GetLogger() << "\n*************************************\n";
        structure->GetLogger() << "       Start Interface Problem           ";
        structure->GetLogger() << "\n*************************************\n\n";
    }


    world.barrier();

    int iteration = CPG(projection,x,rhs);
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

    VectorXd tmp = B * mSolver.solve(Btrans*x);
    boost::mpi::all_reduce(world,boost::mpi::inplace(tmp.data()), tmp.size(),std::plus<double>());

    alphaGlobal = GtransGinv * G.transpose() * (displacementGap - tmp);

    world.barrier();

    std::vector<int> recvCount;
    std::vector<int> displs;
    MpiGatherRecvCountAndDispls(recvCount, displs, numRigidBodyModesLocal);
    VectorXd alphaLocal = alphaGlobal.segment(displs[structure->mRank], numRigidBodyModesLocal);

    StructureOutputBlockVector  delta_dof_dt0(structure->GetDofStatus());
    VectorXd delta_dof_active = mSolver.solve(residual_mod.Export() - Btrans * deltaLambda);


    if (structure->mNumRigidBodyModes > 0)
        delta_dof_active -= rigidBodyModes * alphaLocal;

int offset = 0;
for (const auto& dof : activeDofSet)
{
    int numActiveDofs = structure->GetNumActiveDofs(dof);
    delta_dof_dt0.J[dof] = delta_dof_active.segment(offset, structure->GetNumActiveDofs(dof));
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

        structure->AssembleBoundaryDofIds();
        const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& boundaryDofIds       = structure->mBoundaryDofIds;

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


        /*---------------------------------*\
            |        Allocate Variables         |
            \*---------------------------------*/

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

        /*---------------------------------*\
            |    Declare and fill Output Maps   |
            \*---------------------------------*/

        // Declare output maps

        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradient;
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradientAndHessian0;
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalHessian0;

        evalInternalGradient                [eStructureOutput::INTERNAL_GRADIENT]   = &intForce;

        evalInternalGradientAndHessian0     [eStructureOutput::INTERNAL_GRADIENT]   = &intForce;
        evalInternalGradientAndHessian0     [eStructureOutput::HESSIAN0]            = &hessian0;

        evalHessian0                        [eStructureOutput::HESSIAN0]            = &hessian0;

        /*---------------------------------*\
            |    Declare and fill Input map     |
            \*---------------------------------*/

        ConstitutiveInputMap inputMap;
        inputMap[Constitutive::eInput::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(
                    eCalculateStaticData::EULER_BACKWARD);

        ExtractDofValues(lastConverged_dof_dt0, lastConverged_dof_dt1, lastConverged_dof_dt2);

        UpdateAndGetAndMergeConstraintRHS(curTime, lastConverged_dof_dt0);

        PostProcess(residual);

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

                residual    = extForce - intForce;
                residual.J -= hessian0.JK * deltaBRHS;
                residual.J  = BlockFullVector<double>(residual.J.Export() - Btrans * lambda, dofStatus);

                SparseMatrix tangentStiffnessMatrix = hessian0.JJ.ExportToEigenSparseMatrix();

                mSolver.compute(tangentStiffnessMatrix);
                mLocalPreconditioner = B * tangentStiffnessMatrix * Btrans;

                delta_dof_dt0   = FetiSolve(residual.J, activeDofSet, deltaLambda);
                delta_dof_dt0.K = deltaBRHS;

                dof_dt0 = lastConverged_dof_dt0 + delta_dof_dt0;
                lambda  = lastConverged_lambda + deltaLambda;

                mStructure->NodeMergeDofValues(dof_dt0);

                // ******************************************************
                mStructure->Evaluate(inputMap, evalInternalGradient);
                // ******************************************************


                residual    = extForce - intForce;
                residual.J  = BlockFullVector<double>(residual.J.Export() - Btrans * lambda, dofStatus);

                normResidual = residual.J.CalculateInfNorm();
                for (const auto& dof : activeDofSet)
                    boost::mpi::all_reduce(world,boost::mpi::inplace(normResidual[dof]),std::plus<double>());


                if (structure->mRank == 0)
                {
                    mStructure->GetLogger() << "normResidual= \t" << normResidual << "\n" <<"tolerance= \t" << mToleranceResidual <<  "\n";

                }

                world.barrier();


                int iteration = 0;
                while( not(normResidual < mToleranceResidual) and iteration < mMaxNumIterations)
                {

                    // ******************************************************
                    mStructure->Evaluate(inputMap, evalHessian0);
                    // ******************************************************

                    mSolver.compute(hessian0.JJ.ExportToEigenSparseMatrix());

                    delta_dof_dt0 = FetiSolve(residual.J, activeDofSet, deltaLambda);
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
        throw e;
    }


    return NuTo::eError::SUCCESSFUL;

}
