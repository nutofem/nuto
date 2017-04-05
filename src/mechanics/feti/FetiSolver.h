//
// Created by phuschke on 2/8/17.
//

#pragma once

#include <mpi.h>
#include <boost/mpi.hpp>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>


class FetiSolver {
private:

    using MatrixXd = Eigen::MatrixXd;
    using VectorXd = Eigen::VectorXd;

public:
    ///
    /// \brief GatherInterfaceRigidBodyModes
    /// \param interfaceRigidBodyModes
    /// \param numRigidBodyModesGlobal
    /// \return
    ///
    Eigen::MatrixXd GatherInterfaceRigidBodyModes(Eigen::MatrixXd& interfaceRigidBodyModes, const int numRigidBodyModesGlobal)
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


    ///
    /// \brief GatherRigidBodyForceVector
    /// \param rigidBodyForceVectorLocal
    /// \param numRigidBodyModesGlobal
    /// \return
    ///
    Eigen::VectorXd GatherRigidBodyForceVector(Eigen::VectorXd &rigidBodyForceVectorLocal, const int numRigidBodyModesGlobal)
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
    ///
    /// \brief MpiGatherRecvCountAndDispls
    /// \param recvCount
    /// \param displs
    /// \param numValues
    ///
    void MpiGatherRecvCountAndDispls(std::vector<int> &recvCount, std::vector<int> &displs, const int numValues)
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

    void AssembleGMatrix(MatrixXd interfaceRigidModes)
    {

    }

    void Solve()
    {

    }
private:

    MatrixXd mG;
    MatrixXd mGtransGinv;


};

