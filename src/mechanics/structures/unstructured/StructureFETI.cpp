#include <mpi/mpi.h>

#include <boost/mpi.hpp>
#include <json/json.h>


#include "mechanics/structures/unstructured/StructureFETI.h"
#include "mechanics/nodes/NodeBase.h"

#include "mechanics/MechanicsException.h"



#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"

#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/sections/SectionEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "base/ErrorEnum.h"

#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/structures/StructureBaseEnum.h"

using std::cout;
using std::endl;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;
using NuTo::Constitutive::eConstitutiveType;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Node::eDof;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Interpolation::eShapeType;
using NuTo::eGroupId;
using NuTo::eVisualizeWhat;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

NuTo::StructureFETI::StructureFETI(int rDimension):
    Structure(rDimension),
    mRank(boost::mpi::communicator().rank()),
    mNumProcesses(boost::mpi::communicator().size())
{

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//void NuTo::StructureFETI::AssembleBoundaryDofIds()
//{

//    const auto& dofTypes = GetDofStatus().GetDofTypes();

//    int numActiveDofs = 0;

//    for (const auto& dofType : dofTypes)
//        numActiveDofs           += GetNumActiveDofs(dofType);

//    mBoundaryDofIds.setZero(numActiveDofs);

//    int offset = 0;
//    for (const auto& dofType : dofTypes)
//    {
//        for (const auto& nodeId : mSubdomainBoundaryNodeIds)
//        {
//            const std::vector<int> dofIds = NodeGetDofIds(nodeId, dofType);

//            for (const auto& dofId : dofIds)
//            {
//                if(dofId < GetNumActiveDofs(dofType))
//                    mBoundaryDofIds.diagonal()(dofId + offset) = 1;
//            }

//        }
//        offset += GetNumActiveDofs(dofType);
//    }



//}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::AssembleConnectivityMatrix()
{

    const auto& dofTypes    = GetDofStatus().GetDofTypes();
    const int numTotalDofs  = GetNumTotalDofs();

    int numLagrangeMultipliers  = 0;

    for (const auto& dofType : dofTypes)
    {
        numLagrangeMultipliers  += GetDofDimension(dofType) * mNumInterfaceNodesTotal;
    }

    numLagrangeMultipliers  += mNumTotalBoundaryDofIds;
    mConnectivityMatrix.resize(numLagrangeMultipliers, numTotalDofs);

    int offsetRows = 0;
    int offsetCols = 0;
    for (const auto& dofType : dofTypes)
    {
        for (const auto& interface : mInterfaces)
            for (const auto& nodePair : interface.mNodeIdsMap)
            {

                const std::vector<int> dofVector    = NodeGetDofIds(nodePair.second, dofType);
                const int globalIndex               = nodePair.first * dofVector.size();

                for (unsigned i = 0; i < dofVector.size(); ++i)
                {
                    // remove the mNumRigidBodyModes because it is only associated with displacements
                    if (dofVector[i] < GetNumActiveDofs(dofType) + mNumRigidBodyModes)
                        mConnectivityMatrix.insert(globalIndex + i + offsetRows , dofVector[i] + offsetCols) = interface.mValue;
                }
            }

        offsetRows += GetDofDimension(dofType) * mNumInterfaceNodesTotal;
        // remove the mNumRigidBodyModes because it is only associated with displacements
        offsetCols += GetNumActiveDofs(dofType) + mNumRigidBodyModes;
    }

    int globalIndex = numLagrangeMultipliers - mNumTotalBoundaryDofIds + mGlobalStartIndexBoundaryDofIds;
    for (const int& id :  mBoundaryDofIds)
    {
        mConnectivityMatrix.insert(globalIndex,id) = 1;
        ++globalIndex;
    }



    GetLogger() << "Number of Lagrange multipliers:  \t"        << numLagrangeMultipliers       << "\n\n";
    GetLogger() << "Total number of interface nodes: \t"        << mNumInterfaceNodesTotal      << "\n\n";
    GetLogger() << "Total number of boundary dofs:   \t"        << mNumTotalBoundaryDofIds      << "\n\n";
    GetLogger() << "Total number of dofs:            \t"        << GetNumTotalDofs()            << "\n\n";
    GetLogger() << "Total number of active dofs:     \t"        << GetNumTotalActiveDofs()      << "\n\n";
    GetLogger() << "Total number of dependent dofs:  \t"        << GetNumTotalDependentDofs()   << "\n\n";



}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::ImportMeshJson(std::string rFileName, const int interpolationTypeId)
{

    Json::Value root;
    Json::Reader reader;

    std::ifstream file(rFileName.c_str(), std::ios::in);

    if(not reader.parse(file,root, false))
        throw MechanicsException(__PRETTY_FUNCTION__, "Error parsing mesh file.");


    // only supports nodes.size() == 1
    for (auto const& nodes : root["Nodes"])
    {
        mNodes.resize(nodes["Coordinates"].size());
        for (unsigned i = 0; i < mNodes.size(); ++i)
        {
            mNodes[i].mCoordinates[0] = nodes["Coordinates"][i][0].asDouble();
            mNodes[i].mCoordinates[1] = nodes["Coordinates"][i][1].asDouble();
            mNodes[i].mCoordinates[2] = nodes["Coordinates"][i][2].asDouble();
            mNodes[i].mId             = nodes["Indices"][i].asInt();
        }


    }


    // only supports elements.size() == 1
    for (auto const& elements : root["Elements"])
    {
        mElements.resize(elements["NodalConnectivity"].size());
        const int elementType = elements["Type"].asInt();


        for (unsigned i = 0; i < mElements.size(); ++i)
        {
            if (elementType == 1)
            {
                mSubdomainBoundaryNodeIds.insert(elements["NodalConnectivity"][i][0].asInt());
                mSubdomainBoundaryNodeIds.insert(elements["NodalConnectivity"][i][1].asInt());

            }
            else if (elementType == 2) // 3 node tri element
            {
                mElements[i].mNodeIds.resize(3);

                mElements[i].mNodeIds[0] = elements["NodalConnectivity"][i][0].asInt();
                mElements[i].mNodeIds[1] = elements["NodalConnectivity"][i][1].asInt();
                mElements[i].mNodeIds[2] = elements["NodalConnectivity"][i][2].asInt();
                mElements[i].mId         = elements["Indices"][i].asInt();

            }
            else if (elementType == 3) // 4 node quad element
            {

                mElements[i].mNodeIds.resize(4);

                mElements[i].mNodeIds[0] = elements["NodalConnectivity"][i][0].asInt();
                mElements[i].mNodeIds[1] = elements["NodalConnectivity"][i][1].asInt();
                mElements[i].mNodeIds[2] = elements["NodalConnectivity"][i][2].asInt();
                mElements[i].mNodeIds[3] = elements["NodalConnectivity"][i][3].asInt();
                mElements[i].mId         = elements["Indices"][i].asInt();
            }
            else
            {
                throw MechanicsException(__PRETTY_FUNCTION__, "Import of element type not implemented. Element type id = " +std::to_string(elementType));
            }
        }
    }


    mInterfaces.resize(root["Interface"].size());
    for (unsigned i = 0; i < mInterfaces.size(); ++i)
    {

        int globalId = root["Interface"][i]["GlobalStartId"][0].asInt();

        mInterfaces[i].mValue = root["Interface"][i]["Value"][0].asInt();

        for (unsigned k = 0; k < root["Interface"][i]["NodeIds"][0].size(); ++k)
        {
            mInterfaces[i].mNodeIdsMap.emplace(globalId, root["Interface"][i]["NodeIds"][0][k].asInt());
            mSubdomainBoundaryNodeIds.insert(root["Interface"][i]["NodeIds"][0][k].asInt());
            globalId++;
        }


    }


    mNumInterfaceNodesTotal = root["NumInterfaceNodes"][0].asInt();

    file.close();

    for (const auto& node : mNodes)
        NodeCreate(node.mId, node.mCoordinates.head(2));

    for (const auto& element : mElements)
        ElementCreate(element.mId,interpolationTypeId, element.mNodeIds);

    ElementTotalConvertToInterpolationType();

    NodeBuildGlobalDofs();

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::CalculateRigidBodyModesTotalFETI()
{

    switch (GetDimension())
    {
    case 1:
    {
        mNumRigidBodyModes = 1;
        const int numTotalDofs = GetNumTotalDofs();
        mRigidBodyModes.setOnes(numTotalDofs,mNumRigidBodyModes);
    }
        break;
    case 2:
    {
        mNumRigidBodyModes = 3;

        const int numTotalDofs = GetNumTotalDofs();

        mRigidBodyModes.setZero(numTotalDofs,mNumRigidBodyModes);

        for (const auto& nodePair : mNodeMap)
        {
            const std::vector<int> dofIds = NodeGetDofIds(nodePair.first, NuTo::Node::eDof::DISPLACEMENTS);

            const Eigen::Matrix<double, 2, 1> coordinates = nodePair.second->Get(NuTo::Node::eDof::COORDINATES);

            mRigidBodyModes.row(dofIds[0]) << 1.,   0., -coordinates[1];
            mRigidBodyModes.row(dofIds[1]) << 0.,   1.,  coordinates[0];
        }


    }
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__,"Structural dimension not supported yet.");
    }


    boost::mpi::communicator world;
    mNumRigidBodyModesTotal       = boost::mpi::all_reduce(world,mNumRigidBodyModes, std::plus<int>());

    GetLogger() << "Number of rigid body modes:        \t"        << mNumRigidBodyModes        << "\n\n";
    GetLogger() << "Total number of rigid body modes:  \t"        << mNumRigidBodyModesTotal   << "\n\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::CreateDummy1D()
{
    mInterfaces.resize(1);

    if (mRank == 0)
    {
        mInterfaces[0].mValue = 1;
        mInterfaces[0].mNodeIdsMap.emplace(0,3);
    }

    if (mRank == 1)
    {
        mInterfaces[0].mValue = -1;
        mInterfaces[0].mNodeIdsMap.emplace(0,0);
    }


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::CalculateRigidBodyModes()
{
    const int numActiveDofs     = GetNumTotalActiveDofs();
    const int numRigidBodyModes = mNumRigidBodyModes;
    const int numTotalDofs      = GetNumTotalDofs();


    ConstitutiveIOMap<Constitutive::eInput> inputMap;
    inputMap[Constitutive::eInput::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(
                eCalculateStaticData::EULER_BACKWARD);

    StructureOutputBlockMatrix  hessian0(GetDofStatus(), true);

    std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalHessian0;
    evalHessian0                        [eStructureOutput::HESSIAN0]            = &hessian0;

    Evaluate(inputMap, evalHessian0);

    SparseMatrix hessian0_JJ = hessian0.JJ.ExportToEigenSparseMatrix();
    SparseMatrix hessian0_JK = hessian0.JK.ExportToEigenSparseMatrix();

    // check if K.JJ is invertible
    Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> testSolver;
    testSolver.compute(hessian0_JJ);

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "structure->mRank \t" << mRank << "\n testSolver.rank(); \t" << testSolver.rank() << "\n testSolver.row \t" << testSolver.rows() << std::endl << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> mSolver;
    mSolver.compute(hessian0_JJ);

    // The dimension of the rigid body modes depends on the number of types
    mRigidBodyModes.setZero(numTotalDofs, numRigidBodyModes);

    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << "structure->mRank \t"                          << mRank
              << "\n structure->mRigidBodyModes.rows() \t"      << mRigidBodyModes.rows()
              << "\n structure->mRigidBodyModes.cols() \t"      << mRigidBodyModes.cols() << std::endl << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    if (IsFloating())
    {
        mRigidBodyModes.topRows(numActiveDofs) = mSolver.solve(hessian0_JK);
        mRigidBodyModes *= -1.; // this is stupid please fix asap
    }

    mRigidBodyModes.bottomRows(numRigidBodyModes) = Eigen::MatrixXd::Identity(numRigidBodyModes,numRigidBodyModes);


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::CheckRigidBodyModes(  const StructureOutputBlockMatrix hessian0,
                                                const double tolerance) const
{

    SparseMatrix hessian0_JJ = hessian0.JJ.ExportToEigenSparseMatrix();
    SparseMatrix hessian0_JK = hessian0.JK.ExportToEigenSparseMatrix();


    Eigen::MatrixXd zeroMatrix0      = hessian0_JJ * mRigidBodyModes.topRows(GetNumTotalActiveDofs());
    zeroMatrix0                     += hessian0_JK * mRigidBodyModes.bottomRows(mNumRigidBodyModes);
    const double norm0 = std::max( zeroMatrix0.maxCoeff(), std::abs(zeroMatrix0.minCoeff()) );

    SparseMatrix hessian0_KJ = hessian0.KJ.ExportToEigenSparseMatrix();
    SparseMatrix hessian0_KK = hessian0.KK.ExportToEigenSparseMatrix();

    Eigen::MatrixXd zeroMatrix1      = hessian0_KJ * mRigidBodyModes.topRows(GetNumTotalActiveDofs());
    zeroMatrix1                     += hessian0_KK * mRigidBodyModes.bottomRows(mNumRigidBodyModes);
    const double norm1 = std::max( zeroMatrix1.maxCoeff(), std::abs(zeroMatrix1.minCoeff()) );

    MPI_Barrier(MPI_COMM_WORLD);

    GetLogger() << "Subdomain: \t"          << mRank
                << "\t norm of K*R: \t"     << norm0   << "\t and \t"   << norm1
                << "\t tolerance: \t"       << tolerance                << "\n\n";

    MPI_Barrier(MPI_COMM_WORLD);

    assert(     (norm0 < tolerance)
            and (norm1 < tolerance)
            and "Calculated rigid body modes are not in the null space of K");

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::CheckStiffnessPartitioning(  const StructureOutputBlockMatrix hessian0,
                                                        const double tolerance) const
{

    SparseMatrix hessian0_JJ = hessian0.JJ.ExportToEigenSparseMatrix();
    Eigen::SparseLU<SparseMatrix> hessian0_JJ_solver(hessian0_JJ);

    SparseMatrix hessian0_JK = hessian0.JK.ExportToEigenSparseMatrix();

    SparseMatrix hessian0_KJ = hessian0.KJ.ExportToEigenSparseMatrix();
    SparseMatrix hessian0_KK = hessian0.KK.ExportToEigenSparseMatrix();


    Eigen::MatrixXd zeroMatrix0      = hessian0_KK  - hessian0_KJ * hessian0_JJ_solver.solve(hessian0_JK);
    const double norm = std::max( zeroMatrix0.maxCoeff(), std::abs(zeroMatrix0.minCoeff()) );


    MPI_Barrier(MPI_COMM_WORLD);

    GetLogger() << "Subdomain: \t"          << mRank
                << "\t norm of ( K_kk - K_kj * inv(K_jj) * K_jk: \t"    << norm
                << "\t tolerance: \t"       << tolerance                << "\n\n";

    MPI_Barrier(MPI_COMM_WORLD);

    assert(     (norm < tolerance)
            and "Stiffness matrix is not partitioned correctly. Check constraints.");

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::CheckProjectionMatrix( const double tolerance ) const
{

    const Eigen::MatrixXd zeroMatrix = mProjectionMatrix - mProjectionMatrix*mProjectionMatrix;

    const double norm = std::max( zeroMatrix.maxCoeff(), std::abs(zeroMatrix.minCoeff()) );
    GetLogger() << "Subdomain: \t"              << mRank
                << "\t norm of P-PP: \t \t"     << norm
                << "\t tolerance: \t"           << tolerance                << "\n\n";

    assert(     (norm < tolerance)
            and "Projection matrix does not satisfy: P - PP = 0");

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::CheckProjectionOfCoarseGrid( const double tolerance ) const
{

    const Eigen::MatrixXd zeroMatrix = mProjectionMatrix.transpose() * mG;

    const double norm = std::max( zeroMatrix.maxCoeff(), std::abs(zeroMatrix.minCoeff()) );

    GetLogger() << "Subdomain: \t"                  << mRank
                << "\t norm of Ptrans * G: \t"      << norm
                << "\t tolerance: \t"               << tolerance                << "\n\n";

    assert(     (norm < tolerance)
            and "Projection of coarse space does not satisfy: Ptrans*G = 0");

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::CalculateInterfaceRigidBodyModes()
{
    mInterfaceRigidBodyModes = mConnectivityMatrix * mRigidBodyModes;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::ApplyConstraintsTotalFeti(const int nodeGroupId)
{
    boost::mpi::communicator world;

    Eigen::VectorXi boundaryNodes = GroupGetMemberIds(nodeGroupId);

    for (int i = 0; i < boundaryNodes.rows(); ++i)
    {
        const int nodeId = boundaryNodes[i];
        std::vector<int> dofIds = NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);

        for (const auto& id : dofIds)
            mBoundaryDofIds.push_back(id);
    }


    // determine the global ids for the constraints

    // recvCount:
    // Contais the number of elements that are received from each process.
    std::vector<int> recvCount(mNumProcesses, 0);

    boost::mpi::all_gather<int>(world,mBoundaryDofIds.size(),recvCount);

    // displs:
    // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i.
    std::vector<int> displs;
    displs.resize(mNumProcesses, 0);
    for (int i = 1; i < mNumProcesses; ++i)
        displs[i] = displs[i-1] + recvCount[i-1];


    const int numLocalBoundaryDofIds = mBoundaryDofIds.size();
    MPI_Allreduce(&numLocalBoundaryDofIds,
                  &mNumTotalBoundaryDofIds,
                  1,
                  MPI_INT,
                  MPI_SUM,
                  MPI_COMM_WORLD);



    GetLogger() << "mBoundaryDofIds.size() \t" << mBoundaryDofIds.size() << "\n \n";
    GetLogger() << "mNumTotalBoundaryDofIds \t" << mNumTotalBoundaryDofIds << "\n \n";


    mGlobalStartIndexBoundaryDofIds = displs[mRank];

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::CalculateProjectionMatrix()
{
    mProjectionMatrix = Eigen::MatrixXd::Identity(mG.rows(), mG.rows()) - mG * (mG.transpose() * mG).inverse() * mG.transpose();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFETI::CalculateG()
{

    std::vector<int> recvCount;
    std::vector<int> displs;

    boost::mpi::communicator world;

    // recvCount:
    // Contais the number of elements that are received from each process.
    recvCount.clear();
    recvCount.resize(mNumProcesses, 0);

    boost::mpi::all_gather<int>(world,mInterfaceRigidBodyModes.size(),recvCount);

    // displs:
    // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i.
    displs.clear();
    displs.resize(mNumProcesses, 0);
    for (int i = 1; i < mNumProcesses; ++i)
        displs[i] = displs[i-1] + recvCount[i-1];


    const int numInterfaceEqs                       = mInterfaceRigidBodyModes.rows();

    mG.setZero(numInterfaceEqs,mNumRigidBodyModesTotal);

    MPI_Allgatherv(mInterfaceRigidBodyModes.data(),
                   mInterfaceRigidBodyModes.size(),
                   MPI_DOUBLE,
                   mG.data(),
                   recvCount.data(),
                   displs.data(),
                   MPI_DOUBLE,
                   MPI_COMM_WORLD);

}
