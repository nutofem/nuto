#include <mpi/mpi.h>

#include <boost/mpi.hpp>
#include "json.hpp"


#include "mechanics/feti/StructureFeti.h"
#include "mechanics/nodes/NodeBase.h"

#include "base/Exception.h"


#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"

#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"

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

NuTo::StructureFeti::StructureFeti(int rDimension)
    : Structure(rDimension)
    , mRank(boost::mpi::communicator().rank())
    , mNumProcesses(boost::mpi::communicator().size())
{
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::AssembleConnectivityMatrix()
{
    boost::mpi::communicator world;
    const auto& dofTypes = GetDofStatus().GetDofTypes();
    const int numTotalDofs = GetNumTotalDofs();

    mNumLagrangeMultipliers = 0;

    for (const auto& dofType : dofTypes)
    {
        mNumLagrangeMultipliers += NuTo::Node::GetNumComponents(dofType, GetDimension()) * mNumInterfaceNodesTotal;
    }

    mNumLagrangeMultipliers += mNumTotalBoundaryDofIds;
    mNumLagrangeMultipliers += mNumTotalPrescribedDisplacementDofIds;
    mConnectivityMatrix.resize(mNumLagrangeMultipliers, numTotalDofs);
    mPrescribedDofVector.setZero(mNumLagrangeMultipliers);

    int offsetRows = 0;
    int offsetCols = 0;
    for (const auto& dofType : dofTypes)
    {
        for (const auto& interface : mInterfaces)
            for (const auto& nodePair : interface.mNodeIdsMap)
            {

                const std::vector<int> dofVector = NodeGetDofIds(nodePair.second, dofType);
                const int globalIndex = nodePair.first * dofVector.size();

                for (unsigned i = 0; i < dofVector.size(); ++i)
                {
                    if (dofVector[i] < GetNumActiveDofs(dofType))
                    {
                        const int lagrangeMultiplierId = globalIndex + i + offsetRows;
                        const int localDofId = dofVector[i] + offsetCols;

                        mConnectivityMatrix.insert(lagrangeMultiplierId, localDofId) = interface.mValue;

                        mLagrangeMultipliersGlobalIdToLocalId.emplace(lagrangeMultiplierId, localDofId);
                    }
                    else
                        throw Exception(__PRETTY_FUNCTION__, "All DOFs in the connectivity matrix should be active.");
                }
            }

        offsetRows += NuTo::Node::GetNumComponents(dofType, mDimension) * mNumInterfaceNodesTotal;
        offsetCols += GetNumActiveDofs(dofType);
    }

    int globalIndex = mNumLagrangeMultipliers - mNumTotalBoundaryDofIds - mNumTotalPrescribedDisplacementDofIds +
                      mGlobalStartIndexBoundaryDofIds;
    for (const int& id : mBoundaryDofIds)
    {
        mConnectivityMatrix.insert(globalIndex, id) = 1;
        mLagrangeMultipliersGlobalIdToLocalId.emplace(globalIndex, id);
        ++globalIndex;
    }


    globalIndex = mNumLagrangeMultipliers - mNumTotalPrescribedDisplacementDofIds +
                  mGlobalStartIndexPrescribedDisplacementDofIds;
    for (const auto& idAndValue : mPrescribedDisplacementDofIdToValue)
    {
        mConnectivityMatrix.insert(globalIndex, idAndValue.first) = 1;
        mLagrangeMultipliersGlobalIdToLocalId.emplace(globalIndex, idAndValue.first);
        mPrescribedDofVector[globalIndex] = idAndValue.second;

        ++globalIndex;
    }

    // sum up all prescribed dof vectors
    MPI_Allreduce(MPI_IN_PLACE, mPrescribedDofVector.data(), mPrescribedDofVector.size(), MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    GetLogger() << "Number of Lagrange multipliers:              \t" << mNumLagrangeMultipliers << "\n\n";
    GetLogger() << "Total number of interface nodes:             \t" << mNumInterfaceNodesTotal << "\n\n";
    GetLogger() << "Total number of boundary dofs:               \t" << mNumTotalBoundaryDofIds << "\n\n";
    GetLogger() << "Total number of prescribed displacement dofs:\t" << mNumTotalPrescribedDisplacementDofIds << "\n\n";
    GetLogger() << "Total number of dofs:                        \t" << GetNumTotalDofs() << "\n\n";
    GetLogger() << "Total number of active dofs:                 \t" << GetNumTotalActiveDofs() << "\n\n";
    GetLogger() << "Total number of dependent dofs:              \t" << GetNumTotalDependentDofs() << "\n\n";
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::ApplyVirtualConstraints(const std::vector<int>& nodeIdsBoundaries,
                                                  const std::vector<int>& nodeIdsLoads)
{

    std::set<int> setNodeIdsBoundaries(nodeIdsBoundaries.begin(), nodeIdsBoundaries.end());
    std::set<int> setNodeIdsLoads(nodeIdsLoads.begin(), nodeIdsLoads.end());

    std::set<int> setNodeIdsInterfaces;
    for (const auto& interface : mInterfaces)
        for (const auto& nodePair : interface.mNodeIdsMap)
            setNodeIdsInterfaces.insert(nodePair.second);

    std::vector<int> nodeIdsVirtualConstraints;
    for (const auto& nodePair : NodeGetNodeMap())
    {
        if (setNodeIdsBoundaries.find(nodePair.first) == setNodeIdsBoundaries.end() and
            setNodeIdsLoads.find(nodePair.first) == setNodeIdsLoads.end() and
            setNodeIdsInterfaces.find(nodePair.first) == setNodeIdsInterfaces.end())
        {
            nodeIdsVirtualConstraints.push_back(nodePair.first);
        }
    }

    switch (GetDimension())
    {
    case 2:
    {
        auto& firstNode = *NodeGetNodePtr(nodeIdsVirtualConstraints.front());
        Constraints().Add(eDof::DISPLACEMENTS, Constraint::Component(firstNode, {eDirection::X, eDirection::Y}));
        auto& lastNode = *NodeGetNodePtr(nodeIdsVirtualConstraints.back());
        Constraints().Add(eDof::DISPLACEMENTS, Constraint::Component(lastNode, {eDirection::Y}));

        GetLogger() << "Applied virtual constraint to node id: \t" << nodeIdsVirtualConstraints.front()
                    << "\t in X \n\n";
        GetLogger() << "Applied virtual constraint to node id: \t" << nodeIdsVirtualConstraints.front()
                    << "\t in X \n\n";
        GetLogger() << "Applied virtual constraint to node id: \t" << nodeIdsVirtualConstraints.back()
                    << "\t in Y \n\n";
        break;
    }
    case 3:
    {
        auto& firstNode = *NodeGetNodePtr(nodeIdsVirtualConstraints.front());
        auto& thirdNode = *NodeGetNodePtr(nodeIdsVirtualConstraints[2]);
        auto& lastNode = *NodeGetNodePtr(nodeIdsVirtualConstraints.back());
        Constraints().Add(eDof::DISPLACEMENTS,
                          Constraint::Component(firstNode, {eDirection::X, eDirection::Y, eDirection::Z}));
        Constraints().Add(eDof::DISPLACEMENTS, Constraint::Component(thirdNode, {eDirection::Z}));
        Constraints().Add(eDof::DISPLACEMENTS, Constraint::Component(lastNode, {eDirection::X, eDirection::Y}));

        GetLogger() << "Applied virtual constraint to node id: \t" << nodeIdsVirtualConstraints.front()
                    << "\t in X \n\n";
        GetLogger() << "Applied virtual constraint to node id: \t" << nodeIdsVirtualConstraints.front()
                    << "\t in Y \n\n";
        GetLogger() << "Applied virtual constraint to node id: \t" << nodeIdsVirtualConstraints.front()
                    << "\t in Z \n\n";
        GetLogger() << "Applied virtual constraint to node id: \t" << nodeIdsVirtualConstraints.back()
                    << "\t in X \n\n";
        GetLogger() << "Applied virtual constraint to node id: \t" << nodeIdsVirtualConstraints.back()
                    << "\t in Y \n\n";
        GetLogger() << "Applied virtual constraint to node id: \t" << nodeIdsVirtualConstraints[2] << "\t in Z \n\n";

        break;
    }
    default:
        throw Exception(__PRETTY_FUNCTION__, "Not implemented for dimension: " + std::to_string(GetDimension()));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::ImportMeshJson(std::string rFileName, const int interpolationTypeId)
{
    nlohmann::json root;

    std::ifstream file(rFileName.c_str(), std::ios::in);

    file >> root;

    // only supports nodes.size() == 1
    for (auto const& nodes : root["Nodes"])
    {
        mNodes.resize(nodes["Coordinates"].size());
        for (unsigned i = 0; i < mNodes.size(); ++i)
        {
            mNodes[i].mCoordinates[0] = nodes["Coordinates"][i][0];
            mNodes[i].mCoordinates[1] = nodes["Coordinates"][i][1];
            mNodes[i].mCoordinates[2] = nodes["Coordinates"][i][2];
            mNodes[i].mId = nodes["Indices"][i];
        }
    }


    // only supports elements.size() == 1
    for (auto const& elements : root["Elements"])
    {
        mElements.resize(elements["NodalConnectivity"].size());
        const int elementType = elements["Type"];


        for (unsigned i = 0; i < mElements.size(); ++i)
        {
            if (elementType == 1)
            {
                mSubdomainBoundaryNodeIds.insert(elements["NodalConnectivity"][i][0].get<int>());
                mSubdomainBoundaryNodeIds.insert(elements["NodalConnectivity"][i][1].get<int>());
            }
            else if (elementType == 2) // 3 node tri element
            {
                mElements[i].mNodeIds.resize(3);

                mElements[i].mNodeIds[0] = elements["NodalConnectivity"][i][0];
                mElements[i].mNodeIds[1] = elements["NodalConnectivity"][i][1];
                mElements[i].mNodeIds[2] = elements["NodalConnectivity"][i][2];
                mElements[i].mId = elements["Indices"][i];
            }
            else if (elementType == 3) // 4 node quad element
            {

                mElements[i].mNodeIds.resize(4);

                mElements[i].mNodeIds[0] = elements["NodalConnectivity"][i][0];
                mElements[i].mNodeIds[1] = elements["NodalConnectivity"][i][1];
                mElements[i].mNodeIds[2] = elements["NodalConnectivity"][i][2];
                mElements[i].mNodeIds[3] = elements["NodalConnectivity"][i][3];
                mElements[i].mId = elements["Indices"][i];
            }
            else if (elementType == 5) // 8 node hexahedron
            {
                const int numNodes = 8;
                mElements[i].mNodeIds.resize(numNodes);
                for (int iNode = 0; iNode < numNodes; ++iNode)
                    mElements[i].mNodeIds[iNode] = elements["NodalConnectivity"][i][iNode];

                mElements[i].mId = elements["Indices"][i];
            }
            else
            {
                throw Exception(__PRETTY_FUNCTION__, "Import of element type not implemented. Element type id = " +
                                                             std::to_string(elementType));
            }
        }
    }


    mInterfaces.resize(root["Interface"].size());
    for (unsigned i = 0; i < mInterfaces.size(); ++i)
    {

        int globalId = root["Interface"][i]["GlobalStartId"][0];

        mInterfaces[i].mValue = root["Interface"][i]["Value"][0];

        for (unsigned k = 0; k < root["Interface"][i]["NodeIds"][0].size(); ++k)
        {
            mInterfaces[i].mNodeIdsMap.emplace(globalId, root["Interface"][i]["NodeIds"][0][k]);
            mSubdomainBoundaryNodeIds.insert(root["Interface"][i]["NodeIds"][0][k].get<int>());
            globalId++;
        }
    }


    mNumInterfaceNodesTotal = root["NumInterfaceNodes"][0];

    file.close();

    for (const auto& node : mNodes)
        NodeCreate(node.mId, node.mCoordinates.head(GetDimension()));

    for (const auto& element : mElements)
        ElementCreate(element.mId, interpolationTypeId, element.mNodeIds);

    ElementTotalConvertToInterpolationType(1.e-15, 1.);

    NodeBuildGlobalDofs();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::CalculateRigidBodyModesTotalFETI()
{

    const auto displacements = NuTo::Node::eDof::DISPLACEMENTS;
    const int numTotalDofs = GetNumTotalDofs();

    switch (GetDimension())
    {
    case 1:
    {
        mNumRigidBodyModes = 1;
        mRigidBodyModes.setOnes(numTotalDofs, mNumRigidBodyModes);
    }
    break;
    case 2:
    {
        mNumRigidBodyModes = 3;

        mRigidBodyModes.setZero(numTotalDofs, mNumRigidBodyModes);

        for (const auto& nodePair : mNodeMap)
        {
            const std::vector<int> dofIds = NodeGetDofIds(nodePair.first, displacements);

            const Eigen::Matrix<double, 2, 1> coordinates = nodePair.second->Get(NuTo::Node::eDof::COORDINATES);

            if (IsActiveDofId(dofIds[0], displacements))
                mRigidBodyModes.row(dofIds[0]) << 1., 0., -coordinates[1];
            else
            {

                const int rowId =
                        numTotalDofs - GetNumActiveDofs(displacements) - GetNumDependentDofs(displacements) + dofIds[0];
                mRigidBodyModes.row(rowId) << 1., 0., -coordinates[1];
            }


            if (IsActiveDofId(dofIds[1], displacements))
            {
                mRigidBodyModes.row(dofIds[1]) << 0., 1., coordinates[0];
            }
            else
            {
                const int rowId =
                        numTotalDofs - GetNumActiveDofs(displacements) - GetNumDependentDofs(displacements) + dofIds[1];
                mRigidBodyModes.row(rowId) << 0., 1., coordinates[0];
            }
        }
    }
    break;
    case 3:
    {
        mNumRigidBodyModes = 6;

        mRigidBodyModes.setZero(numTotalDofs, mNumRigidBodyModes);

        for (const auto& nodePair : mNodeMap)
        {
            const std::vector<int> dofIds = NodeGetDofIds(nodePair.first, displacements);

            const Eigen::Matrix<double, 3, 1> coordinates = nodePair.second->Get(NuTo::Node::eDof::COORDINATES);

            const double x = coordinates[0];
            const double y = coordinates[1];
            const double z = coordinates[2];

            if (IsActiveDofId(dofIds[0], displacements))
                mRigidBodyModes.row(dofIds[0]) << 1., 0., 0., 0., -z, y;
            else
            {

                const int rowId =
                        numTotalDofs - GetNumActiveDofs(displacements) - GetNumDependentDofs(displacements) + dofIds[0];
                mRigidBodyModes.row(rowId) << 1., 0., 0., 0., -z, y;
            }


            if (IsActiveDofId(dofIds[1], displacements))
            {
                mRigidBodyModes.row(dofIds[1]) << 0., 1., 0., z, 0, -x;
            }
            else
            {
                const int rowId =
                        numTotalDofs - GetNumActiveDofs(displacements) - GetNumDependentDofs(displacements) + dofIds[1];
                mRigidBodyModes.row(rowId) << 0., 1., 0., z, 0, -x;
            }

            if (IsActiveDofId(dofIds[2], displacements))
            {
                mRigidBodyModes.row(dofIds[2]) << 0., 0., 1., -y, x, 0.;
            }
            else
            {
                const int rowId =
                        numTotalDofs - GetNumActiveDofs(displacements) - GetNumDependentDofs(displacements) + dofIds[2];
                mRigidBodyModes.row(rowId) << 0., 0., 1., -y, x, 0.;
            }
        }
    }
    break;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Structural dimension not supported yet.");
    }

    MPI_Allreduce(&mNumRigidBodyModes, &mNumRigidBodyModesTotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    GetLogger() << "Number of rigid body modes:        \t" << mNumRigidBodyModes << "\n\n";
    GetLogger() << "Total number of rigid body modes:  \t" << mNumRigidBodyModesTotal << "\n\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::CreateDummy1D()
{
    mInterfaces.resize(1);

    if (mRank == 0)
    {
        mInterfaces[0].mValue = 1;
        mInterfaces[0].mNodeIdsMap.emplace(0, 3);
    }

    if (mRank == 1)
    {
        mInterfaces[0].mValue = -1;
        mInterfaces[0].mNodeIdsMap.emplace(0, 0);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::CheckRigidBodyModes(const StructureOutputBlockMatrix& hessian0, const double tolerance) const
{

    SparseMatrix hessian0_JJ = hessian0.JJ.ExportToEigenSparseMatrix();
    SparseMatrix hessian0_JK = hessian0.JK.ExportToEigenSparseMatrix();


    Eigen::MatrixXd zeroMatrix0 = hessian0_JJ * mRigidBodyModes.topRows(GetNumTotalActiveDofs());

    zeroMatrix0 += hessian0_JK * mRigidBodyModes.bottomRows(mNumRigidBodyModes);
    const double norm0 = std::max(zeroMatrix0.maxCoeff(), std::abs(zeroMatrix0.minCoeff()));

    MPI_Barrier(MPI_COMM_WORLD);

    SparseMatrix hessian0_KJ = hessian0.KJ.ExportToEigenSparseMatrix();
    SparseMatrix hessian0_KK = hessian0.KK.ExportToEigenSparseMatrix();

    Eigen::MatrixXd zeroMatrix1 = hessian0_KJ * mRigidBodyModes.topRows(GetNumTotalActiveDofs());
    zeroMatrix1 += hessian0_KK * mRigidBodyModes.bottomRows(mNumRigidBodyModes);
    const double norm1 = std::max(zeroMatrix1.maxCoeff(), std::abs(zeroMatrix1.minCoeff()));

    MPI_Barrier(MPI_COMM_WORLD);

    GetLogger() << "Norm of K*R:  \t" << norm0 << "\t and \t" << norm1 << "\t tolerance: \t" << tolerance << "\n\n";
    //                << "Rigid body modes      : \n"       << mRigidBodyModes                << "\n\n";

    MPI_Barrier(MPI_COMM_WORLD);

    assert((norm0 < tolerance) and (norm1 < tolerance) and
           "Calculated rigid body modes are not in the null space of K");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::CheckProjectionMatrix(const double tolerance) const
{

    const Eigen::MatrixXd zeroMatrix = mProjectionMatrix - mProjectionMatrix * mProjectionMatrix;

    const double norm = std::max(zeroMatrix.maxCoeff(), std::abs(zeroMatrix.minCoeff()));
    GetLogger() << "Subdomain: \t" << mRank << "\t norm of P-PP: \t \t" << norm << "\t tolerance: \t" << tolerance
                << "\n\n";

    assert((norm < tolerance) and "Projection matrix does not satisfy: P - PP = 0");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::CheckProjectionOfCoarseGrid(const double tolerance) const
{

    const Eigen::MatrixXd zeroMatrix = mProjectionMatrix.transpose() * mG;

    const double norm = std::max(zeroMatrix.maxCoeff(), std::abs(zeroMatrix.minCoeff()));

    GetLogger() << "Subdomain: \t" << mRank << "\t norm of Ptrans * G: \t" << norm << "\t tolerance: \t" << tolerance
                << "\n\n";

    assert((norm < tolerance) and "Projection of coarse space does not satisfy: Ptrans*G = 0");
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::CalculateInterfaceRigidBodyModes()
{
    mInterfaceRigidBodyModes = mConnectivityMatrix * mRigidBodyModes;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::ApplyConstraintsTotalFeti(const std::vector<int>& dofIds)
{
    boost::mpi::communicator world;

    mBoundaryDofIds = dofIds;

    // determine the global ids for the constraints

    // recvCount:
    // Contains the number of elements that are received from each process.
    std::vector<int> recvCount(mNumProcesses, 0);

    boost::mpi::all_gather<int>(world, mBoundaryDofIds.size(), recvCount);


    // displs:
    // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i.
    std::vector<int> displs;
    displs.resize(mNumProcesses, 0);
    for (int i = 1; i < mNumProcesses; ++i)
        displs[i] = displs[i - 1] + recvCount[i - 1];


    int numLocalBoundaryDofIds = mBoundaryDofIds.size();
    MPI_Allreduce(&numLocalBoundaryDofIds, &mNumTotalBoundaryDofIds, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


    GetLogger() << "Number of boundary dof ids:       \t" << mBoundaryDofIds.size() << "\n \n";
    GetLogger() << "Total number of boundary dof ids: \t" << mNumTotalBoundaryDofIds << "\n \n";


    mGlobalStartIndexBoundaryDofIds = displs[mRank];
}

void NuTo::StructureFeti::ApplyConstraintsTotalFeti(const Group<NodeBase>& nodeGroup)
{
    boost::mpi::communicator world;

    std::vector<int> boundaryNodeIds = nodeGroup.GetMemberIds();

    for (const int nodeId : boundaryNodeIds)
    {
        std::vector<int> dofIds = NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);

        for (const auto& id : dofIds)
            mBoundaryDofIds.push_back(id);
    }


    // determine the global ids for the constraints

    // recvCount:
    // Contains the number of elements that are received from each process.
    std::vector<int> recvCount(mNumProcesses, 0);

    boost::mpi::all_gather<int>(world, mBoundaryDofIds.size(), recvCount);

    // displs:
    // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i.
    std::vector<int> displs;
    displs.resize(mNumProcesses, 0);
    for (int i = 1; i < mNumProcesses; ++i)
        displs[i] = displs[i - 1] + recvCount[i - 1];


    int numLocalBoundaryDofIds = mBoundaryDofIds.size();
    MPI_Allreduce(&numLocalBoundaryDofIds, &mNumTotalBoundaryDofIds, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


    GetLogger() << "Number of boundary dof ids:       \t" << mBoundaryDofIds.size() << "\n \n";
    GetLogger() << "Total number of boundary dof ids: \t" << mNumTotalBoundaryDofIds << "\n \n";


    mGlobalStartIndexBoundaryDofIds = displs[mRank];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::ApplyPrescribedDisplacements(const std::map<int, double> dofIdAndPrescribedDisplacementMap)
{
    boost::mpi::communicator world;

    mPrescribedDisplacementDofIdToValue = dofIdAndPrescribedDisplacementMap;

    // recvCount:
    // Contains the number of elements that are received from each process.
    std::vector<int> recvCount(mNumProcesses, 0);

    int numLocalDofIds = mPrescribedDisplacementDofIdToValue.size();
    boost::mpi::all_gather<int>(world, numLocalDofIds, recvCount);

    // displs:
    // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i.
    std::vector<int> displs;
    displs.resize(mNumProcesses, 0);
    for (int i = 1; i < mNumProcesses; ++i)
        displs[i] = displs[i - 1] + recvCount[i - 1];

    MPI_Allreduce(&numLocalDofIds, &mNumTotalPrescribedDisplacementDofIds, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    GetLogger() << "Number of prescribed displacement dof ids:       \t" << mPrescribedDisplacementDofIdToValue.size()
                << "\n \n";
    GetLogger() << "Total number of prescribed displacement dof ids: \t" << mNumTotalPrescribedDisplacementDofIds
                << "\n \n";


    mGlobalStartIndexPrescribedDisplacementDofIds = displs[mRank];

    GetLogger() << "mGlobalStartIndexPrescribedDisplacementDofIds: \t" << mGlobalStartIndexPrescribedDisplacementDofIds
                << "\n \n";

    for (const auto& idAndValue : mPrescribedDisplacementDofIdToValue)
        GetLogger() << "Prescribed displacement dof ids:       \t" << idAndValue.first << "\n \n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::CalculateProjectionMatrix()
{
    mProjectionMatrix =
            Eigen::MatrixXd::Identity(mG.rows(), mG.rows()) - mG * (mG.transpose() * mG).inverse() * mG.transpose();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NuTo::StructureFeti::CalculateG()
{

    std::vector<int> recvCount;
    std::vector<int> displs;

    boost::mpi::communicator world;

    // recvCount:
    // Contais the number of elements that are received from each process.
    recvCount.clear();
    recvCount.resize(mNumProcesses, 0);

    boost::mpi::all_gather<int>(world, mInterfaceRigidBodyModes.size(), recvCount);

    // displs:
    // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i.
    displs.clear();
    displs.resize(mNumProcesses, 0);
    for (int i = 1; i < mNumProcesses; ++i)
        displs[i] = displs[i - 1] + recvCount[i - 1];


    const int numInterfaceEqs = mInterfaceRigidBodyModes.rows();

    mG.setZero(numInterfaceEqs, mNumRigidBodyModesTotal);

    MPI_Allgatherv(mInterfaceRigidBodyModes.data(), mInterfaceRigidBodyModes.size(), MPI_DOUBLE, mG.data(),
                   recvCount.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
}
