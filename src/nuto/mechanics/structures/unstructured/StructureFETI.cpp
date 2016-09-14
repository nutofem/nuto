#include <mpi/mpi.h>
#include "nuto/mechanics/structures/unstructured/StructureFETI.h"
#include "nuto/mechanics/nodes/NodeBase.h"

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"


#include "nuto/math/SparseMatrixCSR.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"

#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/visualize/VisualizeEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/base/ErrorEnum.h"


NuTo::StructureFETI::StructureFETI(int rDimension, std::string rMeshFile):
    Structure(rDimension),
    mRank(MPI::COMM_WORLD.Get_rank()),
    mNumProcesses(MPI::COMM_WORLD.Get_size())
{

ImportMesh(rMeshFile);

}


NuTo::StructureFETI::NodeList NuTo::StructureFETI::ReadNodeData(std::ifstream &file)
{
    FindKeywordInFile(file, "Nodes");

    int num_nodes = 0;
    file >> num_nodes;

    NodeList nodes(num_nodes);
    for (auto& node : nodes)
        file >> node.mId >> node.mCoordinates[0] >> node.mCoordinates[1] >> node.mCoordinates[2];

    return nodes;

}

NuTo::StructureFETI::ElementList NuTo::StructureFETI::ReadElementData(std::ifstream &file)
{

    FindKeywordInFile(file, "Elements");

    int num_elements = 0;
    file >> num_elements;

    ElementList elements(num_elements);
    for (auto& element : elements)
    {
        file >> element.mId;

        // restricted to linear quads for now
        element.mNodeIds.resize(4);

        for(auto& nodeId: element.mNodeIds)
            file >> nodeId;
    }

    return elements;

}

NuTo::StructureFETI::BoundaryList NuTo::StructureFETI::ReadBoundaryData(std::ifstream &file)
{

    FindKeywordInFile(file, "Boundaries");

    int num_boundaries = 0;
    file >> num_boundaries;

    BoundaryList boundaries(num_boundaries);
    for (auto& boundary : boundaries)
    {
        int num_nodes = 0;
        file >> num_nodes;

        for(int i = 0; i < num_nodes; ++i)
        {
            int globalId  = 0;
            int localId   = 0;
            file >> globalId >> localId;
            boundary.mNodeIdsMap.emplace(globalId, localId);
        }

    }

    return boundaries;

}

NuTo::StructureFETI::InterfaceList NuTo::StructureFETI::ReadInterfaceData(std::ifstream &file)
{
    FindKeywordInFile(file, "Interfaces");

    int num_interfaces = 0;
    file >> num_interfaces;

    InterfaceList interfaces(num_interfaces);
    for (auto& interface : interfaces)
    {
        int num_nodes = 0;
        file >> num_nodes;

        for(int i = 0; i < num_nodes; ++i)
        {
            int globalId  = 0;
            int localId   = 0;
            file >> globalId >> localId;
            interface.mNodeIdsMap.emplace(globalId, localId);
        }

        file >> interface.mValue;


    }

    return interfaces;

}

void NuTo::StructureFETI::FindKeywordInFile(std::ifstream &file, std::string keyword)
{
    std::string line = "initialization";
    while(std::getline(file, line) and line.compare(keyword)!=0);

    if (line.compare(keyword)!=0)
        throw MechanicsException(__PRETTY_FUNCTION__, "Keyword in mesh file not found");
}


void NuTo::StructureFETI::AssembleConnectivityMatrix()
{

    const int num_interface_nodes_global    = 42;
    const int num_boundary_nodes_global     = 42;
    const int num_lagrange_multipliers      = mDimension * (num_interface_nodes_global + num_boundary_nodes_global);
    const int numActiveDofs  = GetNumActiveDofs(NuTo::Node::eDof::DISPLACEMENTS);

    std::cout << "*************************************************** 0 rank "<< numActiveDofs << " " << mRank << std::endl;
    std::cout << "*************************************************** 0 rank "<< num_lagrange_multipliers << " " << mRank << std::endl;
    mConnectivityMatrix.resize(num_lagrange_multipliers, numActiveDofs);

    for (const auto& interface : mInterfaces)
        for (const auto& nodeId : interface.mNodeIdsMap)
        {
            int globalIndex = mDimension * nodeId.first;
            NuTo::FullVector<int, Eigen::Dynamic> displacementDofs;
            NodeGetDisplacementDofs(nodeId.second, displacementDofs);




            if (displacementDofs[0] < numActiveDofs)
                mConnectivityMatrix.insert(globalIndex   , displacementDofs[0]) = interface.mValue;
            if (displacementDofs[1] < numActiveDofs)
                mConnectivityMatrix.insert(globalIndex +1, displacementDofs[1]) = interface.mValue;
        }





    for (const auto& boundary : mBoundaries)
        for (const auto& nodeId : boundary.mNodeIdsMap)
        {
            int globalIndex = mDimension * (nodeId.first + num_interface_nodes_global);
            NuTo::FullVector<int, Eigen::Dynamic> displacementDofs;
            NodeGetDisplacementDofs(nodeId.second, displacementDofs);



            if (displacementDofs[0] < numActiveDofs)
                mConnectivityMatrix.insert(globalIndex   , displacementDofs[0]) = 1.0;

            if (displacementDofs[1] < numActiveDofs)
                mConnectivityMatrix.insert(globalIndex +1, displacementDofs[1]) = 1.0;
        }



}

void NuTo::StructureFETI::AssembleRigidBodyModes()
{

    const int num_interface_nodes_global    = 42;
    const int num_boundary_nodes_global     = 42;
    const int numActiveDofs  = GetNumActiveDofs(NuTo::Node::eDof::DISPLACEMENTS);

    mRigidBodyModes.resize(numActiveDofs,3);
    for (const auto& node : NodeGetNodeMap())
    {
        NuTo::FullVector<int, Eigen::Dynamic> displacementDofs;
        NodeGetDisplacementDofs(node.first, displacementDofs);
        const Eigen::Matrix<double, 2, 1> coordinates = node.second->Get(NuTo::Node::eDof::COORDINATES);

        if (displacementDofs[0] < numActiveDofs)
            mRigidBodyModes.block(displacementDofs[0], 0, 1, 3) << 1., 0., -coordinates.at(1, 0);

        if (displacementDofs[1] < numActiveDofs)
            mRigidBodyModes.block(displacementDofs[1], 0, 1, 3) << 0., 1.,  coordinates.at(0, 0);
    }

    mInterfaceRigidBodyModes = mConnectivityMatrix * mRigidBodyModes;



}

Eigen::MatrixXd NuTo::StructureFETI::GetRigidBodyModes()
{

    return mRigidBodyModes;

}

Eigen::MatrixXd NuTo::StructureFETI::GetInterfaceRigidBodyModes()
{

    return mInterfaceRigidBodyModes;

}


Eigen::SparseMatrix<double>& NuTo::StructureFETI::AssembleStiffnessMatrix()
{
    // assemble stiffness matrix
    NuTo::StructureOutputBlockMatrix stiffnessMatrix = BuildGlobalHessian0();
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrix.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS));

    std::cout << "stiffnessMatrixCSR.GetNumEntries()" << stiffnessMatrixCSR.GetNumEntries() << std::endl;
    std::cout << "stiffnessMatrixCSR.GetNumColumns()" << stiffnessMatrixCSR.GetColumns().size() << std::endl;
    std::cout << "stiffnessMatrixCSR.GetNumRows()" << stiffnessMatrixCSR.GetRowIndex().size() << std::endl;

    std::vector<Eigen::Triplet<double>> tripletList;

    std::vector<double> val = stiffnessMatrixCSR.GetValues();
    std::vector<int> colInd = stiffnessMatrixCSR.GetColumns();
    std::vector<int> rowInd = stiffnessMatrixCSR.GetRowIndex();

    for (unsigned i = 0; i < rowInd.size() - 1; ++i)
    {
        for (int k = rowInd[i]; k < rowInd[i + 1]; ++k)
            tripletList.push_back(Eigen::Triplet<double>(i, colInd[k], val[k]));
    }


    Eigen::SparseMatrix<double> stiffnessMatrixSparse;
    stiffnessMatrixSparse.setFromTriplets(tripletList.begin(), tripletList.end());
    stiffnessMatrixSparse.makeCompressed();
    return stiffnessMatrixSparse;

}

void NuTo::StructureFETI::ImportMesh(std::string rFileName)
{

    std::ifstream file(rFileName.c_str(), std::ios::in);
    if (not file.is_open())
        std::cout << "Mesh file did not open. Check path." << std::endl;


    mNodes         = ReadNodeData              (file);
    mElements      = ReadElementData           (file);
    mBoundaries    = ReadBoundaryData          (file);
    mInterfaces    = ReadInterfaceData         (file);

    file.close();


    for (const auto& node : mNodes)
        NodeCreate(node.mId, node.mCoordinates.head(2));

    int interpolationType = InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    for (const auto& element : mElements)
        ElementCreate(element.mId,interpolationType, element.mNodeIds,NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,NuTo::IpData::eIpDataType::STATICDATA);

    ElementTotalConvertToInterpolationType();

    NodeBuildGlobalDofs();

    const int num_interface_nodes_global    = 42;
    const int num_boundary_nodes_global     = 42;
    const int num_lagrange_multipliers      = mDimension * (num_interface_nodes_global + num_boundary_nodes_global);
    const int numActiveDofs  = GetNumActiveDofs(NuTo::Node::eDof::DISPLACEMENTS);

//    mConnectivityMatrix.resize(num_lagrange_multipliers, numActiveDofs);

//    for (const auto& interface : interfaces)
//        for (const auto& nodeId : interface.mNodeIdsMap)
//        {
//            int globalIndex = mDimension * nodeId.first;
//            NuTo::FullVector<int, Eigen::Dynamic> displacementDofs;
//            NodeGetDisplacementDofs(nodeId.second, displacementDofs);
//            mConnectivityMatrix.insert(globalIndex   , displacementDofs[0]) = interface.mValue;
//            mConnectivityMatrix.insert(globalIndex +1, displacementDofs[1]) = interface.mValue;
//        }

//    for (const auto& boundary : boundaries)
//        for (const auto& nodeId : boundary.mNodeIdsMap)
//        {
//            int globalIndex = mDimension * (nodeId.first + num_interface_nodes_global);
//            NuTo::FullVector<int, Eigen::Dynamic> displacementDofs;
//            NodeGetDisplacementDofs(nodeId.second, displacementDofs);
//            mConnectivityMatrix.insert(globalIndex   , displacementDofs[0]) = 1.0;
//            mConnectivityMatrix.insert(globalIndex +1, displacementDofs[1]) = 1.0;
//        }

//    mRigidBodyModes.resize(numActiveDofs,3);
//    for (const auto& node : NodeGetNodeMap())
//    {
//        NuTo::FullVector<int, Eigen::Dynamic> displacementDofs;
//        NodeGetDisplacementDofs(node.first, displacementDofs);
//        const Eigen::Matrix<double, 2, 1> coordinates = node.second->Get(NuTo::Node::eDof::COORDINATES);
//        mRigidBodyModes.block(displacementDofs[0], 0, 1, 3) << 1., 0., -coordinates.at(1, 0);
//        mRigidBodyModes.block(displacementDofs[1], 0, 1, 3) << 0., 1.,  coordinates.at(0, 0);
//    }

//    mInterfaceRigidBodyModes = mConnectivityMatrix * mRigidBodyModes;



}
