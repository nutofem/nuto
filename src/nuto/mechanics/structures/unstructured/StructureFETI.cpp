#include <mpi/mpi.h>
#include <jsoncpp/json/json.h>
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

using std::cout;
using std::endl;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;
using NuTo::Constitutive::eConstitutiveType;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Node::eDof;
using NuTo::Interpolation::eTypeOrder;
using NuTo::ElementData::eElementDataType;
using NuTo::IpData::eIpDataType;
using NuTo::Interpolation::eShapeType;
using NuTo::eGroupId;
using NuTo::eVisualizeWhat;

NuTo::StructureFETI::StructureFETI(int rDimension):
    Structure(rDimension),
    mRank(MPI::COMM_WORLD.Get_rank()),
    mNumProcesses(MPI::COMM_WORLD.Get_size())
{

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

void NuTo::StructureFETI::AssembleBoundaryDofIds()
{

    const auto& dofTypes = GetDofStatus().GetDofTypes();

    int numActiveDofs = 0;

    for (const auto& dofType : dofTypes)
        numActiveDofs           += GetNumActiveDofs(dofType);

    mBoundaryDofIds.setZero(numActiveDofs);

    int offset = 0;
    for (const auto& dofType : dofTypes)
    {
        for (const auto& nodeId : mSubdomainBoundaryNodeIds)
        {
            const std::vector<int> dofIds = NodeGetDofIds(nodeId, dofType);

            for (const auto& dofId : dofIds)
            {
                if(dofId < GetNumActiveDofs(dofType))
                    mBoundaryDofIds.diagonal()(dofId + offset) = 1;
            }

        }
        offset += GetNumActiveDofs(dofType);
    }



}

void NuTo::StructureFETI::AssembleConnectivityMatrix()
{

    const auto& dofTypes = GetDofStatus().GetDofTypes();

    int numActiveDofs           = 0;
    int numLagrangeMultipliers  = 0;

    // this should be read from the mesh file
    //    const int numInterfaceNodesTotal    = 42;

    for (const auto& dofType : dofTypes)
    {
        numActiveDofs           += GetNumActiveDofs(dofType);
        numLagrangeMultipliers  += GetDofDimension(dofType) * mNumInterfaceNodesTotal;
    }

    mConnectivityMatrix.resize(numLagrangeMultipliers, numActiveDofs);

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
                    if (dofVector[i] < GetNumActiveDofs(dofType))
                        mConnectivityMatrix.insert(globalIndex + i + offsetRows , dofVector[i] + offsetCols) = interface.mValue;

            }

        offsetRows += GetDofDimension(dofType) * mNumInterfaceNodesTotal;
        offsetCols += GetNumActiveDofs(dofType);
    }

    if (mRank == 0)
    {
        //        std::cout << mConnectivityMatrix << std::endl;
        std::cout << mConnectivityMatrix.rows() << std::endl;
        std::cout << mConnectivityMatrix.cols() << std::endl;
    }
}

void NuTo::StructureFETI::ImportMesh(std::string rFileName, const int interpolationTypeId)
{

    std::ifstream file(rFileName.c_str(), std::ios::in);

    if (not file.is_open())
        throw MechanicsException(__PRETTY_FUNCTION__, "Mesh file did not open. Check path.");

    mNodes         = ReadNodeData              (file);
    mElements      = ReadElementData           (file);
    mBoundaries    = ReadBoundaryData          (file);
    mInterfaces    = ReadInterfaceData         (file);

    file.close();

    for (const auto& node : mNodes)
        NodeCreate(node.mId, node.mCoordinates.head(2));

    for (const auto& element : mElements)
        ElementCreate(element.mId,interpolationTypeId, element.mNodeIds,eElementDataType::CONSTITUTIVELAWIP,eIpDataType::STATICDATA);

    ElementTotalConvertToInterpolationType();

    NodeBuildGlobalDofs();

}


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
            else if (elementType == 3)
            {

                // restricted to linear quads
                mElements[i].mNodeIds.resize(4);

                mElements[i].mNodeIds[0] = elements["NodalConnectivity"][i][0].asInt();
                mElements[i].mNodeIds[1] = elements["NodalConnectivity"][i][1].asInt();
                mElements[i].mNodeIds[2] = elements["NodalConnectivity"][i][2].asInt();
                mElements[i].mNodeIds[3] = elements["NodalConnectivity"][i][3].asInt();
                mElements[i].mId         = elements["Indices"][i].asInt();
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
    {
        //        std::cout << node.mId << std::endl;
        NodeCreate(node.mId, node.mCoordinates.head(2));
    }

    //    NodeInfo(10);
    //std::cout << "element ids"<< std::endl;

    for (const auto& element : mElements)
    {
        //        std::cout << element.mId << std::endl;
        //        std::cout << element.mNodeIds[0] << element.mNodeIds[1] << element.mNodeIds[2] << element.mNodeIds[3] << std::endl;
        ElementCreate(element.mId,interpolationTypeId, element.mNodeIds,eElementDataType::CONSTITUTIVELAWIP,eIpDataType::STATICDATA);
    }


    ElementTotalConvertToInterpolationType();

    NodeBuildGlobalDofs();

}
