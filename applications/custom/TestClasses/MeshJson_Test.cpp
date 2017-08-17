#include <mpi.h>

#include <Epetra_MpiComm.h>

#include "mechanics/feti/StructureFeti.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "json.hpp"


using NuTo::Interpolation::eShapeType;
using NuTo::Interpolation::eTypeOrder;


void ImportMeshJson(std::string rFileName, const int interpolationTypeId)
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


void main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    Epetra_MpiComm comm(MPI_COMM_WORLD);
    int rank = comm.MyPID();

    int dim = 2;
    NuTo::StructureFeti structure(dim);
//    std::string meshFile = "meshes/jsonMesh_Test.mesh" + std::to_string(rank);
    std::string meshFile = "meshes/jsonMesh_Test_0" + std::to_string(rank) + ".mesh";

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile, interpolationTypeId);




    MPI_Finalize();
}

