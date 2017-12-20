#include "MeshGmsh.h"

#include "base/Exception.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"

#include <array>
#include <climits>
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <vector>


struct GmshHeader
{
    double version;
    bool isBinary;
    int doubleSize;
};

struct GmshNode
{
    int id = 0;
    std::array<double, 3> coordinates;
};

struct GmshElement
{
    int id = 0;
    int type = 0;
    std::vector<int> tags;
    std::vector<int> nodes;
};

struct GmshPhysicalNames
{
    int id = 0;
    int dimension = 0;
    std::string name;
};


struct GmshFileContent
{
    GmshHeader header;
    int minNodeId = INT_MAX;
    int maxNodeId = INT_MIN;
    int minElementId = INT_MAX;
    int maxElementId = INT_MIN;
    int dimension = 1;
    std::vector<GmshNode> nodes;
    std::vector<GmshElement> elements;
    std::vector<GmshPhysicalNames> physicalNames;
};

GmshHeader ReadGmshHeader(std::ifstream& file)
{
    std::string line;
    getline(file, line);

    int binary;
    GmshHeader header;
    file >> header.version;
    file >> binary;
    file >> header.doubleSize;

    header.isBinary = (binary == 1);

    std::getline(file, line); // endl
    if (header.version < 2.)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Gmsh version 2.0 or higher requiered. - File version is " +
                                                           std::to_string(header.version));
    }
    getline(file, line); // $EndMeshFormat
    return header;
}


NuTo::MeshGmsh::MeshGmsh(const std::string& fileName)
{
    ReadGmshFile(fileName);
}

const NuTo::Group<NuTo::ElementCollection>& NuTo::MeshGmsh::GetPhysicalGroup(std::string physicalName) const
{
    std::transform(physicalName.begin(), physicalName.end(), physicalName.begin(), ::toupper);
    auto physGroupIt = mPhysicalGroups.find(physicalName);
    if (physGroupIt == mPhysicalGroups.end())
        throw Exception(__PRETTY_FUNCTION__, "Couldn't find group with name " + physicalName);
    return physGroupIt->second;
}

std::vector<NuTo::NodeSimple*> NuTo::MeshGmsh::CreateNodes(const GmshFileContent& fileContent)
{
    std::vector<NodeSimple*> nodePtrVec;
    Eigen::VectorXd coords(fileContent.dimension);
    if (fileContent.nodes.size() == static_cast<unsigned int>(fileContent.maxNodeId - fileContent.minNodeId + 1))
    {
        // contiguous IDs
        nodePtrVec.resize(fileContent.nodes.size());
        for (const GmshNode& gmshNode : fileContent.nodes)
        {
            int vectorPos = gmshNode.id - fileContent.minNodeId;
            for (int i = 0; i < fileContent.dimension; ++i)
                coords[i] = gmshNode.coordinates[i];
            nodePtrVec[vectorPos] = &(mMesh.Nodes.Add(coords));
        }
    }
    else
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Non contiguous node IDs - not implemented");
    }
    return nodePtrVec;
}

const NuTo::InterpolationSimple& CreateElementInterpolation(NuTo::MeshFem& mesh, int gmshType, int dimension)
{
    using namespace NuTo;
    switch (gmshType)
    {
    case 1:
        return mesh.CreateInterpolation(InterpolationTrussLinear(dimension));
    case 2:
        return mesh.CreateInterpolation(InterpolationTriangleLinear(dimension));
    case 3:
        return mesh.CreateInterpolation(InterpolationQuadLinear(dimension));
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Unhandled gmsh element type.");
    }
}

std::string GetPhysicalGroupName(const GmshFileContent& fileContent, int groupId)
{
    std::string physicalName{""};
    for (const GmshPhysicalNames& gmshPhysicalName : fileContent.physicalNames)
        if (gmshPhysicalName.id == groupId)
        {
            physicalName = gmshPhysicalName.name;
            break;
        }
    return physicalName;
}


void NuTo::MeshGmsh::CreateElements(const GmshFileContent& fileContent,
                                    const std::vector<NuTo::NodeSimple*>& nodePtrVec)
{
    std::map<int, const InterpolationSimple*> interpolationPtrMap;
    if (fileContent.elements.size() ==
        static_cast<unsigned int>(fileContent.maxElementId - fileContent.minElementId + 1))
    {
        for (GmshElement gmshElement : fileContent.elements)
        {
            // Skip node elements
            if (gmshElement.type == 15)
                continue;


            auto interpolationIter = interpolationPtrMap.find(gmshElement.type);
            if (interpolationIter == interpolationPtrMap.end())
                interpolationIter =
                        interpolationPtrMap
                                .emplace(gmshElement.type,
                                         &CreateElementInterpolation(mMesh, gmshElement.type, fileContent.dimension))
                                .first;
            std::vector<NodeSimple*> elementNodes(gmshElement.nodes.size());
            for (unsigned int i = 0; i < elementNodes.size(); ++i)
            {
                int nodeVectorPos = gmshElement.nodes[i] - fileContent.minNodeId;
                elementNodes[i] = nodePtrVec[nodeVectorPos];
            }


            NuTo::ElementCollection* elementPtr = &(mMesh.Elements.Add({{elementNodes, *(interpolationIter->second)}}));
            AddElementToPhysicalGroup(*elementPtr, GetPhysicalGroupName(fileContent, gmshElement.tags[0]));
            //        gmshElement.type
        }
    }
    else
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Non contiguous element IDs - not implemented");
    }
}


void NuTo::MeshGmsh::AddElementToPhysicalGroup(NuTo::ElementCollection& element, std::string physicalName)
{
    if (physicalName.length() <= 0)
        return;
    std::transform(physicalName.begin(), physicalName.end(), physicalName.begin(), ::toupper);
    // Regarding the map/iterator stuff: https://stackoverflow.com/questions/97050/stdmap-insert-or-stdmap-find
    auto physGroupIt = mPhysicalGroups.lower_bound(physicalName);
    if (physGroupIt != mPhysicalGroups.end() && !(mPhysicalGroups.key_comp()(physicalName, physGroupIt->first)))
        physGroupIt->second.Add(element);
    else
        mPhysicalGroups.emplace_hint(physGroupIt, physicalName, NuTo::Group<ElementCollection>(element));
}

void NuTo::MeshGmsh::CreateMesh(const GmshFileContent& fileContent)
{
    std::vector<NuTo::NodeSimple*> nodePtrVec = CreateNodes(fileContent);
    CreateElements(fileContent, nodePtrVec);
}

void ReadElementsASCII(std::istream& file, GmshFileContent& fileContent)
{
    const std::array<unsigned int, 32> elementNumNodesLookUp{0,  2,  3,  4,  4, 8, 6,  5,  3,  6, 9,
                                                             10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10,
                                                             12, 15, 15, 21, 4, 5, 6,  20, 26, 56};

    std::string line;
    std::getline(file, line);
    int numElements = std::stoi(line);

    fileContent.elements.resize(numElements);
    for (GmshElement& element : fileContent.elements)
    {
        int numTags;
        file >> element.id;
        file >> element.type;
        file >> numTags;

        fileContent.minElementId = std::min(fileContent.minElementId, element.id);
        fileContent.maxElementId = std::max(fileContent.maxElementId, element.id);

        element.tags.resize(numTags);
        for (int& tag : element.tags)
            file >> tag;

        element.nodes.resize(elementNumNodesLookUp[element.type]);
        for (int& node : element.nodes)
            file >> node;
        std::getline(file, line); // endl
    }
    std::getline(file, line);
    if (line.compare("$EndElements") != 0)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndElements not found!");
}


void ReadNodesASCII(std::istream& file, GmshFileContent& fileContent)
{
    std::string line;
    std::getline(file, line);
    int numNodes = std::stoi(line);

    fileContent.nodes.resize(numNodes);
    bool is2d = false;
    bool is3d = false;
    for (GmshNode& node : fileContent.nodes)
    {
        file >> node.id;
        file >> node.coordinates[0];
        file >> node.coordinates[1];
        file >> node.coordinates[2];
        fileContent.minNodeId = std::min(fileContent.minNodeId, node.id);
        fileContent.maxNodeId = std::max(fileContent.maxNodeId, node.id);
        is2d = (is2d || node.coordinates[1] != 0);
        is3d = (is3d || node.coordinates[2] != 0);
        std::getline(file, line);
    }
    if (is2d)
        fileContent.dimension = 2;
    if (is3d)
        fileContent.dimension = 3;

    std::getline(file, line);
    if (line.compare("$EndNodes") != 0)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndNodes not found!");
}

void ReadPhysicalNamesASCII(std::istream& file, GmshFileContent& fileContent)
{
    std::string line;
    std::getline(file, line);
    int numNames = std::stoi(line);

    fileContent.physicalNames.resize(numNames);
    for (GmshPhysicalNames& name : fileContent.physicalNames)
    {
        file >> name.dimension;
        file >> name.id;
        file >> name.name;

        name.name.erase(std::remove(name.name.begin(), name.name.end(), '\"'), name.name.end());
        std::getline(file, line);
    }
    std::getline(file, line);
    if (line.compare("$EndPhysicalNames") != 0)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndPhysicalNames not found!");
}

void ProcessSectionASCII(std::istream& file, GmshFileContent& fileContent)
{
    std::map<std::string, std::function<void(std::istream&, GmshFileContent&)>> sectionEvalutionMap{
            {"NODES", ReadNodesASCII}, {"PHYSICALNAMES", ReadPhysicalNamesASCII}, {"ELEMENTS", ReadElementsASCII},
            //{"PERIODIC", [](std::istream&, GmshFileContent&) { throw NuTo::Exception("Not Implemented"); }},
            //{"NODEDATA", [](std::istream&, GmshFileContent&) { throw NuTo::Exception("Not Implemented"); }},
            //{"ENDNODEDATA", [](std::istream&, GmshFileContent&) { throw NuTo::Exception("Not Implemented"); }},
            //{"ELEMENTNODEDATA", [](std::istream&, GmshFileContent&) { throw NuTo::Exception("Not Implemented"); }},
            //{"INTERPOLATIONSCHEME", [](std::istream&, GmshFileContent&) { throw NuTo::Exception("Not Implemented"); }}
    };

    std::string line;
    std::getline(file, line);
    if (line[0] == '$')
    {
        // remove $
        std::string sectionName = line.substr(1, line.size() - 1);
        // capitalize
        std::transform(sectionName.begin(), sectionName.end(), sectionName.begin(), ::toupper);
        if (sectionName.find("END") == 0)
            throw NuTo::Exception(__PRETTY_FUNCTION__, "Found unprocessed section closure: " + sectionName);

        auto evalIterator = sectionEvalutionMap.find(sectionName);
        if (evalIterator == sectionEvalutionMap.end())
            throw NuTo::Exception(__PRETTY_FUNCTION__, "Unhandled gmsh section type: " + sectionName);
        evalIterator->second(file, fileContent);
    }
}


void NuTo::MeshGmsh::ReadGmshFile(const std::string& fileName)
{
    std::ifstream file;
    file.open(fileName, std::ios::in);


    if (not file.is_open())
    {
        std::cout << fileName << std::endl;
        throw Exception(__PRETTY_FUNCTION__, "Error opening input file "
                                             "" + fileName +
                                                     ""
                                                     " for read access.");
    }
    GmshFileContent fileContent;


    fileContent.header = ReadGmshHeader(file);

    if (fileContent.header.isBinary)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Not implemented yet");
    }
    else
    {


        while (not file.eof())
        {
            ProcessSectionASCII(file, fileContent);
        }
    }
    file.close();

    CreateMesh(fileContent);
}
