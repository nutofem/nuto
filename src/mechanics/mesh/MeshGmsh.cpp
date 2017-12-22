#include "MeshGmsh.h"

#include "base/Exception.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"

#include <array>
#include <climits>
#include <fstream>
#include <unordered_map>


// Helper structs (cpp only)
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    int dimension = 1;
    std::vector<GmshNode> nodes;
    std::vector<GmshElement> elements;
    std::vector<GmshPhysicalNames> physicalNames;
};

// Helper functions (cpp only)
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GmshHeader ReadGmshHeader(std::ifstream& rFile)
{
    std::string line;
    std::getline(rFile, line);

    int binary;
    GmshHeader header;
    rFile >> header.version;
    rFile >> binary;
    rFile >> header.doubleSize;

    header.isBinary = (binary == 1);

    std::getline(rFile, line); // endl
    if (header.version < 2.)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Gmsh version 2.0 or higher requiered. - File version is " +
                                                           std::to_string(header.version));
    }
    std::getline(rFile, line); // $EndMeshFormat
    return header;
}


std::tuple<std::vector<GmshNode>, int> ReadNodesASCII(std::istream& rFile)
{
    std::string line;
    std::getline(rFile, line);
    int numNodes = std::stoi(line);

    int dimension = 1;
    std::vector<GmshNode> nodes(numNodes);
    bool is2d = false;
    bool is3d = false;
    for (GmshNode& node : nodes)
    {
        rFile >> node.id;
        rFile >> node.coordinates[0];
        rFile >> node.coordinates[1];
        rFile >> node.coordinates[2];
        is2d = (is2d || node.coordinates[1] != 0);
        is3d = (is3d || node.coordinates[2] != 0);
        std::getline(rFile, line);
    }
    if (is2d)
        dimension = 2;
    if (is3d)
        dimension = 3;

    std::getline(rFile, line);
    std::transform(line.begin(), line.end(), line.begin(), ::toupper);
    if (line.compare("$ENDNODES") != 0)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndNodes not found!");
    return {std::move(nodes), dimension};
}

std::vector<GmshElement> ReadElementsASCII(std::istream& rFile)
{
    // first element does not exist. Gmsh uses one based indexing
    const std::array<unsigned int, 32> elementNumNodesLookUp{0,  2,  3,  4,  4, 8, 6,  5,  3,  6, 9,
                                                             10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10,
                                                             12, 15, 15, 21, 4, 5, 6,  20, 26, 56};

    std::string line;
    std::getline(rFile, line);
    int numElements = std::stoi(line);

    std::vector<GmshElement> elements(numElements);
    for (GmshElement& element : elements)
    {
        int numTags;
        rFile >> element.id;
        rFile >> element.type;
        rFile >> numTags;

        element.tags.resize(numTags);
        for (int& tag : element.tags)
            rFile >> tag;

        element.nodes.resize(elementNumNodesLookUp[element.type]);
        for (int& node : element.nodes)
            rFile >> node;
        std::getline(rFile, line); // endl
    }
    std::getline(rFile, line);
    std::transform(line.begin(), line.end(), line.begin(), ::toupper);
    if (line.compare("$ENDELEMENTS") != 0)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndElements not found!");
    return elements;
}


std::vector<GmshPhysicalNames> ReadPhysicalNamesASCII(std::istream& rFile)
{
    std::string line;
    std::getline(rFile, line);
    int numNames = std::stoi(line);

    std::vector<GmshPhysicalNames> physicalNames(numNames);
    for (GmshPhysicalNames& name : physicalNames)
    {
        rFile >> name.dimension;
        rFile >> name.id;
        rFile >> name.name;

        // remove quotation marks from string
        name.name.erase(std::remove(name.name.begin(), name.name.end(), '\"'), name.name.end());
        std::getline(rFile, line);
    }
    std::getline(rFile, line);

    std::transform(line.begin(), line.end(), line.begin(), ::toupper);
    if (line.compare("$ENDPHYSICALNAMES") != 0)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndPhysicalNames not found!");
    return physicalNames;
}

void ProcessSectionASCII(std::istream& rFile, GmshFileContent& rFileContent)
{
    std::string line;
    std::getline(rFile, line);
    if (line[0] == '$')
    {
        // remove $
        std::string sectionName = line.substr(1, line.size() - 1);
        // capitalize
        std::transform(sectionName.begin(), sectionName.end(), sectionName.begin(), ::toupper);

        // Execute section read function
        if (sectionName.compare("NODES") == 0)
            std::tie(rFileContent.nodes, rFileContent.dimension) = ReadNodesASCII(rFile);
        else if (sectionName.compare("ELEMENTS") == 0)
            rFileContent.elements = ReadElementsASCII(rFile);
        else if (sectionName.compare("PHYSICALNAMES") == 0)
            rFileContent.physicalNames = ReadPhysicalNamesASCII(rFile);
        else
            throw NuTo::Exception(__PRETTY_FUNCTION__, "Unhandled gmsh section type: " + sectionName);
    }
}


const NuTo::InterpolationSimple& CreateElementInterpolation(NuTo::MeshFem& rMesh, int gmshType, int dimension)
{
    using namespace NuTo;
    switch (gmshType)
    {
    case 1:
        return rMesh.CreateInterpolation(InterpolationTrussLinear());
    case 2:
        return rMesh.CreateInterpolation(InterpolationTriangleLinear());
    case 3:
        return rMesh.CreateInterpolation(InterpolationQuadLinear());
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Unhandled gmsh element type.");
    }
}


std::string GetPhysicalGroupName(const GmshFileContent& fileContent, int groupId)
{
    for (const GmshPhysicalNames& gmshPhysicalName : fileContent.physicalNames)
        if (gmshPhysicalName.id == groupId)
            return gmshPhysicalName.name;
    return "";
}


// Member functions
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NuTo::MeshGmsh::MeshGmsh(const std::string& fileName)
{
    ReadGmshFile(fileName);
}

const NuTo::Group<NuTo::ElementCollection>& NuTo::MeshGmsh::GetPhysicalGroup(std::string physicalName) const
{
    std::transform(physicalName.begin(), physicalName.end(), physicalName.begin(), ::toupper);
    auto physGroupIt = mNamedPhysicalGroups.find(physicalName);
    if (physGroupIt == mNamedPhysicalGroups.end())
        throw Exception(__PRETTY_FUNCTION__, "Couldn't find group with name " + physicalName);
    return *(physGroupIt->second);
}

const NuTo::Group<NuTo::ElementCollection>& NuTo::MeshGmsh::GetPhysicalGroup(int physicalGroupId) const
{
    auto physGroupIt = mPhysicalGroups.find(physicalGroupId);
    if (physGroupIt == mPhysicalGroups.end())
        throw Exception(__PRETTY_FUNCTION__, "Couldn't find group with ID " + std::to_string(physicalGroupId));
    return physGroupIt->second;
}

std::unordered_map<int, NuTo::NodeSimple*> NuTo::MeshGmsh::CreateNodes(const GmshFileContent& fileContent)
{
    std::unordered_map<int, NodeSimple*> nodePtrs;
    Eigen::VectorXd coords(fileContent.dimension);

    for (const GmshNode& gmshNode : fileContent.nodes)
    {
        for (int i = 0; i < fileContent.dimension; ++i)
            coords[i] = gmshNode.coordinates[i];
        nodePtrs[gmshNode.id] = &(mMesh.Nodes.Add(coords));
    }
    return nodePtrs;
}


void NuTo::MeshGmsh::CreateElements(const GmshFileContent& fileContent,
                                    const std::unordered_map<int, NuTo::NodeSimple*>& nodePtrs)
{
    std::map<int, const InterpolationSimple*> interpolationPtrMap;

    for (GmshElement gmshElement : fileContent.elements)
    {
        // Skip node elements
        if (gmshElement.type == 15)
            continue;

        auto interpolationIter = interpolationPtrMap.find(gmshElement.type);
        if (interpolationIter == interpolationPtrMap.end())
            interpolationIter = interpolationPtrMap
                                        .emplace(gmshElement.type, &CreateElementInterpolation(mMesh, gmshElement.type,
                                                                                               fileContent.dimension))
                                        .first;
        std::vector<NodeSimple*> elementNodes(gmshElement.nodes.size());
        for (unsigned int i = 0; i < elementNodes.size(); ++i)
            elementNodes[i] = nodePtrs.at(gmshElement.nodes[i]);

        NuTo::ElementCollection* elementPtr = &(mMesh.Elements.Add({{elementNodes, *(interpolationIter->second)}}));
        AddElementToPhysicalGroup(fileContent, *elementPtr, gmshElement.tags[0]);
    }
}


void NuTo::MeshGmsh::AddElementToPhysicalGroup(const GmshFileContent& fileContent, NuTo::ElementCollection& rElement,
                                               int physicalGroupId)
{
    // Regarding the map/iterator stuff: https://stackoverflow.com/questions/97050/stdmap-insert-or-stdmap-find
    auto physGroupIt = mPhysicalGroups.lower_bound(physicalGroupId);
    if (physGroupIt != mPhysicalGroups.end() && !(mPhysicalGroups.key_comp()(physicalGroupId, physGroupIt->first)))
        // Add to existing group
        physGroupIt->second.Add(rElement);
    else
    {
        // Create new group
        physGroupIt =
                mPhysicalGroups.emplace_hint(physGroupIt, physicalGroupId, NuTo::Group<ElementCollection>(rElement));

        // Create new named group, if a name is defined
        std::string physicalName = GetPhysicalGroupName(fileContent, physicalGroupId);
        if (physicalName.length() > 0)
        {
            std::transform(physicalName.begin(), physicalName.end(), physicalName.begin(), ::toupper);
            mNamedPhysicalGroups.emplace(physicalName, &(physGroupIt->second));
        }
    }
}

void NuTo::MeshGmsh::CreateMesh(const GmshFileContent& fileContent)
{
    CreateElements(fileContent, CreateNodes(fileContent));
}


void NuTo::MeshGmsh::ReadGmshFile(const std::string& fileName)
{
    std::ifstream file;
    file.open(fileName, std::ios::in);

    if (not file.is_open())
        throw Exception(__PRETTY_FUNCTION__, "Error opening input file " + fileName + " for read access.");

    GmshFileContent fileContent;
    fileContent.header = ReadGmshHeader(file);

    if (fileContent.header.isBinary)
        throw Exception(__PRETTY_FUNCTION__, "Not implemented yet");

    while (not file.eof())
        ProcessSectionASCII(file, fileContent);

    file.close();

    CreateMesh(fileContent);
}
