#include "MeshGmsh.h"

#include "nuto/base/Exception.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include "nuto/mechanics/interpolation/InterpolationQuadLinear.h"
#include "nuto/mechanics/interpolation/InterpolationQuadQuadratic.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTetrahedronQuadratic.h"
#include "nuto/mechanics/interpolation/InterpolationBrickLinear.h"
#include "nuto/mechanics/interpolation/InterpolationPrismLinear.h"
#include "nuto/mechanics/interpolation/InterpolationPrismQuadratic.h"
#include "nuto/mechanics/interpolation/InterpolationPyramidLinear.h"

#include <array>
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
    std::string physicalName;
};


struct GmshFileContent
{
    GmshHeader header;
    int dimension = 1;
    std::vector<GmshNode> nodes;
    std::vector<GmshElement> elements;
    std::vector<GmshPhysicalNames> physicalNames;
};

//! Returns true, if sLong and sShort are equal up to the length of sShort.
//! CompareLeft("ABC", "ABC") --> true
//! CompareLeft("ABCD", "ABC") --> true
//! CompareLeft("ABC   ", "ABC") --> true
//! CompareLeft(" ABC", "ABC") --> false
bool CompareLeft(const std::string& sLong, const std::string& sShort)
{
    return sLong.compare(0, sShort.size(), sShort) == 0;
}

void ExpectNextLineToBe(std::ifstream& rFile, std::string expected)
{
    std::string line;
    std::getline(rFile, line);
    if (not CompareLeft(line, expected))
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Expected " + expected + ", got " + line + "!");
}

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
        throw NuTo::Exception(__PRETTY_FUNCTION__,
                              "Gmsh version 2.0 or higher requiered. - File version is " +
                                      std::to_string(header.version));
    }

    if (binary)
    {
        // The next line of the binary format is "1" in binary to control the binary format.
        // http://gmsh.info/doc/texinfo/gmsh.html#MSH-binary-file-format
        int one;
        rFile.read((char*)&one, sizeof(int));
        if (one != 1)
            throw NuTo::Exception(__PRETTY_FUNCTION__, "Invalid binary format.");
        // read the rest of the line
        getline(rFile, line);
    }

    ExpectNextLineToBe(rFile, "$EndMeshFormat");
    return header;
}

int FindDimension(const std::vector<GmshNode>& nodes)
{
    bool is2d = false;
    bool is3d = false;
    for (const GmshNode& node : nodes)
    {
        is2d = (is2d || node.coordinates[1] != 0);
        is3d = (is3d || node.coordinates[2] != 0);
    }
    if (is3d)
        return 3;
    if (is2d)
        return 2;
    return 1;
}

constexpr int NumNodes(int gmshElementType)
{
    // first element does not exist. Gmsh uses one based indexing
    constexpr std::array<unsigned int, 32> elementNumNodesLookUp{0,  2,  3,  4,  4, 8, 6,  5,  3,  6, 9,
                                                                 10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10,
                                                                 12, 15, 15, 21, 4, 5, 6,  20, 26, 56};
    return elementNumNodesLookUp[gmshElementType];
}

std::tuple<std::vector<GmshNode>, int> ReadNodesASCII(std::ifstream& rFile)
{
    std::string line;
    std::getline(rFile, line);
    int numNodes = std::stoi(line);

    std::vector<GmshNode> nodes(numNodes);
    for (GmshNode& node : nodes)
    {
        rFile >> node.id;
        rFile >> node.coordinates[0];
        rFile >> node.coordinates[1];
        rFile >> node.coordinates[2];
        std::getline(rFile, line);
    }
    ExpectNextLineToBe(rFile, "$EndNodes");
    return std::tuple<std::vector<GmshNode>, int>(std::move(nodes), FindDimension(nodes));
}

std::vector<GmshElement> ReadElementsASCII(std::ifstream& rFile)
{

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

        element.nodes.resize(NumNodes(element.type));
        for (int& node : element.nodes)
            rFile >> node;
        std::getline(rFile, line); // endl
    }
    ExpectNextLineToBe(rFile, "$EndElements");
    return elements;
}

std::tuple<std::vector<GmshNode>, int> ReadNodesBinary(std::ifstream& rFile)
{
    std::string line;
    int numNodes = 0;
    rFile >> numNodes;
    std::getline(rFile, line); // read remaining line ('\n')

    std::vector<GmshNode> nodes(numNodes);

    for (GmshNode& node : nodes)
    {
        rFile.read((char*)&node.id, sizeof(int));
        rFile.read((char*)&node.coordinates, 3 * sizeof(double));
    }
    std::getline(rFile, line); // endl

    ExpectNextLineToBe(rFile, "$EndNodes");
    return std::tuple<std::vector<GmshNode>, int>(std::move(nodes), FindDimension(nodes));
}

std::vector<GmshElement> ReadElementsBinary(std::ifstream& rFile)
{
    std::string line;
    int numElements = 0;
    rFile >> numElements;
    std::getline(rFile, line); // read remaining line ('\n')

    std::vector<GmshElement> elements(numElements);
    // read element data
    int element_type;
    int num_elm_follow;
    int num_tags;

    for (int elemCount = 0; elemCount < numElements; elemCount++)
    {
        /* There is a bit of optimization going on the binary format, better read
         * http://gmsh.info/doc/texinfo/gmsh.html#MSH-binary-file-format
         * if you are interested.
         *
         * TLDR: The element section looks as follows:
         * number-of-elements
         * element-header-binary
         * elements-binary
         * element-header-binary
         * elements-binary
         * ...
         *
         * element-header-binary contains 3 ints
         *  [element_type, num_elm_follow, num_tags]
         *
         * So we have an extra loop (indexH) over all elements with the same header
         */

        rFile.read((char*)&element_type, sizeof(int));
        rFile.read((char*)&num_elm_follow, sizeof(int)); // number of elements with the same header
        rFile.read((char*)&num_tags, sizeof(int));

        int cur_num_elm_nodes = NumNodes(element_type);

        for (int indexH = 0; indexH < num_elm_follow; indexH++)
        {

            // set element type
            elements[elemCount].type = element_type;

            // read element number
            rFile.read((char*)&elements[elemCount].id, sizeof(int));

            elements[elemCount].tags.resize(num_tags);
            elements[elemCount].nodes.resize(cur_num_elm_nodes);

            // read tags
            for (int& tag : elements[elemCount].tags)
                rFile.read((char*)&tag, sizeof(int));

            // read nodes
            for (int& node : elements[elemCount].nodes)
                rFile.read((char*)&node, sizeof(int));

            elemCount += indexH;
        }
    }
    getline(rFile, line); // endl
    ExpectNextLineToBe(rFile, "$EndElements");
    return elements;
}


std::vector<GmshPhysicalNames> ReadPhysicalNames(std::ifstream& rFile)
{
    std::string line;
    std::getline(rFile, line);
    int numNames = std::stoi(line);

    std::vector<GmshPhysicalNames> physicalNames(numNames);
    for (GmshPhysicalNames& physicalName : physicalNames)
    {
        rFile >> physicalName.dimension;
        rFile >> physicalName.id;
        rFile >> physicalName.physicalName;

        // remove quotation marks from string
        physicalName.physicalName.erase(
                std::remove(physicalName.physicalName.begin(), physicalName.physicalName.end(), '\"'),
                physicalName.physicalName.end());
        std::getline(rFile, line);
    }
    ExpectNextLineToBe(rFile, "$EndPhysicalNames");
    return physicalNames;
}

void ProcessSection(std::ifstream& rFile, GmshFileContent& rFileContent)
{
    std::string line;
    std::getline(rFile, line);

    if (line.empty())
        return;

    if (CompareLeft(line, "$PhysicalNames"))
    {
        rFileContent.physicalNames = ReadPhysicalNames(rFile);
        return;
    }

    if (CompareLeft(line, "$Nodes"))
    {
        if (rFileContent.header.isBinary)
            std::tie(rFileContent.nodes, rFileContent.dimension) = ReadNodesBinary(rFile);
        else
            std::tie(rFileContent.nodes, rFileContent.dimension) = ReadNodesASCII(rFile);
        return;
    }

    if (CompareLeft(line, "$Elements"))
    {
        if (rFileContent.header.isBinary)
            rFileContent.elements = ReadElementsBinary(rFile);
        else
            rFileContent.elements = ReadElementsASCII(rFile);
        return;
    }

    throw NuTo::Exception(__PRETTY_FUNCTION__, "Unhandled gmsh section type: " + line);
}


const NuTo::InterpolationSimple& CreateElementInterpolation(NuTo::MeshFem& rMesh, int gmshType)
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
    case 4:
        return rMesh.CreateInterpolation(InterpolationTetrahedronLinear());
    case 5:
        return rMesh.CreateInterpolation(InterpolationBrickLinear());
    case 6:
        return rMesh.CreateInterpolation(InterpolationPrismLinear());
    case 7:
        return rMesh.CreateInterpolation(InterpolationPyramidLinear());
    case 8:
        return rMesh.CreateInterpolation(InterpolationTrussLobatto(2));
    case 9:
        return rMesh.CreateInterpolation(InterpolationTriangleQuadratic());
    case 11:
        return rMesh.CreateInterpolation(InterpolationTetrahedronQuadratic());
    case 13:
        return rMesh.CreateInterpolation(InterpolationPrismQuadratic());
    case 16:
        return rMesh.CreateInterpolation(InterpolationQuadQuadratic());
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Unhandled gmsh element type " + std::to_string(gmshType) + ".");
    }
}

std::vector<NuTo::NodeSimple*> GetElementNodes(const std::unordered_map<int, NuTo::NodeSimple*>& nodePtrs,
                                               const GmshElement& gmshElement)
{
    std::vector<NuTo::NodeSimple*> elementNodes(gmshElement.nodes.size());
    for (unsigned int i = 0; i < elementNodes.size(); ++i)
        elementNodes[i] = nodePtrs.at(gmshElement.nodes[i]);

    switch (gmshElement.type)
    {
    case 8: /* TrussQuadratic */
        std::swap(elementNodes[1], elementNodes[2]);
        break;
    default:
        break;
    }
    return elementNodes;
}


std::string GetPhysicalGroupName(const GmshFileContent& fileContent, int groupId)
{
    for (const GmshPhysicalNames& gmshPhysicalName : fileContent.physicalNames)
        if (gmshPhysicalName.id == groupId)
            return gmshPhysicalName.physicalName;
    return "";
}

// Member functions
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NuTo::MeshGmsh::MeshGmsh(const std::string& fileName)
{
    ReadGmshFile(fileName);
}

const NuTo::Group<NuTo::ElementCollectionFem>& NuTo::MeshGmsh::GetPhysicalGroup(std::string physicalName) const
{
    std::transform(physicalName.begin(), physicalName.end(), physicalName.begin(), ::toupper);
    auto physGroupIt = mNamedPhysicalGroups.find(physicalName);
    if (physGroupIt == mNamedPhysicalGroups.end())
        throw Exception(__PRETTY_FUNCTION__, "Couldn't find group with physicalName " + physicalName);
    return *(physGroupIt->second);
}

const NuTo::Group<NuTo::ElementCollectionFem>& NuTo::MeshGmsh::GetPhysicalGroup(int physicalGroupId) const
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
            interpolationIter =
                    interpolationPtrMap.emplace(gmshElement.type, &CreateElementInterpolation(mMesh, gmshElement.type))
                            .first;
        auto elementNodes = GetElementNodes(nodePtrs, gmshElement);

        NuTo::ElementCollectionFem& element = mMesh.Elements.Add({{elementNodes, *(interpolationIter->second)}});
        AddElementToPhysicalGroup(fileContent, element, gmshElement.tags[0]);
    }
}


void NuTo::MeshGmsh::AddElementToPhysicalGroup(const GmshFileContent& fileContent, NuTo::ElementCollectionFem& rElement,
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
                mPhysicalGroups.emplace_hint(physGroupIt, physicalGroupId, NuTo::Group<ElementCollectionFem>(rElement));

        // Create new named group, if a physicalName is defined
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

    while (not file.eof())
        ProcessSection(file, fileContent);

    file.close();

    CreateMesh(fileContent);
}
