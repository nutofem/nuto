#include "MeshGmsh.h"

#include "base/Exception.h"

#include <array>
#include <fstream>
#include <iostream>
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
    return header;
}


NuTo::MeshGmsh::MeshGmsh(const std::string& fileName)
{
    ReadGmshFile(fileName);
}


std::vector<GmshNode> NuTo::MeshGmsh::ReadNodesASCII(std::ifstream& file)
{
    std::string line;
    do
    {
        std::getline(file, line);
        if (file.eof())
            throw NuTo::Exception(__PRETTY_FUNCTION__, "End of file! - Did not find $Nodes");
    } while (line.compare("$Nodes") != 0);

    std::getline(file, line);
    mNumNodes = std::stoi(line);

    std::vector<GmshNode> nodes(mNumNodes);
    for (GmshNode& node : nodes)
    {
        file >> node.id;
        file >> node.coordinates[0];
        file >> node.coordinates[1];
        file >> node.coordinates[2];
        std::getline(file, line);
    }
    std::getline(file, line);
    if (line.compare("$EndNodes") != 0)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "$EndNodes not found!");

    return nodes;
}


std::vector<GmshElement> NuTo::MeshGmsh::ReadElementsASCII(std::ifstream& file)
{
    const std::array<unsigned int, 32> elementNumNodesLookUp{0,  2,  3,  4,  4, 8, 6,  5,  3,  6, 9,
                                                             10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10,
                                                             12, 15, 15, 21, 4, 5, 6,  20, 26, 56};

    std::string line;
    do
    {
        std::getline(file, line);
        if (file.eof())
            throw NuTo::Exception(__PRETTY_FUNCTION__, "End of file! - Did not find $Elements");
    } while (line.compare("$Elements") != 0);

    std::getline(file, line);
    mNumElements = std::stoi(line);

    std::vector<GmshElement> elements(mNumElements);
    for (GmshElement& element : elements)
    {
        int numTags;
        file >> element.id;
        file >> element.type;
        file >> numTags;

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

    return elements;
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


    GmshHeader header = ReadGmshHeader(file);

    if (header.isBinary)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Not implemented yet");
    }
    else
    {
        std::vector<GmshNode> nodes = ReadNodesASCII(file);
        std::vector<GmshElement> elements = ReadElementsASCII(file);
    }


    file.close();
}
