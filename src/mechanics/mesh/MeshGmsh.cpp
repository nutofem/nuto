#include "MeshGmsh.h"

#include "base/Exception.h"

#include <array>
#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
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
    for (GmshNode& node : fileContent.nodes)
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
}
