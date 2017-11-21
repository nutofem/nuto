#include <string>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <algorithm>

#include "mechanics/mesh/MeshCompanion.h"



class MeshFileGenerator
{
private:
    struct GmshHeader
    {
        double version;
        bool isBinary;
        int double_size;
    };

    struct GmshNode
    {
        int id = 0;
        double Coordinates[3];
    };

    struct GmshElementTags
    {
        int size = 0;
        int physGroup = 0;
        int geomGroup = 0;
        int numDomains = 0;
        int masterDomain = 0;
        std::vector<int> slaveDomains;
    };

    struct GmshElement
    {
        int id = 0;
        int type = 0;
//        std::vector<int> tags;
        GmshElementTags tags;
        std::vector<int> nodes;
    };

    GmshHeader ReadGmshHeader(std::ifstream& rFile)
    {
        // ignore first line
        std::string line;
        getline(rFile, line);


        GmshHeader header;
        int binary;
        rFile >> header.version;
        rFile >> binary;
        rFile >> header.double_size;

        header.isBinary = binary == 1;

        std::getline(rFile, line); // endl

        if (static_cast<int>(std::floor(header.version)) != 2)
            throw std::runtime_error("Incompatible version. Version 2.x required.");

        return header;
    }

    std::vector<GmshNode> ReadNodesASCII(std::ifstream& rFile)
    {
        std::string line;
        getline(rFile, line);
        if (line != "$EndMeshFormat")
            throw std::runtime_error("$EndMeshFormat not found.");

        // begin node section
        getline(rFile, line);
        if (line != "$Nodes")
            throw std::runtime_error("$Nodes not found.");

        // read number of nodes
        getline(rFile, line);
        int numNodes = std::stoi(line);

        // read node data
        std::vector<GmshNode> nodes(numNodes);
        for (GmshNode& node : nodes)
        {
            rFile >> node.id;
            --node.id;
            rFile >> node.Coordinates[0];
            rFile >> node.Coordinates[1];
            rFile >> node.Coordinates[2];
            getline(rFile, line); // endl
        }

        // end node section
        getline(rFile, line);
        if (line != "$EndNodes")
            throw std::runtime_error("$EndNodes not found.");

        return nodes;
    }

    int GetNumNodesPerElementType(int rElementType)
    {
        const std::vector<int> numNodesPerElement(
                {0, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13, 0, 10, 0, 15});
        return numNodesPerElement[rElementType];
    }

    std::vector<GmshElement> ReadElementsASCII(std::ifstream& rFile)
    {
        std::string line;
        getline(rFile, line);
        if (line != "$Elements")
            throw std::runtime_error("$Elements not found.");

        // read number of elements
        getline(rFile, line);
        int numElements = std::stoi(line);
        std::vector<GmshElement> elements(numElements);


        // read element data
        for (GmshElement& element : elements)
        {
            int numTags;
            rFile >> element.id;
            --element.id;
            rFile >> element.type;
            rFile >> numTags;

            int numNodesInThisElement = GetNumNodesPerElementType(element.type);

//            element.tags.resize(numTags);
//            for (auto& tag : element.tags)
//                rFile >> tag;

            std::vector<int> tags(numTags);
            for (auto& tag : tags)
                rFile >> tag;

            element.tags.size = numTags;
            element.tags.physGroup = tags[0];
            element.tags.geomGroup = tags[1];
            element.tags.numDomains = tags[2];
            element.tags.masterDomain = tags[3];
            element.tags.slaveDomains.resize(element.tags.numDomains-1);
            for (int i = 0; i < element.tags.numDomains-1; ++i)
                element.tags.slaveDomains[i] = tags[4+i];

            element.nodes.resize(numNodesInThisElement);
            for (auto& node : element.nodes)
            {
                rFile >> node;
                --node;
            }

            getline(rFile, line); // endl;
        }
        // end element section
        getline(rFile, line);
        if (line != "$EndElements")
            throw std::runtime_error("$EndElements not found.");

        return elements;
    }

    std::vector<GmshNode> ReadNodesBinary(std::ifstream& rFile)
    {
        std::string line;

        // read the first two lines again
        getline(rFile, line);
        getline(rFile, line);

        // check size of integer
        int one;
        rFile.read((char*)&one, sizeof(int));
        if (one != 1)
            throw std::runtime_error("Invalid binary format.");
        rFile.seekg(1, std::ios::cur);

        getline(rFile, line);
        if (line != "$EndMeshFormat")
            throw std::runtime_error("$EndMeshFormat not found.");

        // begin node section
        getline(rFile, line);
        if (line != "$Nodes")
        {
            throw std::runtime_error("$Nodes not found.");
        }

        // read number of nodes
        int numNodes;
        rFile >> numNodes;
        getline(rFile, line); // endl

        // read node data
        std::vector<GmshNode> nodes(numNodes);

        for (auto& node : nodes)
        {
            rFile.read((char*)&node.id, sizeof(int));
            --node.id;
            rFile.read((char*)node.Coordinates, 3 * sizeof(double));
        }
        // endl
        getline(rFile, line);

        getline(rFile, line);
        if (line != "$EndNodes")
            throw std::runtime_error("$EndNodes not found.");

        // begin element section
        getline(rFile, line);
        if (line != "$Elements")
            throw std::runtime_error("$Elements not found.");

        return nodes;
    }

    std::vector<GmshElement> ReadElementsBinary(std::ifstream& rFile)
    {
        std::string line;


        // read number of elements
        int numElements;
        rFile >> numElements;
        getline(rFile, line); // endl

        std::vector<GmshElement> elements(numElements);

        // read element data
        int element_type;
        int num_elm_follow;
        int num_tags;

        for (int elemCount = 0; elemCount < numElements; elemCount++)
        {

            // Read element type
            rFile.read((char*)&element_type, sizeof(int));

            // Read num of Elem with the same header
            rFile.read((char*)&num_elm_follow, sizeof(int));

            // set num_elemt_node
            int cur_num_elm_nodes = GetNumNodesPerElementType(element_type);

            // Read numOfTags
            rFile.read((char*)&num_tags, sizeof(int));

            for (int indexH = 0; indexH < num_elm_follow; indexH++)
            {

                // set element type
                elements[elemCount].type = element_type;

                // read element number
                rFile.read((char*)&elements[elemCount].id, sizeof(int));
                --elements[elemCount].id;

//                elements[elemCount].tags.resize(num_tags);
                elements[elemCount].tags.size = num_tags;
                elements[elemCount].nodes.resize(cur_num_elm_nodes);

                // read tags
//                for (int tagCount = 0; tagCount < num_tags; tagCount++)
//                    rFile.read((char*)&elements[elemCount].tags[tagCount], sizeof(int));

                std::vector<int> tags(num_tags);
                for (int tagCount = 0; tagCount < num_tags; tagCount++)
                    rFile.read((char*)&tags[tagCount], sizeof(int));

                elements[elemCount].tags.physGroup = tags[0];
                elements[elemCount].tags.geomGroup = tags[1];
                elements[elemCount].tags.numDomains = tags[2];
                elements[elemCount].tags.masterDomain = tags[3];
                elements[elemCount].tags.slaveDomains.resize(elements[elemCount].tags.numDomains-1);
                for (int i = 0; i < elements[elemCount].tags.numDomains-1; ++i)
                    elements[elemCount].tags.slaveDomains[i] = tags[4+i];

                // read nodes
                for (int nodeCount = 0; nodeCount < cur_num_elm_nodes; nodeCount++)
                {
                    rFile.read((char*)&elements[elemCount].nodes[nodeCount], sizeof(int));
                    --elements[elemCount].nodes[nodeCount];
                }

                elemCount += indexH;
            }
        }
        getline(rFile, line); // endl

        // end element section
        getline(rFile, line);
        if (line != "$EndElements")
            throw std::runtime_error("$EndElements not found.");
        return elements;
    }




public:
    enum NuTo_DofTypes : int
    {
        DISPLACEMENTS,
        TEMPERATURE
    };

    enum NuTo_InterpolationOrders : int
    {
        EQUIDISTANT1,
        EQUIDISTANT2
    };

    MeshFileGenerator(){}

    void generateMeshFromGmsh(std::string rFileName, std::string rOutputFile)
    {
        std::ifstream file_gmsh(rFileName.c_str(), std::ios::in);
        if (not file_gmsh.is_open())
        {
            std::cout << rFileName << std::endl;
            throw std::runtime_error("Error opening input file for read access.");
        }


        GmshHeader header = ReadGmshHeader(file_gmsh);

        std::vector<GmshNode> nodes;
        std::vector<GmshElement> elements;

        if (header.isBinary)
        {
            if (header.double_size != sizeof(double))
                throw std::runtime_error("Invalid size of double.");

            // close rFile and open as binary
            file_gmsh.close();
            file_gmsh.open(rFileName.c_str(), std::ios::in | std::ios::binary);
            if (file_gmsh.is_open() == false)
                throw std::runtime_error("Error opening input rFile for read access.");
            nodes = ReadNodesBinary(file_gmsh);
            elements = ReadElementsBinary(file_gmsh);
        }
        else
        {
            nodes = ReadNodesASCII(file_gmsh);
            elements = ReadElementsASCII(file_gmsh);
        }

        file_gmsh.close();


        const std::string blank = "   ";
        int elementCount = elements.size();
        int nodeCount = nodes.size();

        std::ofstream file_mesh(rOutputFile);

        file_mesh << "{\n";
        file_mesh << blank << "\"Elements\" : [\n";
        file_mesh << blank << blank << "{\n";
        file_mesh << blank << blank << blank << "\"Indices\" : [ ";

        for (int i = 0; i < elementCount - 1; ++i)
        {
            file_mesh << elements[i].id << ", ";
        }
        file_mesh << elements[elementCount-1].id << "],\n";

        file_mesh << blank << blank << blank << "\"NodalConnectivity\" : [\n";
        for (int i = 0; i < elementCount - 1; ++i)
        {
            file_mesh << blank << blank << blank << blank << "[ ";
            for (int j = 0; j < elements[i].nodes.size()-1; ++j)
            {
                file_mesh << elements[i].nodes[j] << ", ";
            }
            file_mesh << elements[i].nodes[elements[i].nodes.size()-1] << " ],\n";
        }
        file_mesh << blank << blank << blank << blank << "[ ";
        for (int j = 0; j < elements[elementCount-1].nodes.size()-1; ++j)
        {
            file_mesh << elements[elementCount-1].nodes[j] << ", ";
        }
        file_mesh << elements[elementCount-1].nodes[elements[elementCount-1].nodes.size()-1] << " ]\n";
        file_mesh << blank << blank << blank << "],\n";

        file_mesh << blank << blank << blank << "\"Type\" : " << elements[0].type << ",\n";
        file_mesh << blank << blank << blank << "\"TypeName\" : " << "\"QUAD2D\"" << "\n";     //TODO: conversion Type -> TypeName

        file_mesh << blank << blank << "}\n";

        file_mesh << blank << "],\n";

        file_mesh << blank << "\"Nodes\" : [\n";
        file_mesh << blank << blank << "{\n";
        file_mesh << blank << blank << blank << "\"Coordinates\" : [\n";
        for (int i = 0; i < nodeCount - 1; ++i)
        {
            file_mesh << blank << blank << blank << blank << "[ ";
            for (int j = 0; j < 2; ++j)
            {
                file_mesh << nodes[i].Coordinates[j] << ", ";
            }
            file_mesh << nodes[i].Coordinates[2] << " ],\n";
        }
        file_mesh << blank << blank << blank << blank << "[ ";
        for (int j = 0; j < 2; ++j)
        {
            file_mesh << nodes[nodeCount - 1].Coordinates[j] << ", ";
        }
        file_mesh << nodes[nodeCount - 1].Coordinates[2] << " ]\n";

        file_mesh << blank << blank << blank << "],\n";
        file_mesh << blank << blank << blank << "\"Indices\" : [";

        for (int i = 0; i < nodeCount - 1; ++i)
        {
            file_mesh << nodes[i].id << ", ";
        }
        file_mesh << nodes[nodeCount-1].id << " ]\n";

        file_mesh << blank << blank << "}\n";
        file_mesh << blank << "],\n";

        file_mesh << blank << "\"DofNodes\" : [\n";
        file_mesh << blank << blank << "{\n";

        file_mesh << blank << blank << "},\n";

        file_mesh << blank << blank << "{\n";
        file_mesh << blank << blank << "}\n";
        file_mesh << blank << "]\n";



        file_mesh << "}";



        file_mesh.close();

    }


    void generateMeshFilesFromGmsh(std::string rFileName, std::string rOutputFile, std::vector<NuTo_DofTypes> rDofTypes, NuTo_InterpolationOrders rInterpolationOrder)
    {
        std::ifstream file_gmsh(rFileName.c_str(), std::ios::in);
        if (not file_gmsh.is_open())
        {
            std::cout << rFileName << std::endl;
            throw std::runtime_error("Error opening input file for read access.");
        }


        GmshHeader header = ReadGmshHeader(file_gmsh);

        std::vector<GmshNode> nodes;
        std::vector<GmshElement> elements;

        if (header.isBinary)
        {
            if (header.double_size != sizeof(double))
                throw std::runtime_error("Invalid size of double.");

            // close rFile and open as binary
            file_gmsh.close();
            file_gmsh.open(rFileName.c_str(), std::ios::in | std::ios::binary);
            if (file_gmsh.is_open() == false)
                throw std::runtime_error("Error opening input rFile for read access.");
            nodes = ReadNodesBinary(file_gmsh);
            elements = ReadElementsBinary(file_gmsh);
        }
        else
        {
            nodes = ReadNodesASCII(file_gmsh);
            elements = ReadElementsASCII(file_gmsh);
        }

        file_gmsh.close();


        const std::string blank = "   ";
        int elementCount = elements.size();
        std::string dofType = "DISPLACEMENTS";
        std::string interpolationOrder = "EQUIDISTANT1";

        int maxDomainID = 1;
        std::vector<std::vector<int>> masterNodeIDs;
        std::vector<std::vector<int>> domainNodeIDs;
//        std::vector<std::vector<std::vector<int>>> interfaceNodeIDs;
        std::vector<std::map<int, std::vector<int>>> interfaceNodeIDs;
        std::vector<int> allNodeIDs;
        std::vector<int> elementCounts(elementCount);
        std::vector<int> nodeCounts(elementCount);
        int subDomainCounter = 0;
        while(maxDomainID > subDomainCounter)
        {
            ++subDomainCounter;
            std::vector<int> newDomainNodeIDs;
            for (GmshElement& element : elements)
            {
                if (subDomainCounter == 1)
                {
                    maxDomainID = std::max(maxDomainID, element.tags.masterDomain);
                }

                if (element.tags.masterDomain == subDomainCounter)
                {
                    elementCounts[element.tags.masterDomain-1] += 1;
                    for (int i = 0; i < element.nodes.size(); ++i)
                    {
                        if (std::find(newDomainNodeIDs.begin(), newDomainNodeIDs.end(), element.nodes[i]) == newDomainNodeIDs.end())
                        {
                            newDomainNodeIDs.push_back(element.nodes[i]);
                            nodeCounts[element.tags.masterDomain-1] += 1;
                        }
                    }
                }
            }
            std::sort(newDomainNodeIDs.begin(), newDomainNodeIDs.end());
            domainNodeIDs.push_back(newDomainNodeIDs);
            std::vector<int> newMasterNodeIDs;
            std::set_difference(newDomainNodeIDs.begin(), newDomainNodeIDs.end(), allNodeIDs.begin(), allNodeIDs.end(), std::back_inserter(newMasterNodeIDs));
            masterNodeIDs.push_back(newMasterNodeIDs);
            allNodeIDs.insert(allNodeIDs.end(), newMasterNodeIDs.begin(), newMasterNodeIDs.end());
            std::sort(allNodeIDs.begin(), allNodeIDs.end());


//            std::vector<std::vector<int>> interfaces;
            std::map<int, std::vector<int>> interfacesMap;
            for (int j = 0; j < subDomainCounter-1; ++j)
            {
                std::vector<int> newInterfaceNodeIDs;
                std::set_intersection(newDomainNodeIDs.begin(), newDomainNodeIDs.end(), masterNodeIDs[j].begin(), masterNodeIDs[j].end(), std::back_inserter(newInterfaceNodeIDs));
                std::sort(newInterfaceNodeIDs.begin(), newInterfaceNodeIDs.end());
                if (newInterfaceNodeIDs.size() > 0)
                {
//                    interfaces.push_back(newInterfaceNodeIDs);
//                    interfacesMap[j+1] = newInterfaceNodeIDs;
                    interfacesMap[j] = newInterfaceNodeIDs;
                }
            }
//            interfaceNodeIDs.push_back(interfaces);
            interfaceNodeIDs.push_back(interfacesMap);

        }

        elementCounts.resize(maxDomainID);
        nodeCounts.resize(maxDomainID);


        for (int subDomain = 1; subDomain < maxDomainID+1; ++subDomain)
        {
            std::string outputFileName =  rOutputFile + ".mff_" + std::to_string(subDomain);

            std::ofstream file_mesh(outputFileName);

            file_mesh << "{\n";
            file_mesh << blank << "\"Elements\" : [\n";
            file_mesh << blank << blank << "{\n";
            file_mesh << blank << blank << blank << "\"Indices\" : [ ";

//            for (int i = 0; i < elementCounts[subDomain-1] - 1; ++i)
            int i = 0;
            for (GmshElement& element : elements)
            {
                if (element.tags.masterDomain == subDomain)
                {
                    if (i < elementCounts[subDomain-1]-1)
                    {
                        file_mesh << element.id << ", ";
                        ++i;
                    }
                    else
                    {
                        file_mesh << element.id << " ],\n";
                        break;
                    }
                }
            }

            file_mesh << blank << blank << blank << "\"NodalConnectivity\" : [\n";

//            for (int i = 0; i < elementCounts[subDomain-1] - 1; ++i)
            i = 0;
            for (GmshElement& element : elements)
            {
                if (element.tags.masterDomain == subDomain)
                {
                    if (i < elementCounts[subDomain-1]-1)
                    {
                        file_mesh << blank << blank << blank << blank << "[ ";
                        for (int j = 0; j < element.nodes.size()-1; ++j)
                        {
                            file_mesh << element.nodes[j] << ", ";
                        }
                        file_mesh << element.nodes[element.nodes.size()-1] << " ],\n";

                        ++i;
                    }
                    else
                    {
                        file_mesh << blank << blank << blank << blank << "[ ";
                        for (int j = 0; j < element.nodes.size()-1; ++j)
                        {
                            file_mesh << element.nodes[j] << ", ";
                        }
                        file_mesh << element.nodes[element.nodes.size()-1] << " ]\n";
                        break;
                    }
                }
            }

            file_mesh << blank << blank << blank << "],\n";

            file_mesh << blank << blank << blank << "\"Type\" : " << elements[0].type << ",\n";     //TODO: different element types possible?
            file_mesh << blank << blank << blank << "\"TypeName\" : " << "\"QUAD2D\"" << "\n";     //TODO: conversion Type -> TypeName

            file_mesh << blank << blank << "}\n";

            file_mesh << blank << "],\n";

            file_mesh << blank << "\"Nodes\" : [\n";
            file_mesh << blank << blank << "{\n";
            file_mesh << blank << blank << blank << "\"Coordinates\" : [\n";
//            for (int i = 0; i < nodeCounts[subDomain-1] - 1; ++i)
            i = 0;
            for (GmshNode& node : nodes)
            {
                if (std::find(domainNodeIDs[subDomain-1].begin(), domainNodeIDs[subDomain-1].end(), node.id) != domainNodeIDs[subDomain-1].end())
                {
                    if (i < nodeCounts[subDomain-1]-1)
                    {
                        file_mesh << blank << blank << blank << blank << "[ ";
                        for (int j = 0; j < 2; ++j)
                        {
                            file_mesh << node.Coordinates[j] << ", ";
                        }
                        file_mesh << node.Coordinates[2] << " ],\n";

                        ++i;
                    }
                    else
                    {
                        file_mesh << blank << blank << blank << blank << "[ ";
                        for (int j = 0; j < 2; ++j)
                        {
                            file_mesh << node.Coordinates[j] << ", ";
                        }
                        file_mesh << node.Coordinates[2] << " ]\n";
                        break;
                    }
                }
            }

            file_mesh << blank << blank << blank << "],\n";
            file_mesh << blank << blank << blank << "\"Indices\" : [ ";

//            for (int i = 0; i < nodeCounts[subDomain-1] - 1; ++i)
            i = 0;
            for (GmshNode& node : nodes)
            {
                if (std::find(domainNodeIDs[subDomain-1].begin(), domainNodeIDs[subDomain-1].end(), node.id) != domainNodeIDs[subDomain-1].end())
                {
                    if (i < nodeCounts[subDomain-1]-1)
                    {
                        file_mesh << node.id << ", ";
                        ++i;
                    }
                    else
                    {
                        file_mesh << node.id << " ]\n";
                        break;
                    }
                }
            }

            file_mesh << blank << blank << "}\n";
            file_mesh << blank << "],\n";

            file_mesh << blank << "\"DofNodes\" : [\n";
            file_mesh << blank << blank << "{\n";

            file_mesh << blank << blank << blank << "\"DofType\" : \"" << dofType << "\",\n";
            file_mesh << blank << blank << blank << "\"InterpolationOrder\": \"" << interpolationOrder << "\",\n";

            file_mesh << blank << blank << blank << "\"MasterNodeIDs\" : [ ";
            for (int j = 0; j < masterNodeIDs[subDomain-1].size()-1; ++j)
            {
                file_mesh << masterNodeIDs[subDomain-1][j] << ", ";
            }
            file_mesh << masterNodeIDs[subDomain-1][masterNodeIDs[subDomain-1].size()-1] << " ],\n";

            file_mesh << blank << blank << blank << "\"Interfaces\" : [ ";
            int interfaceCount = interfaceNodeIDs[subDomain-1].size();
            int interfaceNodeCount = 0;
            i = 0;
//            for (int j = 0; j < interfaceCount-1; ++j)
            for (std::pair<int, std::vector<int>> interfaceMap : interfaceNodeIDs[subDomain-1])
            {
                file_mesh << "\n" << blank << blank << blank << blank << "{\n";
//                file_mesh << blank << blank << blank << blank << blank << "\"Master\" : " << j << ",\n";
                file_mesh << blank << blank << blank << blank << blank << "\"Master\" : " << interfaceMap.first << ",\n";
                file_mesh << blank << blank << blank << blank << blank << "\"NodeIDs\" : [ ";

//                interfaceNodeCount = interfaceNodeIDs[subDomain-1][j].size();
                interfaceNodeCount = interfaceMap.second.size();
                for (int k = 0; k < interfaceNodeCount-1; ++k)
                {
//                    file_mesh << interfaceNodeIDs[subDomain-1][j][k] << ", ";
                    file_mesh << interfaceMap.second[k] << ", ";
                }
//                file_mesh << interfaceNodeIDs[subDomain-1][j][interfaceNodeCount-1] << " ]\n";
                file_mesh << interfaceMap.second[interfaceNodeCount-1] << " ]\n";
                file_mesh << blank << blank << blank << blank << "}";
                ++i;
                if (i < interfaceCount)
                    file_mesh << ",";
                else
                {
                    file_mesh << "\n";
                    file_mesh << blank << blank << blank;
                }
            }
//            if (interfaceCount > 0)
//            {
//                file_mesh << "\n" << blank << blank << blank << blank << "{\n";
//                file_mesh << blank << blank << blank << blank << blank << "\"Master\" : " << interfaceCount-1 << ",\n";
//                file_mesh << blank << blank << blank << blank << blank << "\"NodeIDs\" : [ ";

//                interfaceNodeCount = interfaceNodeIDs[subDomain-1][interfaceCount-1].size();
//                for (int k = 0; k < interfaceNodeCount-1; ++k)
//                {
//                    file_mesh << interfaceNodeIDs[subDomain-1][interfaceCount-1][k] << ", ";
//                }
//                file_mesh << interfaceNodeIDs[subDomain-1][interfaceCount-1][interfaceNodeCount-1] << " ]\n";
//                file_mesh << blank << blank << blank << blank << "}\n";
//                file_mesh << blank << blank << blank;
//            }

            file_mesh << "]\n";

//            file_mesh << blank << blank << "},\n";

//            file_mesh << blank << blank << "{\n";
            file_mesh << blank << blank << "}\n";
            file_mesh << blank << "]\n";



            file_mesh << "}";



            file_mesh.close();

        }
    }



};
