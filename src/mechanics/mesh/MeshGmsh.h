#pragma once

#include <fstream>
#include <string>
#include <vector>

struct GmshNode;
struct GmshElement;

namespace NuTo
{


class MeshGmsh
{
public:
    explicit MeshGmsh(const std::string& fileName);


    unsigned int GetNumNodes() const
    {
        return mNumNodes;
    }

    unsigned int GetNumElements() const
    {
        return mNumElements;
    }

private:
    void ReadGmshFile(const std::string& fileName);

    std::vector<GmshNode> ReadNodesASCII(std::ifstream& file);
    std::vector<GmshElement> ReadElementsASCII(std::ifstream& file);

    unsigned int mNumNodes = 0;
    unsigned int mNumElements = 0;
};
}
