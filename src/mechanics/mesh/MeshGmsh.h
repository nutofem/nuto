#pragma once

#include "MeshFem.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class GmshFileContent;

namespace NuTo
{


class MeshGmsh
{
public:
    explicit MeshGmsh(const std::string& fileName);

private:
    std::vector<NodeSimple*> CreateNodes(const GmshFileContent& fileContent);
    std::vector<ElementCollection*> CreateElements(const GmshFileContent& fileContent,
                                                   const std::vector<NodeSimple*>& nodePtrVec);
    void CreatePhysicalGroups(const std::vector<ElementCollection*>& elementPtrVec);
    void CreateMesh(const GmshFileContent& fileContent);
    void ReadGmshFile(const std::string& fileName);

    MeshFem mMesh;
};
}
