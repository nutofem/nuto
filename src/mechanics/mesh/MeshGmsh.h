#pragma once

#include "MeshFem.h"
#include "base/Group.h"
#include <string>
#include <vector>
#include <map>

class GmshFileContent;

namespace NuTo
{


class MeshGmsh
{
public:
    MeshGmsh() = delete;
    MeshGmsh(const MeshGmsh&) = delete;
    MeshGmsh(MeshGmsh&&) = default;
    MeshGmsh& operator=(const MeshGmsh&) = delete;
    MeshGmsh& operator=(MeshGmsh&&) = delete;
    ~MeshGmsh() = default;

    explicit MeshGmsh(const std::string& fileName);

    const NuTo::Group<ElementCollection>& GetPhysicalGroup(std::string groupName) const;

private:
    void AddElementToPhysicalGroup(ElementCollection& element, std::string physicalName);
    std::vector<NodeSimple*> CreateNodes(const GmshFileContent& fileContent);
    void CreateElements(const GmshFileContent& fileContent, const std::vector<NodeSimple*>& nodePtrVec);
    void CreateMesh(const GmshFileContent& fileContent);
    void ReadGmshFile(const std::string& fileName);

    MeshFem mMesh;
    std::map<std::string, Group<ElementCollection>> mPhysicalGroups;
};
}
