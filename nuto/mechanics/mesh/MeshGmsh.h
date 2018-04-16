#pragma once

#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/base/Group.h"
#include <string>
#include <map>
#include <unordered_map>

struct GmshFileContent;

namespace NuTo
{

//! @brief Reads a gmsh msh-file and converts it to a MeshFEM. Additionally physical groups are created and can be
//! accessed via ID or name
class MeshGmsh
{
public:
    //! @brief ctor
    //! @remark deleted: The class is always constructed with its associated files name
    MeshGmsh() = delete;

    //! @brief copy ctor
    //! @remark deleted: Member MeshFEM is not copyable
    MeshGmsh(const MeshGmsh&) = delete;

    //! @brief move ctor
    MeshGmsh(MeshGmsh&&) = default;

    //! @brief copy assignment operator
    //! @remark deleted: Member MeshFEM is not copyable
    MeshGmsh& operator=(const MeshGmsh&) = delete;

    //! @brief move assignment operator
    MeshGmsh& operator=(MeshGmsh&&) = default;

    //! @brief dtor
    ~MeshGmsh() = default;

    //! @brief ctor
    //! @param fileName Name of the file (including path)
    explicit MeshGmsh(const std::string& fileName);

    //! @brief Gets an element collection group associated to the physical name defined in gmsh
    //! @param physicalName Group name (defined in gmsh)
    //! @return Desired element collection group
    const NuTo::Group<ElementCollectionFem>& GetPhysicalGroup(std::string physicalName) const;


    //! @brief Gets an element collection group associated to the physical ID defined in gmsh
    //! @param physicalGroupId gmsh ID
    //! @return Desired element collection group
    const NuTo::Group<ElementCollectionFem>& GetPhysicalGroup(int physicalGroupId) const;

    //! @brief Gets the MeshFem
    MeshFem& GetMeshFEM()
    {
        return mMesh;
    }

private:
    //! @brief Adds an element to its corresponding physical group. If the element is part of named group and the group
    //! does not exist, then a named physical group is also created.
    //! @param fileContent Special structure (see MeshGmsh.cpp) that holds the read file content
    //! @param rElement Current element
    //! @param physicalGroupId Id of the elements physical group
    void AddElementToPhysicalGroup(const GmshFileContent& fileContent, ElementCollectionFem& rElement,
                                   int physicalGroupId);

    //! @brief Creates all nodes defined in the gmsh file
    //! @param fileContent Special structure (see MeshGmsh.cpp) that holds the read file content
    //! @return Vector containing pointers to all nodes
    //! @remark The node dimension is determined automatically. If at least one node has a y-coordinate the dimension is
    //! set to 2d (or 3d). If at least one node has a z-coordinate the dimension is set to 3d.
    //! @remark The returned collection is used in CreateElements to easily find the nodes connected to a specific gmsh
    //! ID.
    std::unordered_map<int, CoordinateNode*> CreateNodes(const GmshFileContent& fileContent);

    //! @brief Creates all elements defined in the gmsh file
    //! @param fileContent Special structure (see MeshGmsh.cpp) that holds the read file content
    //! @param nodePtrs Collection that holds pointers to all created elements. (See CreateNodes)
    void CreateElements(const GmshFileContent& fileContent, const std::unordered_map<int, CoordinateNode*>& nodePtrs);

    //! @brief Fills the MeshFem member from gmsh file content
    //! @param fileContent Special structure (see MeshGmsh.cpp) that holds the read file content
    void CreateMesh(const GmshFileContent& fileContent);

    //! @brief Reads the file content from a gmsh file and calls CreateMesh
    //! @param fileName Name of the file (including path)
    void ReadGmshFile(const std::string& fileName);

    //! @brief Internal mesh
    MeshFem mMesh;

    //! @brief Physical groups of the mesh
    std::map<int, Group<ElementCollectionFem>> mPhysicalGroups;

    //! @brief Named physical groups of the mesh
    std::map<std::string, Group<ElementCollectionFem>*> mNamedPhysicalGroups;
};
}
