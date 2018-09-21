#pragma once

#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/base/Group.h"
#include <map>

namespace NuTo
{

class MeshTopology
{
public:
    MeshTopology(Group<ElementCollectionFem> elements);

    int GetNumEdges()
    {
        return mEdgesToVertices.size();
    }

    int GetNumFaces()
    {
        return mFacesToEdges.size();
    }

    int GetNumVolumes()
    {
        return mVolumesToFaces.size();
    }

    //! @brief Returns edge elements that are shared with given element
    //! @param elm an element
    //! @return Group of edge elements (elements of local dimension 1)
    //! @remark If elm itself is an edge it is included in the result
    Group<ElementCollectionFem> GetAdjacentEdges(ElementCollectionFem& elm);

    //! @brief Returns face elements that are shared with given element
    //! @param elm an element
    //! @return Group of face elements (elements of local dimension 2)
    //! @remark If elm itself is a face it is included in the result
    Group<ElementCollectionFem> GetAdjacentFaces(ElementCollectionFem& elm);

    //! @brief Returns volume elements that are shared with given element
    //! @param elm an element
    //! @return Group of volume elements (elements of local dimension 3)
    //! @remark If elm itself is a volume it is included in the result
    Group<ElementCollectionFem> GetAdjacentVolumes(ElementCollectionFem& elm);

private:
    std::map<NodeSimple*, std::vector<ElementCollectionFem*>> mVerticesToEdges;
    std::map<ElementCollectionFem*, std::vector<NodeSimple*>> mEdgesToVertices;

    std::map<ElementCollectionFem*, std::vector<ElementCollectionFem*>> mFacesToEdges;
    std::map<ElementCollectionFem*, std::vector<ElementCollectionFem*>> mEdgesToFaces;

    std::map<ElementCollectionFem*, std::vector<ElementCollectionFem*>> mVolumesToFaces;
    std::map<ElementCollectionFem*, std::vector<ElementCollectionFem*>> mFacesToVolumes;
};
} /* NuTo */
