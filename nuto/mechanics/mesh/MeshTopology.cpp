#include "nuto/mechanics/mesh/MeshTopology.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/base/Group.h"
#include <map>
#include <set>


namespace NuTo
{

MeshTopology::MeshTopology(Group<ElementCollectionFem> elements)
{
    Group<ElementCollectionFem> edgeElements;
    Group<ElementCollectionFem> faceElements;
    Group<ElementCollectionFem> volumeElements;

    std::map<NodeSimple*, std::vector<ElementCollectionFem*>> verticesToFaces;
    std::map<NodeSimple*, std::vector<ElementCollectionFem*>> verticesToVolumes;

    for (ElementCollectionFem& elmColl : elements)
    {
        int dim = elmColl.CoordinateElement().Interpolation().GetLocalCoords(0).size();
        switch (dim)
        {
        case 1:
        {
            edgeElements.Add(elmColl);
            for (int i = 0; i < elmColl.CoordinateElement().GetNumNodes(); i++)
            {
                mVerticesToEdges[&elmColl.CoordinateElement().GetNode(i)].push_back(&elmColl);
                mEdgesToVertices[&elmColl].push_back(&elmColl.CoordinateElement().GetNode(i));
            }
            break;
        }
        case 2:
        {
            faceElements.Add(elmColl);
            for (int i = 0; i < elmColl.CoordinateElement().GetNumNodes(); i++)
            {
                verticesToFaces[&elmColl.CoordinateElement().GetNode(i)].push_back(&elmColl);
            }
            break;
        }
        case 3:
        {
            volumeElements.Add(elmColl);
            for (int i = 0; i < elmColl.CoordinateElement().GetNumNodes(); i++)
            {
                verticesToVolumes[&elmColl.CoordinateElement().GetNode(i)].push_back(&elmColl);
            }
            break;
        }
        default:
            throw Exception(__PRETTY_FUNCTION__, "Dimension out of range (1-3).");
        }
    }

    for (ElementCollectionFem& face : faceElements)
    {
        for (int i = 0; i < face.CoordinateElement().GetNumNodes(); i++)
        {
            std::vector<ElementCollectionFem*> connectedEdges = mVerticesToEdges[&face.CoordinateElement().GetNode(i)];
            for (ElementCollectionFem* edge : connectedEdges)
            {
                bool edgeIsContained = true;
                for (int i = 0; i < edge->CoordinateElement().GetNumNodes(); i++)
                {
                    auto& otherFaces = verticesToFaces[&(edge->CoordinateElement().GetNode(i))];
                    if (std::find(otherFaces.begin(), otherFaces.end(), &face) == otherFaces.end())
                    {
                        edgeIsContained = false;
                        break;
                    }
                }
                if (edgeIsContained)
                {
                    mFacesToEdges[&face].push_back(edge);
                    mEdgesToFaces[edge].push_back(&face);
                }
            }
        }
    }

    for (ElementCollectionFem& volume : volumeElements)
    {
        for (int i = 0; i < volume.CoordinateElement().GetNumNodes(); i++)
        {
            std::vector<ElementCollectionFem*> connectedFaces = verticesToFaces[&volume.CoordinateElement().GetNode(i)];
            for (ElementCollectionFem* face : connectedFaces)
            {
                bool faceIsContained = true;
                for (int i = 0; i < face->CoordinateElement().GetNumNodes(); i++)
                {
                    auto& otherVolumes = verticesToVolumes[&(face->CoordinateElement().GetNode(i))];
                    if (std::find(otherVolumes.begin(), otherVolumes.end(), &volume) == otherVolumes.end())
                    {
                        faceIsContained = false;
                        break;
                    }
                }
                if (faceIsContained)
                {
                    mVolumesToFaces[&volume].push_back(face);
                    mFacesToVolumes[face].push_back(&volume);
                }
            }
        }
    }
}

Group<ElementCollectionFem> MeshTopology::GetAdjacentEdges(ElementCollectionFem& elm)
{
    Group<ElementCollectionFem> result;
    int dim = elm.CoordinateElement().Interpolation().GetLocalCoords(0).size();
    switch (dim)
    {
    case 1:
    {
        for (NodeSimple* node : mEdgesToVertices[&elm])
            for (ElementCollectionFem* edge : mVerticesToEdges[node])
                result.Add(*edge);
        return result;
    }
    case 2:
    {
        for (ElementCollectionFem* edge : mFacesToEdges[&elm])
            result.Add(*edge);
        return result;
    }
    case 3:
    {
        for (ElementCollectionFem* face : mVolumesToFaces[&elm])
            for (ElementCollectionFem* edge : mFacesToEdges[face])
                result.Add(*edge);
        return result;
    }
    default:
        throw Exception(__PRETTY_FUNCTION__, "Dimension out of range (1-3).");
    }
}

Group<ElementCollectionFem> MeshTopology::GetAdjacentFaces(ElementCollectionFem& elm)
{
    Group<ElementCollectionFem> result;
    int dim = elm.CoordinateElement().Interpolation().GetLocalCoords(0).size();
    switch (dim)
    {
    case 1:
    {
        for (ElementCollectionFem* face : mEdgesToFaces[&elm])
            result.Add(*face);
        return result;
    }
    case 2:
    {
        for (ElementCollectionFem* edge : mFacesToEdges[&elm])
            for (ElementCollectionFem* face : mEdgesToFaces[edge])
                result.Add(*face);
        return result;
    }
    case 3:
    {
        for (ElementCollectionFem* face : mVolumesToFaces[&elm])
            result.Add(*face);
        return result;
    }
    default:
        throw Exception(__PRETTY_FUNCTION__, "Dimension out of range (1-3).");
    }
}

Group<ElementCollectionFem> MeshTopology::GetAdjacentVolumes(ElementCollectionFem& elm)
{
    Group<ElementCollectionFem> result;
    int dim = elm.CoordinateElement().Interpolation().GetLocalCoords(0).size();
    switch (dim)
    {
    case 1:
    {
        for (ElementCollectionFem* face : mEdgesToFaces[&elm])
            for (ElementCollectionFem* volume : mFacesToVolumes[face])
                result.Add(*volume);
        return result;
    }
    case 2:
    {
        for (ElementCollectionFem* volume : mFacesToVolumes[&elm])
            result.Add(*volume);
        return result;
    }
    case 3:
    {
        for (ElementCollectionFem* face : mVolumesToFaces[&elm])
            for (ElementCollectionFem* volume : mFacesToVolumes[face])
                result.Add(*volume);
        return result;
    }
    default:
        throw Exception(__PRETTY_FUNCTION__, "Dimension out of range (1-3).");
    }
}
} /* NuTo */
