#pragma once

#include "nuto/mechanics/mesh/MeshFem.h"
#include <map>
#include <set>


namespace NuTo
{

class MeshTopology
{
public:
    MeshTopology(MeshFem* rMesh)
    {
        Group<ElementCollectionFem> edgeElements;
        Group<ElementCollectionFem> faceElements;
        Group<ElementCollectionFem> volumeElements;

        for (ElementCollectionFem& elmColl : rMesh->Elements)
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
                    mVerticesToFaces[&elmColl.CoordinateElement().GetNode(i)].push_back(&elmColl);
                }
                break;
            }
            case 3:
            {
                volumeElements.Add(elmColl);
                for (int i = 0; i < elmColl.CoordinateElement().GetNumNodes(); i++)
                {
                    mVerticesToVolumes[&elmColl.CoordinateElement().GetNode(i)].push_back(&elmColl);
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
                std::vector<ElementCollectionFem*> connectedEdges =
                        mVerticesToEdges[&face.CoordinateElement().GetNode(i)];
                for (ElementCollectionFem* edge : connectedEdges)
                {
                    bool edgeIsContained = true;
                    for (int i = 0; i < edge->CoordinateElement().GetNumNodes(); i++)
                    {
                        auto& otherFaces = mVerticesToFaces[&(edge->CoordinateElement().GetNode(i))];
                        if (std::find(otherFaces.begin(), otherFaces.end(), &face) == otherFaces.end())
                        {
                            edgeIsContained = false;
                            break;
                        }
                    }
                    if (edgeIsContained)
                        mFacesToEdges[&face].push_back(edge);
                }
            }
        }

        for (ElementCollectionFem& volume : volumeElements)
        {
            for (int i = 0; i < volume.CoordinateElement().GetNumNodes(); i++)
            {
                std::vector<ElementCollectionFem*> connectedFaces =
                        mVerticesToFaces[&volume.CoordinateElement().GetNode(i)];
                for (ElementCollectionFem* face : connectedFaces)
                {
                    bool faceIsContained = true;
                    for (int i = 0; i < face->CoordinateElement().GetNumNodes(); i++)
                    {
                        auto& otherVolumes = mVerticesToVolumes[&(face->CoordinateElement().GetNode(i))];
                        if (std::find(otherVolumes.begin(), otherVolumes.end(), &volume) == otherVolumes.end())
                        {
                            faceIsContained = false;
                            break;
                        }
                    }
                    if (faceIsContained)
                        mVolumesToFaces[&volume].push_back(face);
                }
            }
        }
    }

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

private:
    std::map<NodeSimple*, std::vector<ElementCollectionFem*>> mVerticesToEdges;
    std::map<NodeSimple*, std::vector<ElementCollectionFem*>> mVerticesToFaces;
    std::map<NodeSimple*, std::vector<ElementCollectionFem*>> mVerticesToVolumes;

    std::map<ElementCollectionFem*, std::vector<NodeSimple*>> mEdgesToVertices;
    std::map<ElementCollectionFem*, std::vector<ElementCollectionFem*>> mFacesToEdges;
    std::map<ElementCollectionFem*, std::vector<ElementCollectionFem*>> mVolumesToFaces;
};
} /* NuTo */
