#include "mechanics/mesh/MeshFemDofConvert.h"
#include "math/SpatialContainer.h"

struct DummyNode
{
    Eigen::VectorXd coords;
    NuTo::NodeSimple* existingNodePtr = nullptr;
};

struct DummyNodeCoordinate
{
    Eigen::VectorXd operator()(const DummyNode& node) const
    {
        return node.coords;
    }
};


// create a vector containing all coordinates of the new nodes. Contains a lot of duplicate dummy nodes.
std::vector<DummyNode> GetNewDummyNodes(const NuTo::MeshFem& mesh, const NuTo::InterpolationSimple& interpolation)
{
    std::vector<DummyNode> dummyNodes;
    for (auto& elementCollection : mesh.Elements)
    {
        const auto& coordinateElement = elementCollection.CoordinateElement();
        for (int iNode = 0; iNode < interpolation.GetNumNodes(); ++iNode)
        {
            DummyNode n;
            n.coords = Interpolate(coordinateElement, interpolation.GetLocalCoords(iNode));
            dummyNodes.push_back(n);
        }
    }
    return dummyNodes;
}

void NuTo::AddDofInterpolation(MeshFem* rMesh, DofType dofType, const InterpolationSimple& interpolation)
{
    // We will build a spatial container containing many duplicate DummyNodes at a
    // given coordinate. 
    auto dummyNodes = GetNewDummyNodes(*rMesh, interpolation);
    NuTo::SpatialContainer<DummyNode, DummyNodeCoordinate> dummyNodesSpatial(dummyNodes);

    // Plan:
    // If those dummynodes have _NO_ existingNodePtr, a new node
    // is added to the rMesh->Nodes container. A reference to rMesh newly created node
    // is saved at all of those duplicate dummy nodes.

    // so here we go:
    for (auto& elementCollection : rMesh->Elements)
    {
        std::vector<NodeSimple*> nodesForTheNewlyCreatedElement;

        const auto& coordinateElement = elementCollection.CoordinateElement();
        for (int iNode = 0; iNode < interpolation.GetNumNodes(); ++iNode)
        {
            Eigen::VectorXd coord = Interpolate(coordinateElement, interpolation.GetLocalCoords(iNode));
            auto idsAtCoord = dummyNodesSpatial.FindIDsWithinRadius(coord, 1.e-6);

            // Either _ALL_ dummy nodes at rMesh coordinate have an existingNodePtr or _NONE_. So pick one.
            int randomIdAtCoord = idsAtCoord.front();
            NodeSimple* existingNode = dummyNodes[randomIdAtCoord].existingNodePtr;

            if (existingNode == nullptr)
            {
                // create new node and ...
                existingNode = &rMesh->Nodes.Add(Eigen::VectorXd::Zero(dofType.GetNum()));

                // ... put it at _ALL_ dummy node at that coordinate
                for (int id : idsAtCoord)
                    dummyNodes[id].existingNodePtr = existingNode;
            }
            // existingNode now contains a pointer to either
            //     a newly created node or
            //     an existing node pointer at coord, extracted from the spatial container
            //
            //  Boom, that's it. 

            nodesForTheNewlyCreatedElement.push_back(existingNode);
        }
        elementCollection.AddDofElement(dofType, ElementFem(nodesForTheNewlyCreatedElement, interpolation));
    }

    // This main reason why this is still ulgy: The SpatialContainer (because of ANN) does not allow adding
    // entities once it is created. If this restriction is gone, it'll look much cleaner and will only require
    // one loop over all elements and their nodes. If the syntax is similar to STL container, the inner part 
    // could look like:
    //
    //  SomeIterator it = superSpatialContainer.at(coord)
    //  if (it != superSpatialContainer.end())
    //      nodesForTheNewlyCreatedElement[iNode] = *it;
    //  else
    //      nodesForTheNewlyCreatedElement[iNode] = superSpatialContaine.insert(rMesh->Nodes.Add(...));
    //
    // Not much shorter, but this does not require load of comments to get.

}
