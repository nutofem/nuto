#include "mechanics/mesh/MeshFem.h"

#include <sstream>
#include "math/SpatialContainer.h"
#include "base/Exception.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"

using namespace NuTo;

InterpolationSimple& MeshFem::CreateInterpolation(const InterpolationSimple& interpolation)
{
    mInterpolations.push_back(interpolation.Clone().release());
    return *mInterpolations.rbegin();
}

NodeSimple& MeshFem::NodeAtCoordinate(Eigen::VectorXd coords, DofType dofType, double tol /* = 1.e-10 */)
{
    for (auto& element : this->Elements)
    {
        auto& dofElement = element.DofElement(dofType);
        const auto& dofInterpolation = dofElement.Interpolation();
        for (int iNode = 0; iNode < dofInterpolation.GetNumNodes(); ++iNode)
        {
            NaturalCoords dofNodeCoords = dofInterpolation.GetLocalCoords(iNode);
            Eigen::VectorXd globalNodeCoords = Interpolate(element.CoordinateElement(), dofNodeCoords);
            if ((globalNodeCoords - coords).isMuchSmallerThan(tol))
                return dofElement.GetNode(iNode);
        }
    }
    std::stringstream coordsString;
    coordsString << coords.transpose();
    throw NuTo::Exception(__PRETTY_FUNCTION__,
                          "There is no node for dof type " + dofType.GetName() + " at " + coordsString.str());
}

NodeSimple& MeshFem::NodeAtCoordinate(Eigen::VectorXd coords, double tol /* = 1.e-10 */)
{
    for (auto& element : this->Elements)
    {
        auto& coordinateElement = element.CoordinateElement();
        for (int iNode = 0; iNode < coordinateElement.Interpolation().GetNumNodes(); ++iNode)
        {
            Eigen::VectorXd globalNodeCoords = coordinateElement.GetNode(iNode).GetValues();
            if ((globalNodeCoords - coords).isMuchSmallerThan(tol))
                return coordinateElement.GetNode(iNode);
        }
    }
    std::stringstream coordsString;
    coordsString << coords.transpose();
    throw NuTo::Exception(__PRETTY_FUNCTION__, "There is no coordinate node at " + coordsString.str());
}

Groups::Group<NodeSimple> MeshFem::NodesAtAxis(eDirection direction, DofType dofType, double axisOffset /* = 0.*/,
                                               double tol /* = 1.e-10 */)
{
    Groups::Group<NodeSimple> group;
    const int directionComponent = ToComponentIndex(direction);
    for (auto& element : this->Elements)
    {
        auto& dofElement = element.DofElement(dofType);
        const auto& dofInterpolation = dofElement.Interpolation();
        for (int iNode = 0; iNode < dofInterpolation.GetNumNodes(); ++iNode)
        {
            NaturalCoords dofNodeCoords = dofInterpolation.GetLocalCoords(iNode);
            Eigen::VectorXd globalNodeCoords = Interpolate(element.CoordinateElement(), dofNodeCoords);

            if (std::abs(globalNodeCoords[directionComponent] - axisOffset) < tol)
                group.Add(dofElement.GetNode(iNode));
        }
    }
    return group;
}

Groups::Group<NodeSimple> MeshFem::NodesTotal()
{
    Groups::Group<NodeSimple> group;
    for (auto& element : this->Elements)
        for (int iNode = 0; iNode < element.CoordinateElement().Interpolation().GetNumNodes(); ++iNode)
            group.Add(element.CoordinateElement().GetNode(iNode));
    return group;
}

Groups::Group<NodeSimple> MeshFem::NodesTotal(DofType d)
{
    Groups::Group<NodeSimple> group;
    for (auto& element : this->Elements)
        for (int iNode = 0; iNode < element.DofElement(d).Interpolation().GetNumNodes(); ++iNode)
            group.Add(element.DofElement(d).GetNode(iNode));
    return group;
}

void MeshFem::AddDofInterpolation(DofType dofType, const InterpolationSimple& interpolation)
{
    // We will build a spatial container containing many duplicate DummyNodes at a 
    // given coordinate. If those dummynodes have _NO_ existingNodePtr, a new node 
    // is added to the this->Nodes container. A reference to this newly created node
    // is saved at all of those duplicate dummy nodes. 
    struct DummyNode
    {
        Eigen::VectorXd coords;
        NodeSimple* existingNodePtr = nullptr;
    };

    struct DummyNodeCoordinate
    {
        Eigen::VectorXd operator()(const DummyNode& node) const
        {
            return node.coords;
        }
    };

    // create a vector containing all coordinates of the new nodes. This contains a lot of duplicate dummy nodes.
    std::vector<DummyNode> dummyNodes;
    for (auto& elementCollection : this->Elements)
    {
        const auto& coordinateElement = elementCollection.CoordinateElement();
        for (int iNode = 0; iNode < interpolation.GetNumNodes(); ++iNode)
        {
            DummyNode n;
            n.coords = Interpolate(coordinateElement, interpolation.GetLocalCoords(iNode));
            dummyNodes.push_back(n);
        }
    }

    NuTo::SpatialContainer<DummyNode, DummyNodeCoordinate> dummyNodesSpatial(dummyNodes);
    for (auto& elementCollection : this->Elements)
    {
        const auto& coordinateElement = elementCollection.CoordinateElement();
        std::vector<NodeSimple*> dofNodes;
        for (int iNode = 0; iNode < interpolation.GetNumNodes(); ++iNode)
        {
            Eigen::VectorXd coord = Interpolate(coordinateElement, interpolation.GetLocalCoords(iNode));
            auto idsAtCoord = dummyNodesSpatial.FindIDsWithinRadius(coord, 1.e-6);

            // Either _ALL_ dummy nodes at this coordinate have an existingNodePtr or _NONE_. So pick one.
            NodeSimple* existingNode = dummyNodes[idsAtCoord.front()].existingNodePtr;
            if (existingNode == nullptr)
            {
                // create new node and ...
                existingNode = &this->Nodes.Add(Eigen::VectorXd::Zero(dofType.GetNum()));

                // ... put it at _ALL_ dummy node at that coordinate
                for (int id : idsAtCoord)
                    dummyNodes[id].existingNodePtr = existingNode;
            }
            dofNodes.push_back(existingNode);
        }
        elementCollection.AddDofElement(dofType, ElementFem(dofNodes, interpolation));
    }
}


// anonymous helper functions for mesh creation and transformation

MeshFem CreateNodes2D(int numX, int numY)
{
    MeshFem mesh;
    for (int iY = 0; iY < numY + 1; ++iY)
        for (int iX = 0; iX < numX + 1; ++iX)
        {
            const double x = static_cast<double>(iX) / numX;
            const double y = static_cast<double>(iY) / numY;
            mesh.Nodes.Add(Eigen::Vector2d(x, y));
        }
    return mesh;
}

void TransformNodes(Groups::Group<NodeSimple> nodes, std::function<Eigen::VectorXd(Eigen::VectorXd)> f)
{
    for (auto& node : nodes)
        node.SetValues(f(node.GetValues()));
}

// UnitMeshFem functions

MeshFem UnitMeshFem::CreateTriangles(int numX, int numY)
{
    MeshFem mesh = CreateNodes2D(numX, numY);
    const auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear(2));
    for (int iY = 0; iY < numY; ++iY)
        for (int iX = 0; iX < numX; ++iX)
        {
            auto& node0 = mesh.Nodes[iX + iY * (numX + 1)];
            auto& node1 = mesh.Nodes[iX + 1 + iY * (numX + 1)];
            auto& node2 = mesh.Nodes[iX + 1 + (iY + 1) * (numX + 1)];
            auto& node3 = mesh.Nodes[iX + (iY + 1) * (numX + 1)];
            mesh.Elements.Add({{{node0, node1, node2}, interpolation}});
            mesh.Elements.Add({{{node0, node2, node3}, interpolation}});
        }
    return mesh;
}

MeshFem UnitMeshFem::CreateQuads(int numX, int numY)
{
    MeshFem mesh = CreateNodes2D(numX, numY);
    const auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationQuadLinear(2));
    for (int iY = 0; iY < numY; ++iY)
        for (int iX = 0; iX < numX; ++iX)
        {
            auto& node0 = mesh.Nodes[iX + iY * (numX + 1)];
            auto& node1 = mesh.Nodes[iX + 1 + iY * (numX + 1)];
            auto& node2 = mesh.Nodes[iX + 1 + (iY + 1) * (numX + 1)];
            auto& node3 = mesh.Nodes[iX + (iY + 1) * (numX + 1)];
            mesh.Elements.Add({{{node0, node1, node2, node3}, interpolation}});
        }
    return mesh;
}

void UnitMeshFem::Transform(MeshFem* rMesh, std::function<Eigen::VectorXd(Eigen::VectorXd)> f)
{
    auto nodes = rMesh->NodesTotal();
    TransformNodes(nodes, f);
}
