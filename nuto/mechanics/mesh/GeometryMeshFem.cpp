#include "nuto/mechanics/mesh/GeometryMeshFem.h"

#include <sstream>
#include "nuto/base/Exception.h"

using namespace NuTo;

InterpolationSimple& GeometryMeshFem::CreateInterpolation(const InterpolationSimple& interpolation)
{
    mInterpolations.push_back(interpolation.Clone());
    return *mInterpolations.rbegin()->get();
}

CoordinateNode& GeometryMeshFem::NodeAtCoordinate(Eigen::VectorXd coords, double tol /* = 1.e-10 */)
{
    for (auto& element : this->Elements)
    {
        auto& coordinateElement = element.CoordinateElement();
        for (int iNode = 0; iNode < coordinateElement.Interpolation().GetNumNodes(); ++iNode)
        {
            Eigen::VectorXd globalNodeCoords = coordinateElement.GetNode(iNode).GetCoordinates();
            if ((globalNodeCoords - coords).isMuchSmallerThan(tol, 1))
                return coordinateElement.GetNode(iNode);
        }
    }
    std::stringstream coordsString;
    coordsString << coords.transpose();
    throw NuTo::Exception(__PRETTY_FUNCTION__, "There is no coordinate node at " + coordsString.str());
}

Group<CoordinateNode> GeometryMeshFem::NodesAtAxis(eDirection direction, double axisOffset /* = 0.*/,
                                                   double tol /* = 1.e-10 */)
{
    Group<CoordinateNode> group;
    const int directionComponent = ToComponentIndex(direction);
    for (auto& element : this->Elements)
    {
        auto& coordinateElement = element.CoordinateElement();
        for (int iNode = 0; iNode < coordinateElement.GetNumNodes(); ++iNode)
        {
            Eigen::VectorXd globalNodeCoords = coordinateElement.GetNode(iNode).GetCoordinates();
            if (std::abs(globalNodeCoords[directionComponent] - axisOffset) < tol)
                group.Add(coordinateElement.GetNode(iNode));
        }
    }
    return group;
}

Group<CoordinateNode> GeometryMeshFem::NodesTotal()
{
    Group<CoordinateNode> group;
    for (auto& element : this->Elements)
        for (int iNode = 0; iNode < element.CoordinateElement().Interpolation().GetNumNodes(); ++iNode)
            group.Add(element.CoordinateElement().GetNode(iNode));
    return group;
}

Group<ElementCollectionFem> GeometryMeshFem::ElementsTotal()
{
    Group<ElementCollectionFem> elements;
    for (auto& element : this->Elements)
        elements.Add(element);
    return elements;
}
