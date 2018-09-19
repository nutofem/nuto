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
    for (auto& node : this->CoordinateNodes)
    {
        Eigen::VectorXd globalNodeCoords = node.GetCoordinates();
        if ((globalNodeCoords - coords).isMuchSmallerThan(tol, 1))
            return node;
    }
    std::stringstream message;
    message << "There is no coordinate node at " << coords.transpose();
    throw NuTo::Exception(__PRETTY_FUNCTION__, message.str());
}

Group<CoordinateNode> GeometryMeshFem::NodesAtAxis(eDirection direction, double axisOffset /* = 0.*/,
                                                   double tol /* = 1.e-10 */)
{
    Group<CoordinateNode> group;
    const int directionComponent = ToComponentIndex(direction);
    for (auto& node : this->CoordinateNodes)
    {
        Eigen::VectorXd globalNodeCoords = node.GetCoordinates();
        if (std::abs(globalNodeCoords[directionComponent] - axisOffset) < tol)
            group.Add(node);
    }
    return group;
}

Group<CoordinateNode> GeometryMeshFem::NodesTotal()
{
    Group<CoordinateNode> group;
    for (auto& node : this->CoordinateNodes)
        group.Add(node);
    return group;
}

Group<CoordinateElementFem> GeometryMeshFem::ElementsTotal()
{
    Group<CoordinateElementFem> elements;
    for (auto& element : this->Elements)
        elements.Add(element);
    return elements;
}
