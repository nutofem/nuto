#include "nuto/mechanics/mesh/MeshFem.h"

#include "nuto/mechanics/mesh/GeometryMeshFem.h"

#include <sstream>
#include "nuto/base/Exception.h"


using namespace NuTo;

MeshFem::MeshFem(GeometryMeshFem& geometryMesh)
    : mGeometryMesh(geometryMesh)
{
    for (size_t i = 0; i < geometryMesh.NumElements(); i++)
    {
        auto& cElm = geometryMesh.GetElement(i);
        Elements.Add(ElementCollectionFem(cElm));
    }
}

InterpolationSimple& MeshFem::CreateInterpolation(const InterpolationSimple& interpolation)
{
    mInterpolations.push_back(interpolation.Clone());
    return *mInterpolations.rbegin()->get();
}

DofNode& MeshFem::NodeAtCoordinate(Eigen::VectorXd coords, DofType dofType, double tol /* = 1.e-10 */)
{
    for (auto& element : this->Elements)
    {
        if (!element.Has(dofType))
            continue;

        auto& dofElement = element.DofElement(dofType);
        const auto& dofInterpolation = dofElement.Interpolation();
        for (int iNode = 0; iNode < dofInterpolation.GetNumNodes(); ++iNode)
        {
            NaturalCoords dofNodeCoords = dofInterpolation.GetLocalCoords(iNode);
            Eigen::VectorXd globalNodeCoords = Interpolate(element.CoordinateElement(), dofNodeCoords);
            if ((globalNodeCoords - coords).isMuchSmallerThan(tol, 1))
                return dofElement.GetNode(iNode);
        }
    }
    std::stringstream coordsString;
    coordsString << coords.transpose();
    throw NuTo::Exception(__PRETTY_FUNCTION__,
                          "There is no node for dof type " + dofType.GetName() + " at " + coordsString.str());
}

Group<DofNode> MeshFem::NodesAtAxis(eDirection direction, DofType dofType, double axisOffset /* = 0.*/,
                                    double tol /* = 1.e-10 */)
{
    Group<DofNode> group;
    const int directionComponent = ToComponentIndex(direction);
    for (auto& element : this->Elements)
    {
        if (!element.Has(dofType))
            continue;

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

Group<DofNode> MeshFem::NodesTotal(DofType d)
{
    Group<DofNode> group;
    for (auto& element : this->Elements)
    {
        if (!element.Has(d))
            continue;
        for (int iNode = 0; iNode < element.DofElement(d).Interpolation().GetNumNodes(); ++iNode)
            group.Add(element.DofElement(d).GetNode(iNode));
    }
    return group;
}

Group<ElementCollectionFem> MeshFem::ElementsTotal()
{
    Group<ElementCollectionFem> elements;
    for (auto& element : this->Elements)
        elements.Add(element);
    return elements;
}

void MeshFem::AllocateDofInstances(DofType dofType, int numInstances)
{
    for (auto& node : NodesTotal(dofType))
        node.AllocateInstances(numInstances);
}
