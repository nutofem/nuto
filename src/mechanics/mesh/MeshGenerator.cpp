#include "mechanics/mesh/MeshGenerator.h"
#include <cassert>

#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/structures/unstructured/Structure.h"

using namespace NuTo::Interpolation;

std::pair<int, int> CreateInterpolationTypeAndGroup(NuTo::Structure& s, eShapeType elementShape)
{
    int interpolationType = s.InterpolationTypeCreate(elementShape);
    s.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    int elementGroup = s.GroupCreate(NuTo::eGroupId::Elements);
    return std::make_pair(elementGroup, interpolationType);
}

std::vector<std::vector<int>> GetElementNodeIds3D(const std::vector<int>& cornerNodes, eShapeType elementShape)
{
    switch (elementShape)
    {
    case eShapeType::BRICK3D:
    {
        return {cornerNodes};
    }
    case eShapeType::TETRAHEDRON3D:
    {
        std::vector<int> nodesTet0({cornerNodes[0], cornerNodes[1], cornerNodes[3], cornerNodes[7]});
        std::vector<int> nodesTet1({cornerNodes[0], cornerNodes[1], cornerNodes[7], cornerNodes[4]});
        std::vector<int> nodesTet2({cornerNodes[5], cornerNodes[4], cornerNodes[7], cornerNodes[1]});
        std::vector<int> nodesTet3({cornerNodes[6], cornerNodes[5], cornerNodes[7], cornerNodes[1]});
        std::vector<int> nodesTet4({cornerNodes[2], cornerNodes[7], cornerNodes[1], cornerNodes[6]});
        std::vector<int> nodesTet5({cornerNodes[2], cornerNodes[3], cornerNodes[1], cornerNodes[7]});
        return {nodesTet0, nodesTet1, nodesTet2, nodesTet3, nodesTet4, nodesTet5};
    }
    case eShapeType::PRISM3D:
    {
        std::vector<int> nodes0(
                {cornerNodes[0], cornerNodes[1], cornerNodes[2], cornerNodes[4], cornerNodes[5], cornerNodes[6]});
        std::vector<int> nodes1(
                {cornerNodes[0], cornerNodes[2], cornerNodes[3], cornerNodes[4], cornerNodes[6], cornerNodes[7]});
        return {nodes0, nodes1};
    }
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, ShapeTypeToString(elementShape) + " not supported as 3D element");
    }
}

std::vector<std::vector<int>> GetElementNodeIds2D(const std::vector<int>& cornerNodes, eShapeType elementShape)
{
    switch (elementShape)
    {
    case eShapeType::QUAD2D:
    {
        return {cornerNodes};
    }
    case eShapeType::TRIANGLE2D:
    {
        std::vector<int> e1{cornerNodes[0], cornerNodes[1], cornerNodes[2]};
        std::vector<int> e2{cornerNodes[0], cornerNodes[2], cornerNodes[3]};
        return {e1, e2};
    }
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, ShapeTypeToString(elementShape) + " not supported as 3D element");
    }
}

std::vector<int> CreateNodes(NuTo::Structure& s, std::vector<int> numNodes, std::vector<double> delta,
                             std::vector<double> start)
{
    std::vector<int> nodeIds;
    switch (s.GetDimension())
    {
    case 1:
    {
        nodeIds.reserve(numNodes[0]);
        for (int iX = 0; iX < numNodes[0]; ++iX)
        {
            const double x = iX * delta[0] + start[0];
            nodeIds.push_back(s.NodeCreate(Eigen::Matrix<double, 1, 1>::Constant(x)));
        }
        break;
    }
    case 2:
    {
        nodeIds.reserve(numNodes[0] * numNodes[1]);
        for (int iY = 0; iY < numNodes[1]; ++iY)
            for (int iX = 0; iX < numNodes[0]; ++iX)
            {
                const double x = iX * delta[0] + start[0];
                const double y = iY * delta[1] + start[1];
                nodeIds.push_back(s.NodeCreate(Eigen::Vector2d(x, y)));
            }

        break;
    }
    case 3:
    {
        nodeIds.reserve(numNodes[0] * numNodes[1] * numNodes[2]);
        for (int iZ = 0; iZ < numNodes[2]; ++iZ)
            for (int iY = 0; iY < numNodes[1]; ++iY)
                for (int iX = 0; iX < numNodes[0]; ++iX)
                {
                    const double x = iX * delta[0] + start[0];
                    const double y = iY * delta[1] + start[1];
                    const double z = iZ * delta[2] + start[2];
                    nodeIds.push_back(s.NodeCreate(Eigen::Vector3d(x, y, z)));
                }

        break;
    }
    default:
        throw;
    }
    return nodeIds;
}

std::vector<int> CreateNodes(NuTo::Structure& s, std::vector<int> numNodes, std::vector<double> delta)
{
    std::vector<double> start(s.GetDimension(), 0.);
    return CreateNodes(s, numNodes, delta, start);
}

std::pair<int, int> CreateElements(NuTo::Structure& s, std::vector<int> nodeIds, std::vector<int> numDivisions,
                                   NuTo::Interpolation::eShapeType elementShape)
{
    auto info = CreateInterpolationTypeAndGroup(s, elementShape);
    switch (s.GetDimension())
    {
    case 1:
    {
        for (int iX = 0; iX < numDivisions[0]; ++iX)
        {
            std::vector<int> cornerNodes(2);
            cornerNodes[0] = nodeIds[iX];
            cornerNodes[1] = nodeIds[iX + 1];
            s.GroupAddElement(info.first, s.ElementCreate(info.second, cornerNodes));
        }
        break;
    }
    case 2:
    {
        for (int iY = 0; iY < numDivisions[1]; ++iY)
            for (int iX = 0; iX < numDivisions[0]; ++iX)
            {
                std::vector<int> cornerNodes(4);

                cornerNodes[0] = nodeIds[iX + iY * (numDivisions[0] + 1)];
                cornerNodes[1] = nodeIds[iX + 1 + iY * (numDivisions[0] + 1)];
                cornerNodes[2] = nodeIds[iX + 1 + (iY + 1) * (numDivisions[0] + 1)];
                cornerNodes[3] = nodeIds[iX + (iY + 1) * (numDivisions[0] + 1)];

                auto elementNodes = GetElementNodeIds2D(cornerNodes, elementShape);
                for (auto& nodes : elementNodes)
                    s.GroupAddElement(info.first, s.ElementCreate(info.second, nodes));
            }
        break;
    }
    case 3:
    {
        for (int iZ = 0; iZ < numDivisions[2]; ++iZ)
            for (int iY = 0; iY < numDivisions[1]; ++iY)
                for (int iX = 0; iX < numDivisions[0]; ++iX)
                {
                    std::vector<int> cornerNodes(8);
                    int numX = numDivisions[0] + 1;
                    int numY = numDivisions[1] + 1;

                    cornerNodes[0] = nodeIds[iX + iY * numX + iZ * numX * numY];
                    cornerNodes[1] = nodeIds[iX + 1 + iY * numX + iZ * numX * numY];
                    cornerNodes[2] = nodeIds[iX + 1 + (iY + 1) * numX + iZ * numX * numY];
                    cornerNodes[3] = nodeIds[iX + (iY + 1) * numX + iZ * numX * numY];
                    cornerNodes[4] = nodeIds[iX + iY * numX + (iZ + 1) * numX * numY];
                    cornerNodes[5] = nodeIds[iX + 1 + iY * numX + (iZ + 1) * numX * numY];
                    cornerNodes[6] = nodeIds[iX + 1 + (iY + 1) * numX + (iZ + 1) * numX * numY];
                    cornerNodes[7] = nodeIds[iX + (iY + 1) * numX + (iZ + 1) * numX * numY];

                    auto elementNodes = GetElementNodeIds3D(cornerNodes, elementShape);
                    for (auto& nodes : elementNodes)
                        s.GroupAddElement(info.first, s.ElementCreate(info.second, nodes));
                }
        break;
    }
    default:
        throw;
    }
    return info;
}

std::pair<int, int> NuTo::MeshGenerator::Grid(Structure& s, Eigen::VectorXd start, Eigen::VectorXd end,
                                              Eigen::VectorXi numDivisions, Interpolation::eShapeType elementShape)
{
    std::vector<double> startVec(start.data(), start.data() + start.size());
    std::vector<double> endVec(end.data(), end.data() + end.size());
    std::vector<int> numDivisionsVec(numDivisions.data(), numDivisions.data() + numDivisions.size());

    return Grid(s, startVec, endVec, numDivisionsVec, elementShape);
}


std::pair<int, int> NuTo::MeshGenerator::Grid(Structure& s, std::vector<double> start, std::vector<double> end,
                                              std::vector<int> numDivisions, Interpolation::eShapeType elementShape)
{
    const int dimension = s.GetDimension();

    assert(dimension == static_cast<int>(start.size()) and "Dimensions mismatch");
    assert(dimension == static_cast<int>(end.size()) and "Dimensions mismatch");
    assert(dimension == static_cast<int>(numDivisions.size()) and "Dimensions mismatch");

    std::vector<int> numNodes(dimension);
    std::vector<double> delta(dimension);

    for (int i = 0; i < dimension; ++i)
    {
        assert(start[i] < end[i]);

        numNodes[i] = numDivisions[i] + 1;
        delta[i] = (end[i] - start[i]) / numDivisions[i];
    }

    std::vector<int> nodeIds = CreateNodes(s, numNodes, delta, start);
    return CreateElements(s, nodeIds, numDivisions, elementShape);
}

std::pair<int, int> NuTo::MeshGenerator::Grid(Structure& s, std::vector<double> end, std::vector<int> numDivisions,
                                              Interpolation::eShapeType elementShape)
{
    std::vector<double> start(s.GetDimension(), 0.);
    return Grid(s, start, end, numDivisions, elementShape);
}

std::pair<int, int> NuTo::MeshGenerator::Grid(Structure& s, std::vector<double> end, std::vector<int> numDivisions)
{
    std::vector<double> start(s.GetDimension(), 0.);
    return Grid(s, start, end, numDivisions);
}


std::pair<int, int> NuTo::MeshGenerator::Grid(Structure& s, std::vector<double> start, std::vector<double> end,
                                              std::vector<int> numDivisions)
{
    switch (s.GetDimension())
    {
    case 1:
        return MeshGenerator::Grid(s, start, end, numDivisions, eShapeType::TRUSS1D);
    case 2:
        return MeshGenerator::Grid(s, start, end, numDivisions, eShapeType::QUAD2D);
    case 3:
        return MeshGenerator::Grid(s, start, end, numDivisions, eShapeType::BRICK3D);
    default:
        throw;
    }
}

std::function<Eigen::Vector3d(Eigen::Vector3d)> NuTo::MeshGenerator::GetCylinderMapping(double radius, double height)
{
    return [radius, height](Eigen::Vector3d v) -> Eigen::VectorXd {
        v.x() = v.x() * 2 - 1;
        v.y() = v.y() * 2 - 1;
        v.z() = v.z() * 2 - 1;
        v.x() *= 1. + (1. - std::abs(v.x())) / 2.;
        v.y() *= 1. + (1. - std::abs(v.y())) / 2.;
        v.z() *= 1. + (1. - std::abs(v.z())) / 2.;
        Eigen::VectorXd CoordVec(3);
        CoordVec << v.x() * sqrt(1 - (v.y() * v.y()) / 2.0) * radius / 2.0,
                v.y() * sqrt(1 - (v.x() * v.x()) / 2.0) * radius / 2.0, v.z() * height / 2.0;
        return CoordVec;
    };
}
