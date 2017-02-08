#include "mechanics/MechanicsEnums.h"
#include "mechanics/mesh/MeshGenerator.h"

#include "mechanics/structures/unstructured/Structure.h"

using namespace NuTo::Interpolation;

std::pair<int, int> CreateInterpolationTypeAndGroup(NuTo::Structure& rS, eShapeType rElementShape)
{
    int interpolationType = rS.InterpolationTypeCreate(rElementShape);
    rS.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    int elementGroup = rS.GroupCreate(NuTo::eGroupId::Elements);
    return std::make_pair(elementGroup, interpolationType);
}

std::vector<std::vector<int>> GetElementNodeIds3D(const std::vector<int>& rCornerNodes, eShapeType rElementShape)
{
    switch (rElementShape)
    {
    case eShapeType::BRICK3D:
    {
        return {rCornerNodes};
    }
    case eShapeType::TETRAHEDRON3D:
    {
        std::vector<int> nodesTet0({rCornerNodes[0], rCornerNodes[1], rCornerNodes[3], rCornerNodes[7]});
        std::vector<int> nodesTet1({rCornerNodes[0], rCornerNodes[1], rCornerNodes[7], rCornerNodes[4]});
        std::vector<int> nodesTet2({rCornerNodes[5], rCornerNodes[4], rCornerNodes[7], rCornerNodes[1]});
        std::vector<int> nodesTet3({rCornerNodes[6], rCornerNodes[5], rCornerNodes[7], rCornerNodes[1]});
        std::vector<int> nodesTet4({rCornerNodes[2], rCornerNodes[7], rCornerNodes[1], rCornerNodes[6]});
        std::vector<int> nodesTet5({rCornerNodes[2], rCornerNodes[3], rCornerNodes[1], rCornerNodes[7]});
        return {nodesTet0, nodesTet1, nodesTet2, nodesTet3, nodesTet4, nodesTet5};
    }
    case eShapeType::PRISM3D:
    {
        std::vector<int> nodes0({rCornerNodes[0], rCornerNodes[1], rCornerNodes[2],
                                 rCornerNodes[4], rCornerNodes[5], rCornerNodes[6]});
        std::vector<int> nodes1({rCornerNodes[0], rCornerNodes[2], rCornerNodes[3],
                                 rCornerNodes[4], rCornerNodes[6], rCornerNodes[7]});
        return {nodes0, nodes1};
    }
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, ShapeTypeToString(rElementShape) + " not supported as 3D element");
    }
}

std::vector<std::vector<int>> GetElementNodeIds2D(const std::vector<int>& rCornerNodes, eShapeType rElementShape)
{
    switch (rElementShape)
    {
    case eShapeType::QUAD2D:
    {
        return {rCornerNodes};
    }
    case eShapeType::TRIANGLE2D:
    {
        std::vector<int> e1{rCornerNodes[0], rCornerNodes[1], rCornerNodes[2]};
        std::vector<int> e2{rCornerNodes[0], rCornerNodes[2], rCornerNodes[3]};
        return {e1, e2};
    }
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, ShapeTypeToString(rElementShape) + " not supported as 3D element");
    }
}


namespace NuTo
{



template<>
std::pair<int, int> MeshGenerator::Grid<1>(Structure& rS,
                                           std::array<double, 1> rEnd,
                                           std::array<int, 1> rNumDivisions,
                                           eShapeType rElementShape)
{
    int numNodes = rNumDivisions[0] + 1;
    double delta = rEnd[0] / rNumDivisions[0];

    std::vector<int> nodeIds;
    nodeIds.reserve(numNodes);

    for (int iX = 0; iX < numNodes; ++iX)
        nodeIds.push_back(rS.NodeCreate(Eigen::Matrix<double, 1, 1>::Constant(iX * delta)));

    auto info = CreateInterpolationTypeAndGroup(rS, rElementShape);

    for (int iX = 0; iX < rNumDivisions[0]; ++iX)
    {
        std::vector<int> cornerNodes(2);
        cornerNodes[0] = nodeIds[iX  ];
        cornerNodes[1] = nodeIds[iX+1];

        rS.GroupAddElement(info.first, rS.ElementCreate(info.second, cornerNodes));
    }
    return info;
}

template<>
std::pair<int, int> MeshGenerator::Grid<2>(Structure& rS,
                                           std::array<double, 2> rEnd,
                                           std::array<int, 2> rNumDivisions,
                                           Interpolation::eShapeType rElementShape)
{
    std::array<int, 2> numNodes;
    std::array<double, 2> delta;

    for (int i = 0; i < 2; ++i)
    {
        numNodes[i] = rNumDivisions[i] + 1;
        delta[i] = rEnd[i] / rNumDivisions[i];
    }

    std::vector<int> nodeIds;
    nodeIds.reserve(numNodes[0] * numNodes[1]);

    for (int iY = 0; iY < numNodes[1]; ++iY)
        for (int iX = 0; iX < numNodes[0]; ++iX)
            nodeIds.push_back(rS.NodeCreate(Eigen::Vector2d(iX * delta[0], iY * delta[1])));

    auto info = CreateInterpolationTypeAndGroup(rS, rElementShape);

    for (int iY = 0; iY < rNumDivisions[1]; ++iY)
        for (int iX = 0; iX < rNumDivisions[0]; ++iX)
        {
            std::vector<int> cornerNodes(4);

            cornerNodes[0] = nodeIds[iX   +  iY    * numNodes[0]];
            cornerNodes[1] = nodeIds[iX+1 +  iY    * numNodes[0]];
            cornerNodes[2] = nodeIds[iX+1 + (iY+1) * numNodes[0]];
            cornerNodes[3] = nodeIds[iX   + (iY+1) * numNodes[0]];

            auto elementNodes = GetElementNodeIds2D(cornerNodes, rElementShape);
            for (auto& nodes : elementNodes)
                rS.GroupAddElement(info.first, rS.ElementCreate(info.second, nodes));
        }
    return info;
}

template<>
std::pair<int, int> MeshGenerator::Grid<3>(Structure& rS,
                                           std::array<double, 3> rEnd,
                                           std::array<int, 3> rNumDivisions,
                                           eShapeType rElementShape)
{
    std::array<int, 3> numNodes;
    std::array<double, 3> delta;

    for (int i = 0; i < 3; ++i)
    {
        numNodes[i] = rNumDivisions[i] + 1;
        delta[i] = rEnd[i] / rNumDivisions[i];
    }

    std::vector<int> nodeIds;
    nodeIds.reserve(numNodes[0] * numNodes[1] * numNodes[2]);

    for (int iZ = 0; iZ < numNodes[2]; ++iZ)
        for (int iY = 0; iY < numNodes[1]; ++iY)
            for (int iX = 0; iX < numNodes[0]; ++iX)
                nodeIds.push_back(rS.NodeCreate(Eigen::Vector3d(iX * delta[0], iY * delta[1], iZ * delta[2])));

    auto info = CreateInterpolationTypeAndGroup(rS, rElementShape);

    for (int iZ = 0; iZ < rNumDivisions[2]; ++iZ)
        for (int iY = 0; iY < rNumDivisions[1]; ++iY)
            for (int iX = 0; iX < rNumDivisions[0]; ++iX)
            {
                std::vector<int> cornerNodes(8);

                cornerNodes[0] = nodeIds[iX   +  iY    * numNodes[0] +  iZ    * numNodes[0] * numNodes[1]];
                cornerNodes[1] = nodeIds[iX+1 +  iY    * numNodes[0] +  iZ    * numNodes[0] * numNodes[1]];
                cornerNodes[2] = nodeIds[iX+1 + (iY+1) * numNodes[0] +  iZ    * numNodes[0] * numNodes[1]];
                cornerNodes[3] = nodeIds[iX   + (iY+1) * numNodes[0] +  iZ    * numNodes[0] * numNodes[1]];
                cornerNodes[4] = nodeIds[iX   +  iY    * numNodes[0] + (iZ+1) * numNodes[0] * numNodes[1]];
                cornerNodes[5] = nodeIds[iX+1 +  iY    * numNodes[0] + (iZ+1) * numNodes[0] * numNodes[1]];
                cornerNodes[6] = nodeIds[iX+1 + (iY+1) * numNodes[0] + (iZ+1) * numNodes[0] * numNodes[1]];
                cornerNodes[7] = nodeIds[iX   + (iY+1) * numNodes[0] + (iZ+1) * numNodes[0] * numNodes[1]];

                auto elementNodes = GetElementNodeIds3D(cornerNodes, rElementShape);
                for (auto& nodes : elementNodes)
                    rS.GroupAddElement(info.first, rS.ElementCreate(info.second, nodes));
            }

    return info;
}

template<>
std::pair<int, int> MeshGenerator::Grid<1>(Structure& rS, std::array<double, 1> rEnd, std::array<int, 1> rNumDivisions)
{
    return MeshGenerator::Grid<1>(rS, rEnd, rNumDivisions, eShapeType::TRUSS1D);
}
template<>
std::pair<int, int> MeshGenerator::Grid<2>(Structure& rS, std::array<double, 2> rEnd, std::array<int, 2> rNumDivisions)
{
    return MeshGenerator::Grid<2>(rS, rEnd, rNumDivisions, eShapeType::QUAD2D);
}
template<>
std::pair<int, int> MeshGenerator::Grid<3>(Structure& rS, std::array<double, 3> rEnd, std::array<int, 3> rNumDivisions)
{
    return MeshGenerator::Grid<3>(rS, rEnd, rNumDivisions, eShapeType::BRICK3D);
}


} // namespace NuTo

std::function<Eigen::Vector3d(Eigen::Vector3d)> NuTo::MeshGenerator::GetCylinderMapping(
    double rRadius,
    double rHeight)
{
    return [&rRadius,rHeight](Eigen::Vector3d v) -> Eigen::VectorXd
                {
                    v.x() = v.x() * 2 -1;
                    v.y() = v.y() * 2 -1;
                    v.z() = v.z() * 2 -1;
                    v.x() *= 1. + (1. - std::abs(v.x())) / 2.;
                    v.y() *= 1. + (1. - std::abs(v.y())) / 2.;
                    v.z() *= 1. + (1. - std::abs(v.z())) / 2.;
                    Eigen::VectorXd CoordVec(3);
                    CoordVec << v.x() * sqrt(1 - (v.y() * v.y()) / 2.0 ) * rRadius / 2.0,
                                v.y() * sqrt(1 - (v.x() * v.x()) / 2.0 ) * rRadius / 2.0,
                                v.z() * rHeight / 2.0;
                    return CoordVec;
                };
}
