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



template<>
std::pair<int, int> MeshGenerator::CreateGrid<1>(Structure& rS,
                                              Eigen::Matrix<double, 1, 1> rEnd,
                                              Eigen::Matrix<int, 1, 1> rNumDivisions,
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
std::pair<int, int> MeshGenerator::CreateGrid<2>(Structure& rS,
                                              Eigen::Matrix<double, 2, 1> rEnd,
                                              Eigen::Matrix<int, 2, 1> rNumDivisions,
                                              Interpolation::eShapeType rElementShape)
{
    Eigen::Vector2i numNodes = rNumDivisions.array() + 1;
    Eigen::Vector2d delta = rEnd.cwiseQuotient(rNumDivisions.cast<double>());

    std::vector<int> nodeIds;
    nodeIds.reserve(numNodes.prod());


    for (int iY = 0; iY < numNodes.y(); ++iY)
        for (int iX = 0; iX < numNodes.x(); ++iX)
            nodeIds.push_back(rS.NodeCreate(Eigen::Vector2d(iX * delta.x(), iY * delta.y())));

    auto info = CreateInterpolationTypeAndGroup(rS, rElementShape);

    for (int iY = 0; iY < rNumDivisions.y(); ++iY)
        for (int iX = 0; iX < rNumDivisions.x(); ++iX)
        {
            std::vector<int> cornerNodes(4);

            cornerNodes[0] = nodeIds[iX   +  iY    * numNodes.x()];
            cornerNodes[1] = nodeIds[iX+1 +  iY    * numNodes.x()];
            cornerNodes[2] = nodeIds[iX+1 + (iY+1) * numNodes.x()];
            cornerNodes[3] = nodeIds[iX   + (iY+1) * numNodes.x()];

            auto elementNodes = GetElementNodeIds2D(cornerNodes, rElementShape);
            for (auto& nodes : elementNodes)
                rS.GroupAddElement(info.first, rS.ElementCreate(info.second, nodes));
        }
    return info;
}

template<>
std::pair<int, int> MeshGenerator::CreateGrid<3>(Structure& rS,
                                              Eigen::Matrix<double, 3, 1> rEnd,
                                              Eigen::Matrix<int, 3, 1> rNumDivisions,
                                              eShapeType rElementShape)
{
    Eigen::Vector3i numNodes = rNumDivisions.array() + 1;
    Eigen::Vector3d delta = rEnd.cwiseQuotient(rNumDivisions.cast<double>());

    std::vector<int> nodeIds;
    nodeIds.reserve(numNodes.prod());


    for (int iZ = 0; iZ < numNodes.z(); ++iZ)
        for (int iY = 0; iY < numNodes.y(); ++iY)
            for (int iX = 0; iX < numNodes.x(); ++iX)
                nodeIds.push_back(rS.NodeCreate(Eigen::Vector3d(iX * delta.x(), iY * delta.y(), iZ * delta.z())));

    auto info = CreateInterpolationTypeAndGroup(rS, rElementShape);

    for (int iZ = 0; iZ < rNumDivisions.z(); ++iZ)
        for (int iY = 0; iY < rNumDivisions.y(); ++iY)
            for (int iX = 0; iX < rNumDivisions.x(); ++iX)
            {
                std::vector<int> cornerNodes(8);

                cornerNodes[0] = nodeIds[iX   +  iY    * numNodes.x() +  iZ    * numNodes.x() * numNodes.y()];
                cornerNodes[1] = nodeIds[iX+1 +  iY    * numNodes.x() +  iZ    * numNodes.x() * numNodes.y()];
                cornerNodes[2] = nodeIds[iX+1 + (iY+1) * numNodes.x() +  iZ    * numNodes.x() * numNodes.y()];
                cornerNodes[3] = nodeIds[iX   + (iY+1) * numNodes.x() +  iZ    * numNodes.x() * numNodes.y()];
                cornerNodes[4] = nodeIds[iX   +  iY    * numNodes.x() + (iZ+1) * numNodes.x() * numNodes.y()];
                cornerNodes[5] = nodeIds[iX+1 +  iY    * numNodes.x() + (iZ+1) * numNodes.x() * numNodes.y()];
                cornerNodes[6] = nodeIds[iX+1 + (iY+1) * numNodes.x() + (iZ+1) * numNodes.x() * numNodes.y()];
                cornerNodes[7] = nodeIds[iX   + (iY+1) * numNodes.x() + (iZ+1) * numNodes.x() * numNodes.y()];

                auto elementNodes = GetElementNodeIds3D(cornerNodes, rElementShape);
                for (auto& nodes : elementNodes)
                    rS.GroupAddElement(info.first, rS.ElementCreate(info.second, nodes));
            }

    return info;
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
