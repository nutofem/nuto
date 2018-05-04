#pragma once
#include "nuto/base/Group.h"
#include "nuto/base/ValueVector.h"
#include "nuto/mechanics/DirectionEnum.h"
#include "nuto/mechanics/nodes/NodeSimple.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/iga/Nurbs.h"

#include <memory>
#include <vector>

namespace NuTo
{
template <int TDimParameter>
class MeshIga
{
public:
    MeshIga() = default;

    MeshIga(const MeshIga&) = delete;
    MeshIga& operator=(const MeshIga&) = delete;

    MeshIga(MeshIga&&) = default;
    MeshIga& operator=(MeshIga&&) = default;

    // construct mesh with refinement
    void MeshWithRefinement(const DofType& coordinates, const ValueVector<DofType>& dofVector,
                            const std::array<std::vector<double>, TDimParameter>& knots,
                            const std::vector<Eigen::VectorXd>& controlPoints, const std::vector<double>& weights,
                            const std::array<int, TDimParameter>& degree,
                            const std::array<std::vector<double>, TDimParameter>& knotsInserted)
    {
        std::array<std::vector<double>, TDimParameter> rKnots;
        std::vector<Eigen::VectorXd> rControlPoints;
        std::vector<double> rWeights;

        Nurbs<TDimParameter>::Refinement(knots, controlPoints, weights, degree, knotsInserted, rControlPoints, rWeights,
                                         rKnots);

        for (int i = 0; i < rControlPoints.size(); i++)
            std::cout << rWeights[i] << " ::: " << rControlPoints[i].transpose() << std::endl;

        std::vector<NodeSimple*> controlPointsPtrsCorrdinates;
        for (Eigen::VectorXd& coordinate : rControlPoints)
        {
            auto& node = mNodes.Add(coordinate);
            controlPointsPtrsCorrdinates.push_back(&node);
        }

        Nurbs<TDimParameter> nurbsCoordinate(rKnots, controlPointsPtrsCorrdinates, rWeights, degree);
        mDofInterpolations.Insert(coordinates, nurbsCoordinate);

        for (auto& dof : dofVector)
        {
            std::vector<NodeSimple*> controlPointsPtrsDofs;
            for (int i = 0; i < controlPointsPtrsCorrdinates.size(); i++)
            {
                auto& node = mNodes.Add(Eigen::Vector2d({0., 0.}));
                controlPointsPtrsDofs.push_back(&node);
            }

            mDofInterpolations.Insert(dof, Nurbs<TDimParameter>(rKnots, controlPointsPtrsDofs, rWeights, degree));
        }

        int numElementsNurbs = nurbsCoordinate.GetNumElementsTotal();

        for (int i = 0; i < numElementsNurbs; i++)
        {
            std::array<int, TDimParameter> knotIDs = nurbsCoordinate.GetElement(i);
            ElementIga<2> igaCoordinates(knotIDs, mDofInterpolations.At(coordinates));
            ElementCollectionIga<TDimParameter>& element = mElements.Add(igaCoordinates);
            for (auto& dof : dofVector)
            {
                ElementIga<2> igaDispl(knotIDs, mDofInterpolations.At(dof));
                element.AddDofElement(dof, igaDispl);
            }
        }
    }

    Group<NodeSimple> NodesAtAxis(eDirection direction, DofType dofType, double axisOffset = 0., double tol = 1.e-10)
    {
        Group<NodeSimple> group;
        const int directionComponent = ToComponentIndex(direction);
        for (auto& element : this->mElements)
        {
            if (!element.Has(dofType))
                continue;

            auto& dofElement = element.DofElement(dofType);
            auto& coordinateElement = element.CoordinateElement();

            if (coordinateElement.GetNumNodes() != dofElement.GetNumNodes())
                throw Exception(__PRETTY_FUNCTION__,
                                "Coordinate and dof elements must have the same amount of control points");

            Eigen::VectorXd nodeCoordinates = coordinateElement.ExtractNodeValues();
            int dim = coordinateElement.GetDofDimension();

            for (int iNode = 0; iNode < coordinateElement.GetNumNodes(); iNode++)
            {
                Eigen::VectorXd coordinate = nodeCoordinates.segment(iNode * dim, dim);
                if (std::abs(coordinate[directionComponent] - axisOffset) < tol)
                    group.Add(*dofElement.GetNode(iNode));
            }
        }
        return group;
    }

    void CreateBoundary(const DofType& coordinates, const ValueVector<DofType>& dofVector, eDirection dirNotIncluded,
                        int locationBoundary)
    {
        const int dirNotIncludedComponent = ToComponentIndex(dirNotIncluded);
        std::array<eDirection, TDimParameter - 1> directions;
        int count = 0;
        for (int i = 0; i < 3; i++)
        {
            if (dirNotIncludedComponent != i)
            {
                switch (i)
                {
                case 0:
                    directions[count++] = eDirection::X;
                    break;
                case 1:
                    directions[count++] = eDirection::Y;
                    break;
                case 2:
                    directions[count++] = eDirection::Z;
                    break;
                }
            }
            if (count == TDimParameter - 1)
                break;
        }

        // coordinates....
        Nurbs<TDimParameter>& nurbsCoordinates = mDofInterpolations.At(coordinates);

        std::array<std::vector<double>, TDimParameter - 1> boundaryKnotsCoordinates;
        std::vector<NodeSimple*> boundaryControlPointsCoordinates;
        std::vector<double> boundaryWeightsCoordinates;
        std::array<int, TDimParameter - 1> boundaryDegreeCoordinates;

        count = 0;
        for (eDirection dir : directions)
        {
            boundaryKnotsCoordinates[count] = nurbsCoordinates.GetKnotVectorDirection(dir);
            boundaryDegreeCoordinates[count] = nurbsCoordinates.GetDegreeDirection(dir);
            count++;
        }

        nurbsCoordinates.GetControlPointsAnsWeightsBoundary(
                dirNotIncluded, locationBoundary, boundaryControlPointsCoordinates, boundaryWeightsCoordinates);

        Nurbs<TDimParameter - 1> nurbsCoordinatesBoundary(boundaryKnotsCoordinates, boundaryControlPointsCoordinates,
                                                          boundaryWeightsCoordinates, boundaryDegreeCoordinates);

        mDofInterpolationsBoundary.Insert(coordinates, nurbsCoordinatesBoundary);

        int numElementsNurbsBoundary = nurbsCoordinatesBoundary.GetNumElementsTotal();

        // other dofs ....
        for (auto& dof : dofVector)
        {
            Nurbs<TDimParameter>& nurbsDofType = mDofInterpolations.At(dof);

            std::array<std::vector<double>, TDimParameter - 1> boundaryKnots;
            std::array<int, TDimParameter - 1> boundaryDegree;
            std::vector<double> boundaryWeights;
            std::vector<NodeSimple*> boundaryControlPoints;

            int count = 0;
            for (eDirection dir : directions)
            {
                boundaryKnots[count] = nurbsDofType.GetKnotVectorDirection(dir);
                boundaryDegree[count] = nurbsDofType.GetDegreeDirection(dir);
                count++;
            }

            nurbsDofType.GetControlPointsAnsWeightsBoundary(dirNotIncluded, locationBoundary, boundaryControlPoints,
                                                            boundaryWeights);

            Nurbs<TDimParameter - 1> nurbsDofTypeBoundary(boundaryKnots, boundaryControlPoints, boundaryWeights,
                                                          boundaryDegree);

            int numElementsBoundaryDof = nurbsDofTypeBoundary.GetNumElementsTotal();
            assert(numElementsNurbsBoundary == numElementsBoundaryDof);

            mDofInterpolationsBoundary.Insert(dof, nurbsDofTypeBoundary);
        }

        for (int i = 0; i < numElementsNurbsBoundary; i++)
        {
            std::array<int, TDimParameter - 1> knotIDs = mDofInterpolationsBoundary.At(coordinates).GetElement(i);
            ElementIga<TDimParameter - 1> igaCoordinates(knotIDs, mDofInterpolationsBoundary.At(coordinates));
            ElementCollectionIga<TDimParameter - 1>& element = mElementsBoundary.Add(igaCoordinates);
            for (auto& dof : dofVector)
            {
                ElementIga<TDimParameter - 1> igaDispl(knotIDs, mDofInterpolationsBoundary.At(dof));
                element.AddDofElement(dof, igaDispl);
            }
        }
    }

    std::vector<ElementCollectionIga<TDimParameter - 1>*>
    AddElementsBoundaryAtAxis(const DofType& coordinates, const ValueVector<DofType>& dofVector,
                              eDirection dirNotIncluded, double dirNotIncludedValue)
    {
        const int dirNotIncludedComponent = ToComponentIndex(dirNotIncluded);

        Nurbs<TDimParameter - 1> nurbsDofTypeBoundary = mDofInterpolationsBoundary.At(coordinates);
        int numElementsNurbsBoundary = nurbsDofTypeBoundary.GetNumElementsTotal();

        for (int i = 0; i < numElementsNurbsBoundary; i++)
        {
            std::array<int, TDimParameter - 1> knotIDs = nurbsDofTypeBoundary.GetElement(i);
            bool insert = false;
            for (int cp = 0; cp < nurbsDofTypeBoundary.GetNumControlPointsElement(); cp++)
            {
                const NodeSimple* node = GetControlPointElement(knotIDs, cp);
                const Eigen::VectorXd& nodeCoord = node->GetValues();
                if (nodeCoord(dirNotIncludedComponent) != dirNotIncludedValue)
                    break;
            }

            if (insert)
            {
                ElementIga<TDimParameter - 1> igaCoordinates(knotIDs, mDofInterpolationsBoundary.At(coordinates));
                ElementCollectionIga<TDimParameter - 1>& element = mElementsBoundary.Add(igaCoordinates);
                for (auto& dof : dofVector)
                {
                    ElementIga<TDimParameter - 1> igaDispl(knotIDs, mDofInterpolationsBoundary.At(dof));
                    element.AddDofElement(dof, igaDispl);
                }
            }
        }
    }


    Group<NodeSimple> NodesTotal()
    {
        throw Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
    }

    Group<NodeSimple> NodesTotal(DofType d)
    {
        throw Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
    }

    Group<ElementCollectionIga<TDimParameter>> ElementsTotal()
    {
        throw Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
    }

public:
    // the structure itself
    //! @brief For each dof a nurbs is stored (its the interpolation)
    DofContainer<Nurbs<TDimParameter>> mDofInterpolations;

    //! @brief Isogeometric elements
    ValueVector<ElementCollectionIga<TDimParameter>> mElements;

    // the boundary
    //! @brief For each dof a nurbs is stored (its the interpolation)
    DofContainer<Nurbs<TDimParameter - 1>> mDofInterpolationsBoundary;

    //! @brief Isogeometric elements
    ValueVector<ElementCollectionIga<TDimParameter - 1>> mElementsBoundary;

    // nodes or control points
    //! @brief the control points
    ValueVector<NodeSimple> mNodes;
};
} /* NuTo */
