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

    void Refinement(DofType d, const std::array<std::vector<double>, TDimParameter>& knots,
                    const std::vector<Eigen::VectorXd>& controlPoints, const std::vector<double>& weights,
                    const std::array<int, TDimParameter>& degree, const std::array<int, TDimParameter>& refinementLevel)
    {
        std::array<std::vector<double>, TDimParameter> rKnots;
        std::vector<Eigen::VectorXd> rControlPoints;
        std::vector<double> rWeights;
        std::array<int, TDimParameter> rDegree;

        Nurbs<TDimParameter>::Refinement(knots, controlPoints, weights, degree, refinementLevel, rKnots, rControlPoints,
                                         rWeights, rDegree);

        std::vector<NodeSimple*> controlPointsPtrs;
        for (Eigen::VectorXd& coordinate : rControlPoints)
        {
            auto& node = mNodes.Add(coordinate);
            controlPointsPtrs.push_back(&node);
        }
        mDofInterpolations.Insert(d, Nurbs<TDimParameter>(rKnots, controlPointsPtrs, rWeights, rDegree));
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
    //! @brief For each dof a nurbs is stored (its the interpolation)
    DofContainer<Nurbs<TDimParameter>> mDofInterpolations;

    //! @brief Isogeometric elements
    ValueVector<ElementCollectionIga<TDimParameter>> mElements;

    //! @brief the control points
    ValueVector<NodeSimple> mNodes;
};
} /* NuTo */
