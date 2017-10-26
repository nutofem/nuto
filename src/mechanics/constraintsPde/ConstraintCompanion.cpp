#include "mechanics/constraintsPde/ConstraintCompanion.h"


namespace NuTo
{
namespace ConstraintPde
{

std::vector<Equation> Component(const NodeSimple& node, std::vector<eDirection> directions, double value)
{
    return Component(node, directions, RhsConstant(value));
}


std::vector<Equation> Component(const NodeSimple& node, std::vector<eDirection> directions, RhsFunction rhs)
{
    std::vector<Equation> eqs;
    for (auto direction : directions)
    {
        int component = ToComponentIndex(direction);
        eqs.push_back(Equation(node, component, rhs));
    }
    return eqs;
}


std::vector<Equation> Component(const Groups::Group<NodeSimple>& nodes, std::vector<eDirection> directions,
                                double value)
{
    return Component(nodes, directions, RhsConstant(value));
}


std::vector<Equation> Component(const Groups::Group<NodeSimple>& nodes, std::vector<eDirection> directions,
                                RhsFunction rhs)
{
    std::vector<Equation> eqs;
    for (const auto& node : nodes)
    {
        auto tmpEqs = Component(node, directions, rhs);
        eqs.insert(eqs.end(), tmpEqs.begin(), tmpEqs.end());
    }
    return eqs;
}


// Equation Direction(const NodeSimple& node, Eigen::VectorXd direction, RhsFunction rhs)
//{
//    if (node.GetNumDofs() < direction.rows())
//        // TODO This check is not meaningful at the moment, since this method returns the
//        // total number of dofs, say 4 (3 disp, 1 temp). This will not find the error
//        // if you try to constrain the Z component of the temperature. Which would be wrong.
//        throw Exception(__PRETTY_FUNCTION__, "Dimension mismatch");

//    direction.normalize();
//    Equation e(rhs);
//    for (int iComponent = 0; iComponent < direction.rows(); ++iComponent)
//        e.AddTerm(Term(node, iComponent, direction[iComponent]));

//    return e;
//}

// Equation Direction(const NodeSimple& node, Eigen::VectorXd direction, double value)
//{
//    return Direction(node, direction, RhsConstant(value));
//}

// std::vector<Equation> Direction(const Groups::Group<NodeSimple>& nodes, Eigen::VectorXd direction, RhsFunction rhs)
//{
//    std::vector<Equation> eqs;
//    for (auto& nodePair : nodes)
//        eqs.push_back(Direction(*nodePair.second, direction, rhs));
//    return eqs;
//}

// std::vector<Equation> Direction(constGroups::Group<NodeSimple>& nodes, Eigen::VectorXd direction, double value)
//{
//    return Direction(nodes, direction, RhsConstant(value));
//}

// Equation Value(const NodeSimple& node, double value)
//{
//    return Value(node, [=](double) { return value; });
//}

// Equation Value(const NodeSimple& node, RhsFunction rhs)
//{
//    return Direction(node, Eigen::VectorXd::Ones(1), rhs);
//}

// std::vector<Equation> Value(constGroups::Group<NodeSimple>& nodes, double value)
//{
//    return Direction(nodes, Eigen::VectorXd::Ones(1), RhsConstant(value));
//}

// std::vector<Equation> Value(constGroups::Group<NodeSimple>& nodes, RhsFunction rhs)
//{
//    return Direction(nodes, Eigen::VectorXd::Ones(1), rhs);
//}

} /* ConstraintPde */
} /* NuTo */
