#pragma once
#include "math/SparseMatrix.h"
#include "mechanics/constraints/Equation.h"
#include <map>

namespace NuTo
{
namespace Constraint
{
class Constraints
{
    using Equations = std::vector<Equation>;

public:
    void Add(Node::eDof dof, Equation equation);

    void Add(Node::eDof dof, std::vector<Equation> equations);

    void SetHasNewConstraints(bool value);

    bool HasNewConstraints() const;

    Eigen::VectorXd GetRhs(Node::eDof dof, double time) const;

    void BuildConstraintMatrix(SparseMatrix<double>& rConstraintMatrix, Node::eDof dof) const;

    int GetNumEquations(Node::eDof dof) const;    

    void ExchangeNodePtr(const NodeBase& oldNode, const NodeBase& newNode);
   

private:
    std::map<Node::eDof, Equations> mEquations;
    bool mHasNewConstraints = false;
};

} /* Constaint */
} /* NuTo */
