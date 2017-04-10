#pragma once
#include <map>
#include <iosfwd>
#include "math/SparseMatrix.h"
#include "mechanics/constraints/Equation.h"

namespace NuTo
{
namespace Constraint
{
class Constraints
{
    typedef std::vector<Equation> Equations;

public:
    void Add(Node::eDof dof, Equation equation);

    void Add(Node::eDof dof, std::vector<Equation> equations);

    void SetHasNewConstraints(bool value);

    bool HasNewConstraints() const;

    Eigen::VectorXd GetRhs(Node::eDof dof, double time) const;

    void BuildConstraintMatrix(SparseMatrix<double>& rConstraintMatrix, Node::eDof dof) const;

    int GetNumEquations(Node::eDof dof) const;    

    void ExchangeNodePtr(const NodeBase& oldNode, const NodeBase& newNode);

    friend std::ostream& operator<<(std::ostream& out, const Constraints& constraints);

private:
    std::map<Node::eDof, Equations> mEquations;
    bool mHasNewConstraints = false;
};

} /* Constaint */
} /* NuTo */
