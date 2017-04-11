#pragma once
#include <map>
#include <iosfwd>
#include "math/SparseMatrix.h"
#include "mechanics/constraints/Equation.h"

namespace NuTo
{
namespace Constraint
{
//! @brief stores constraint equations, separated by their dof type
class Constraints
{
    typedef std::vector<Equation> Equations;

public:
    //! @brief adds an equation
    //! @param dof dof type
    //! @param equation constraint equation
    void Add(Node::eDof dof, Equation equation);

    //! @brief adds multiple equations
    //! @param dof dof type
    //! @param equations constraint equations
    void Add(Node::eDof dof, std::vector<Equation> equations);

    //! @brief builds the time dependent constraint rhs vector for a specific dof type
    //! @param dof dof type
    //! @param time global time
    //! @return constraint rhs vector
    Eigen::VectorXd GetRhs(Node::eDof dof, double time) const;

    //! @brief builds a sparse matrix containing the constraint terms for a specific dof type
    //! @param dof dof type
    //! @param nDofs total number of dofs for the dof type dof, required for a proper resize
    //!        of the sparse matrix
    //! @return sparse matrix containing the constraint terms
    SparseMatrixCSRVector2General<double> BuildConstraintMatrix(Node::eDof dof, int nDofs) const;

    //! @brief calculates the number of constraint equations for a specific dof type
    //! @param dof dof type
    //! @return number of constraint equations
    int GetNumEquations(Node::eDof dof) const;

    //! @brief exchanges the constrainted nodes in all terms
    //! @param oldNode node that is replaced
    //! @param newNode new node
    void ExchangeNodePtr(const NodeBase& oldNode, const NodeBase& newNode);

    //! @brief prints all equations
    friend std::ostream& operator<<(std::ostream& out, const Constraints& constraints);

    //! @brief getter for mHasNewConstraints
    //! @return true if new constraints were added
    bool HasNewConstraints() const;

    //! @brief setter for mHasNewConstraints
    void SetHasNewConstraints(bool value);

private:
    //! @brief dof-wise storage of constraint equatiosn
    std::map<Node::eDof, Equations> mEquations;

    //! @brief flag that indiciates whether or not new constraints were added
    //! @remark the global dof numbering needs a rebuild, if this flag is true.
    //! And it has to be set to false after rebuilding.
    //! Problems with this approach: The behavior of this class is not influenced
    //! by this state variable, currently only NuTo::Assembler is. Thus, it should
    //! be solved in a different approach. Maybe:
    //!  - Assembler keeps track of NumConstraints
    //!  - Assembler keeps track of std::hash(mEquations)
    //!  - Not allowed to add constraints after building the global dof numbering
    //!       = mConstraints is a const member at mAssembler
    bool mHasNewConstraints = false;
};

} /* Constaint */
} /* NuTo */
