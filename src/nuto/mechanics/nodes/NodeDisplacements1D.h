#ifndef NODE_DISPLACEMENTS_1D_H
#define NODE_DISPLACEMENTS_1D_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for nodes having displacement degrees of freedom
class NodeDisplacements1D : public virtual NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeDisplacements1D();

    //! @brief constructor
    NodeDisplacements1D (const double rDisplacements[1]);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief returns the number of displacements of the node
    //! @return number of displacements
    int GetNumDisplacements()const;

    //! @brief gives the global DOF of a displacement component
    //! @param rComponent component
    //! @return global DOF
    int GetDofDisplacement(int rComponent)const;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    void SetDisplacements1D(const double rDisplacements[1]);

    //! @brief writes the displacements of a node to the prescribed pointer
    //! @param rDisplacements displacements
    void GetDisplacements1D(double rDisplacements[1])const;

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF);

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues);

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const;

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering);

    //! @brief returns the type of the node
    //! @return type
    std::string GetNodeTypeStr()const;

protected:
    double mDisplacements[1];
    int mDOF[1];

};
}//namespace NuTo
#endif //NODE_DISPLACEMENTS_1D_H
