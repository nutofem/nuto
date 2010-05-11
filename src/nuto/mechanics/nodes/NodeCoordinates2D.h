#ifndef NODE_COORDINATES_2D_H
#define NODE_COORDINATES_2D_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for 2D reference nodes
class NodeCoordinates2D : public virtual NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief constructor
    NodeCoordinates2D();

    //! @brief constructor
    NodeCoordinates2D(const double rCoordinates[2]);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the number of coordinates of the node
    //! @return number of coordinates
    virtual int GetNumCoordinates()const;

    //! @brief set the coordinates
    //! @param rCoordinates  given coordinates
    virtual void SetCoordinates2D(const double rCoordinates[1]);

    //! @brief writes the coordinates of a node to the prescribed pointer
    //! @param rCoordinates coordinates
    virtual void GetCoordinates2D(double rCoordinates[1])const;

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
    virtual std::string GetNodeTypeStr()const;

protected:
    double mCoordinates[2];

};
}//namespace NuTo
#endif //NODE_COORDINATES_2D_H
