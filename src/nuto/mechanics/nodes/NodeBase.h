#ifndef NODEBASE_H
#define NODEBASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif  // ENABLE_SERIALIZATION
#include <vector>

#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
template <class T> class FullMatrix;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all nodes
class NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    NodeBase();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)=0;

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues) = 0;

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const = 0;

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering) = 0;

    //! @brief returns the number of coordinates of the node
    //! @return number of coordinates
    virtual int GetNumCoordinates()const;

    //! @brief returns the number of displacements of the node
    //! @return number of displacements
    virtual int GetNumDisplacements()const;

    //! @brief returns the number of rotations of the node
    //! @return number of rotations
    virtual int GetNumRotations()const;

    //! @brief returns the number of temperatures of the node
    //! @return number of temperatures
    virtual int GetNumTemperatures()const;

protected:
    //the base class of the nodes must not contain any data

};
}//namespace NuTo
#endif //NODEBASE_H

