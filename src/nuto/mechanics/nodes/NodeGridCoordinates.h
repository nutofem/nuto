#ifndef NODEGRIDCOORDINATES_H
#define NODEGRIDCOORDINATES_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeGridCoordinatesBase.h"

namespace NuTo
{
//! @author Andrea Keßler, ISM
//! @date March 2010
//! @brief ... standard class for nodes without coordinates, grid nodes, ID given
class NodeGridCoordinates : public NodeGridCoordinatesBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeGridCoordinates() : NodeGridCoordinatesBase ()
    {
        mNodeID = -1;
    }

    //! @brief constructor
    NodeGridCoordinates (const unsigned int rNodeID)  : NodeGridCoordinatesBase ()
    {
        mNodeID = rNodeID;
   }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
        & BOOST_SERIALIZATION_NVP(mNodeID);
    }
#endif // ENABLE_SERIALIZATION


    //! @brief returns the number of coordinates of the node
    //! @return number of coordinates
    int GetNumCoordinates()const
    {
        return 3;
        throw MechanicsException("[NuTo::NodeGridCoordinates::GetNumCoordinates] not implemented.");
    }

    //! @brief get the ID
    //! @param rNodeNumber number of the node
    int GetNodeID() const
    {
        return 4;
    }

    //! @brief writes the coordinates of a node to the prescribed pointer
    //! @param rCoordinates coordinates
    void GetCoordinates()const
    {
        throw MechanicsException("[NuTo::NodeGridCoordinates::GetCoordinates] not implemented.");
    }


    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        //empty since coordinates are no DOFs
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        //empty since coordinates are no DOFs
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        //empty since coordinates are no DOFs
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        //empty since coordinates are no DOFs
    }

protected:
    int mNodeID;
};
}//namespace NuTo
#endif //NODE_GRIDCOORDINATES_H
