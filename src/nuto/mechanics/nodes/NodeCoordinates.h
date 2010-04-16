#ifndef NODE_COORDINATES_H
#define NODE_COORDINATES_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeCoordinatesBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for reference nodes
template<int NUMCOORDINATES>
class NodeCoordinates : public NodeCoordinatesBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief constructor
    NodeCoordinates() : NodeCoordinatesBase ()
    {
        for (int count=0; count<NUMCOORDINATES; count++)
            mCoordinates[count]=0;
    }

    //! @brief constructor
    NodeCoordinates(const double rCoordinates[NUMCOORDINATES])  : NodeCoordinatesBase ()
    {
        memcpy(mCoordinates,rCoordinates,NUMCOORDINATES*sizeof(double));
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
        & BOOST_SERIALIZATION_NVP(mCoordinates);
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the number of coordinates of the node
    //! @return number of coordinates
    int GetNumCoordinates()const
    {
        return NUMCOORDINATES;
    }

    //! @brief set the coordinates
    //! @param rCoordinates  given coordinates
    void SetCoordinates(const double rCoordinates[NUMCOORDINATES])
    {
        memcpy(mCoordinates,rCoordinates,NUMCOORDINATES*sizeof(double));
    }

    //! @brief writes the coordinates of a node to the prescribed pointer
    //! @param rCoordinates coordinates
    void GetCoordinates(double rCoordinates[NUMCOORDINATES])const
    {
        memcpy(rCoordinates,mCoordinates,NUMCOORDINATES*sizeof(double));
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
    double mCoordinates[NUMCOORDINATES];

};
}//namespace NuTo
#endif //NODE_COORDINATES_H
