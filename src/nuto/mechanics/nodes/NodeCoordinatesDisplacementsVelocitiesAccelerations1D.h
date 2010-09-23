#ifndef NODE_COORDINATES_DISPLACEMENTS_VELOCITIES_ACCELERATIONS_1D_H_
#define NODE_COORDINATES_DISPLACEMENTS_VELOCITIES_ACCELERATIONS_1D_H_

#include "nuto/mechanics/nodes/NodeCoordinates1D.h"
#include "nuto/mechanics/nodes/NodeDisplacements1D.h"
#include "nuto/mechanics/nodes/NodeVelocities1D.h"
#include "nuto/mechanics/nodes/NodeAccelerations1D.h"

namespace NuTo
{
//! @author Stefan Eckardt, IFF
//! @date July 2010
//! @brief ... class for nodes having coordinates, displacements, velocities and accelerations
class NodeCoordinatesDisplacementsVelocitiesAccelerations1D : public virtual NuTo::NodeCoordinates1D,
    public virtual NuTo::NodeDisplacements1D, public NuTo::NodeVelocities1D, public NuTo::NodeAccelerations1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    NodeCoordinatesDisplacementsVelocitiesAccelerations1D();

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

    //! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues);

    //! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const;

    //! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues);

    //! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const;

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering);

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeCoordinates1D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeDisplacements1D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeVelocities1D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeAccelerations1D);
    }
#endif  // ENABLE_SERIALIZATION
};

}
#endif // NODE_COORDINATES_DISPLACEMENTS_VELOCITIES_ACCELERATIONS_1D_H_
