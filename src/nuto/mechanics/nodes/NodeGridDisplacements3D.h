// $Id$
#ifndef NodeGridDisplacements_3d_H
#define NodeGridDisplacements_3d_H

#include "nuto/mechanics/nodes/NodeGrid3D.h"
#include "nuto/mechanics/nodes/NodeDof.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for grid nodes with displacements
class NodeGridDisplacements3D : public NodeGrid3D, public NodeDof<0,3,0,0,0,0,0>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeGridDisplacements3D() : NodeGrid3D (), NodeDof<0,3,0,0,0,0,0>()
    {}

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    NodeGridDisplacements3D* Clone()const
    {
    	throw MechanicsException("[NodeGridDisplacements3D::Clone] to be implemented.");
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeGrid3D)
           & BOOST_SERIALIZATION_BASE_OBJECT_NVP(BOOST_IDENTITY_TYPE((NodeDof<0,3,0,0,0,0>)));
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        NodeDof<0,3,0,0,0,0,0>::SetGlobalDofs(rDOF);
    }

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(int rTimeDerivative, const FullVector<double,Eigen::Dynamic>& rActiveDofValues, const FullVector<double,Eigen::Dynamic>& rDependentDofValues)
    {
    	NodeDof<0,3,0,0,0,0,0>::SetGlobalDofValues(rTimeDerivative, rActiveDofValues, rDependentDofValues);
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(int rTimeDerivative, FullVector<double,Eigen::Dynamic>& rActiveDofValues, FullVector<double,Eigen::Dynamic>& rDependentDofValues) const
    {
        NodeDof<0,3,0,0,0,0,0>::GetGlobalDofValues(rTimeDerivative, rActiveDofValues, rDependentDofValues);
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
    	NodeDof<0,3,0,0,0,0,0>::RenumberGlobalDofs(rMappingInitialToNewOrdering);
    }

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const
    {
    	return std::string("NodeGridDisplacements3D");
    }
//raus, weil steht in NodeDisplacements3D
//protected:
  //  int mDof[3];

    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    virtual Node::eNodeType GetNodeType()const
    {
        return Node::NodeGridDisplacements3D;
    }

};
}//namespace NuTo

#endif //NodeGridDisplacements_3d_H
