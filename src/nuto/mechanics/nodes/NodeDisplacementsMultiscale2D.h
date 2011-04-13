// $Id: NodeDisplacementsMultiscale2D.h 331 2010-10-06 09:32:11Z arnold2 $
#ifndef NodeDisplacementsMultiscale2D_H
#define NodeDisplacementsMultiscale2D_H

#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
class StructureMultiscale;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for nodes having displacements in a multiscale model (coarse and fine scale displacements)
class NodeDisplacementsMultiscale2D : virtual public NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeDisplacementsMultiscale2D(NuTo::StructureMultiscale* rStructureMultiscale, bool rCrackedDomain);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION
    //! @brief returns the number of displacements of the node
    //! @return number of displacements
    int GetNumFineScaleDisplacements()const;

    //! @brief gives the global DOF of a displacement component
    //! @param rComponent component
    //! @return global DOF
    int GetDofFineScaleDisplacement(int rComponent)const;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    void SetFineScaleDisplacements2D(const double rDisplacements[2]);

    //! @brief writes the displacements of a node to the prescribed pointer
    //! @param rDisplacements displacements
    void GetFineScaleDisplacements2D(double rDisplacements[2])const;

    //! @brief writes the displacements of a node to the prescribed pointer
    //! this is the sum of the fine scale and the coarse scale displacement
    //! @param rDisplacements displacements
    void GetDisplacements2D(double rDisplacements[2])const;

    //! @brief writes the displacements of a node to the prescribed pointer
    //! the difference is e.g. using XFEM, when the nodal degrees of freedom are not identical
    //! @param rDisplacements displacements
    double GetDisplacement(short rIndex)const;

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

    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    virtual Node::eNodeType GetNodeType()const;

    //! @brief returns, if the node is on the boundary of the cracked (true) or homogeneous domain (false)
    bool IsInCrackedDomain()const
    {
    	return mCrackedDomain;
    }

protected:
    //! @brief constructor for serialize
    NodeDisplacementsMultiscale2D(){};

    double mFineScaleDisplacements[2];
    int mDOF[2];
    //! @brief if set to true, the crack has an influence on the displacements, if false, it's just the homogeneous part
    bool mCrackedDomain;
    StructureMultiscale* mStructureMultiscale;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NodeDisplacementsMultiscale2D)
#endif // ENABLE_SERIALIZATION

#endif //NodeDisplacementsMultiscale2D_H
