// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeRotations2D.h"

//! @brief constructor
NuTo::NodeRotations2D::NodeRotations2D()
{
	this->mRotations[0] = 0.;
	this->mDOF[0]=-1;
}

//! @brief constructor
NuTo::NodeRotations2D::NodeRotations2D(const double rRotations[1])  : NodeBase ()
{
	this->mRotations[0] = rRotations[0];
	this->mDOF[0]=-1;
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template void NuTo::NodeRotations2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
    template void NuTo::NodeRotations2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
    template void NuTo::NodeRotations2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
    template void NuTo::NodeRotations2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
    template void NuTo::NodeRotations2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
    template void NuTo::NodeRotations2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
    template<class Archive>
    void NuTo::NodeRotations2D::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
           & BOOST_SERIALIZATION_NVP(mRotations);
    }
    BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeRotations2D)
    BOOST_CLASS_TRACKING(NuTo::NodeRotations2D, track_always)
#endif // ENABLE_SERIALIZATION


//! @brief returns the number of Rotations of the node
//! @return number of Rotations
int NuTo::NodeRotations2D::GetNumRotations()const
{
	return 1;
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
void NuTo::NodeRotations2D::SetRotations2D(const double rRotations[1])
{
	mRotations[0] = rRotations[0];
}

//! @brief writes the Rotations of a node to the prescribed pointer
//! @param rRotations Rotations
void NuTo::NodeRotations2D::GetRotations2D(double rRotations[1])const
{
	rRotations[0] = mRotations[0];
}

//! @brief returns the Rotations of the node
//! @return Rotation
double NuTo::NodeRotations2D::GetRotation(short rComponent)const
{
	if (rComponent!=0)
		throw MechanicsException("[NuTo::NodeRotations2D::GetRotation] Node has only one rotation component.");
	return mRotations[0];
}

//! @brief returns the dof of the rotation of the node
//! @return Rotation
int NuTo::NodeRotations2D::GetDofRotation(int rComponent)const
{
	if (rComponent!=0)
		throw MechanicsException("[NuTo::NodeRotations2D::GetDofRotation] Node has only one rotation component.");
	return mDOF[0];
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeRotations2D::SetGlobalDofs(int& rDOF)
{
    mDOF[0]=rDOF++;
}

//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRotations2D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
	int dof = this->mDOF[0];
	double value;
	if (dof >= rActiveDofValues.GetNumRows())
	{
		dof -= rActiveDofValues.GetNumRows();
		assert(dof < rDependentDofValues.GetNumRows());
		value = rDependentDofValues(dof,0);
	}
	else
	{
		value = rActiveDofValues(dof,0);
	}
	this->mRotations[0] = value;
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRotations2D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
	int dof = this->mDOF[0];
	double value = this->mRotations[0];
	if (dof >= rActiveDofValues.GetNumRows())
	{
		dof -= rActiveDofValues.GetNumRows();
		assert(dof < rDependentDofValues.GetNumRows());
		rDependentDofValues(dof,0) = value;
	}
	else
	{
		rActiveDofValues(dof,0) = value;
    }
}

//! @brief write first time derivative of the dof values (e.g. velocities) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRotations2D::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeRotations2D::SetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative (velocity) dofs.");
}

//! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRotations2D::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeRotations2D::GetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative (velocity) dofs.");
}

//! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRotations2D::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeRotations2D::SetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative (acceleration) dofs.");
}

//! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeRotations2D::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeRotations2D::GetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative (acceleration) dofs.");
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeRotations2D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    mDOF[0]=rMappingInitialToNewOrdering[mDOF[0]];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeRotations2D::GetNodeTypeStr()const
{
	return std::string("NodeRotations2D");
}

//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
NuTo::Node::eNodeType NuTo::NodeRotations2D::GetNodeType()const
{
    return Node::NodeRotations2D;
}

