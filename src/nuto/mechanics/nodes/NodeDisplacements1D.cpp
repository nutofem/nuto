// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeDisplacements1D.h"
#include "nuto/math/FullMatrix.h"

//! @brief constructor
NuTo::NodeDisplacements1D::NodeDisplacements1D() : NodeBase ()
{
    this->mDisplacements[0]=0;
    this->mDOF[0]=-1;
}

//! @brief constructor
NuTo::NodeDisplacements1D::NodeDisplacements1D (const double rDisplacements[1])  : NodeBase ()
{
    this->mDisplacements[0]= rDisplacements[0];
    this->mDOF[0]=-1;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeDisplacements1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeDisplacements1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeDisplacements1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeDisplacements1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeDisplacements1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeDisplacements1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeDisplacements1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeDisplacements1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mDisplacements)
       & BOOST_SERIALIZATION_NVP(mDOF);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeDisplacements1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeDisplacements1D)
BOOST_CLASS_TRACKING(NuTo::NodeDisplacements1D, track_always)
#endif // ENABLE_SERIALIZATION

//! @brief returns the number of displacements of the node
//! @return number of displacements
int NuTo::NodeDisplacements1D::GetNumDisplacements()const
{
    return 1;
}

//! @brief gives the global DOF of a displacement component
//! @param rComponent component
//! @return global DOF
int NuTo::NodeDisplacements1D::GetDofDisplacement(int rComponent)const
{
	if (rComponent!=0)
		throw MechanicsException("[NuTo::NodeDisplacements1D::GetDofDisplacement] Node has only a single displacement component.");
    return mDOF[0];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
void NuTo::NodeDisplacements1D::SetDisplacements1D(const double rDisplacements[1])
{
    mDisplacements[0] = rDisplacements[0];
}

//! @brief writes the displacements of a node to the prescribed pointer
//! @param rDisplacements displacements
void NuTo::NodeDisplacements1D::GetDisplacements1D(double rDisplacements[1])const
{
    rDisplacements[0] = mDisplacements[0];
}

//! @brief writes the displacements of a node to the prescribed pointer
//! the difference is e.g. using XFEM, when the nodal degrees of freedom are not identical
//! @param rDisplacements displacements
double NuTo::NodeDisplacements1D::GetDisplacement(short rIndex)const
{
    if (rIndex==0)
        return mDisplacements[0];
    else
        throw MechanicsException("[NuTo::NodeDisplacements1D::GetDisplacement] node has only one displacements.");
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeDisplacements1D::SetGlobalDofs(int& rDOF)
{
    mDOF[0]=rDOF++;
}

//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacements1D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
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
	this->mDisplacements[0] = value;
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacements1D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
	int dof = this->mDOF[0];
	double value = this->mDisplacements[0];
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
void NuTo::NodeDisplacements1D::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeDisplacements1D::SetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative (velocity) dofs.");
}

//! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacements1D::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeDisplacements1D::GetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative (velocity) dofs.");
}

//! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacements1D::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeDisplacements1D::SetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative (acceleration) dofs.");
}

//! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeDisplacements1D::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeDisplacements1D::GetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative (acceleration) dofs.");
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeDisplacements1D::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    mDOF[0]=rMappingInitialToNewOrdering[mDOF[0]];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeDisplacements1D::GetNodeTypeStr()const
{
	return std::string("NodeDisplacements1D");
}

//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
NuTo::Node::eNodeType NuTo::NodeDisplacements1D::GetNodeType()const
{
    return Node::NodeDisplacements1D;
}
