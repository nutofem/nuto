#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeTemperature.h"

//! @brief constructor
NuTo::NodeTemperature::NodeTemperature() : NodeBase ()
{
	this->mTemperature=0;
	this->mDOF=-1;
}

//! @brief constructor
NuTo::NodeTemperature::NodeTemperature (double rTemperature)  : NodeBase ()
{
	this->mTemperature=rTemperature;
	this->mDOF=-1;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::NodeTemperature::serialize(Archive & ar, const unsigned int version)
{
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
	   & BOOST_SERIALIZATION_NVP(mTemperature);
}
#endif  // ENABLE_SERIALIZATION

//! @brief returns the number of temperatures of the node
//! @return number of temperatures
int NuTo::NodeTemperature::GetNumTemperatures()const
{
	return 1;
}

//! @brief gives the global DOF of a temperature component
//! @param rComponent component
//! @return global DOF
int NuTo::NodeTemperature::GetDofTemperature(int rComponent)const
{
	assert(rComponent==0);
	return mDOF;
}

//! @brief set the rotations
//! @param rRotations  given rotations
void NuTo::NodeTemperature::SetTemperature(double rTemperature)
{
	mTemperature=rTemperature;
}

//! @brief writes the temperature of a node to the prescribed pointer
//! @param rTemperatur temperature
void NuTo::NodeTemperature::GetTemperature(double& rTemperature)const
{
	rTemperature=mTemperature;
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeTemperature::SetGlobalDofs(int& rDOF)
{
	mDOF=rDOF++;
}


//! @brief write dof values to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeTemperature::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	assert(rActiveDofValues.GetNumColumns() == 1);
	assert(rDependentDofValues.GetNumColumns() == 1);
	int dof = this->mDOF;
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
	this->mTemperature = value;
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeTemperature::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	assert(rActiveDofValues.GetNumColumns() == 1);
	assert(rDependentDofValues.GetNumColumns() == 1);
	int dof = this->mDOF;
	double value = this->mTemperature;
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
void NuTo::NodeTemperature::SetGlobalDofFirstTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeTemperature::SetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative dofs.");
}

//! @brief extract first time derivative of the dof values (e.g. velocities) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeTemperature::GetGlobalDofFirstTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeTemperature::GetGlobalDofFirstTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no first time derivative dofs.");
}

//! @brief write second time derivative of the dof values (e.g. accelerations) to the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeTemperature::SetGlobalDofSecondTimeDerivativeValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeTemperature::SetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative dofs.");
}

//! @brief extract second time derivative of the dof values (e.g. accelerations) from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeTemperature::GetGlobalDofSecondTimeDerivativeValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeTemperature::GetGlobalDofSecondTimeDerivativeValues] Node of type " + GetNodeTypeStr() + " has no second time derivative dofs.");
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeTemperature::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
	mDOF=rMappingInitialToNewOrdering[mDOF];
}

//! @brief returns the type of the node
//! @return type
std::string NuTo::NodeTemperature::GetNodeTypeStr()const
{
	return std::string("NodeTemperature");
}
