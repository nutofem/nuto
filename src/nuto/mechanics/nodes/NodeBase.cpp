// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeBase.h"


// For serialization we need to include the classes below
// the export_implement command (part of export_guid) has to be triggered
// only once, else linking errors will occur. This file is able to maintain
// that strategy, thatswhy the inludes below and the export routines at the end of the file
#include "nuto/mechanics/nodes/NodeDof.h"
#include "nuto/mechanics/nodes/NodeCoordinates.h"
#include "nuto/mechanics/nodes/NodeDisplacements.h"
#include "nuto/mechanics/nodes/NodeWaterVolumeFraction.h"
#include "nuto/mechanics/nodes/NodeNonlocalEqPlasticStrain.h"
#include "nuto/mechanics/nodes/NodeNonlocalEqStrain.h"
#include "nuto/mechanics/nodes/NodeRelativeHumidity.h"
#include "nuto/mechanics/nodes/NodeRotations.h"

//! @brief constructor
NuTo::NodeBase::NodeBase()
{
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::NodeBase::SetGlobalDofs(int& rDOF)
{
	throw MechanicsException("[NuTo::NodeBase::SetGlobalDofs] Node of type " + GetNodeTypeStr() + " dofs.");
}

//! @brief sets the global dofs numbers for each dof type
//! @param rDofNumbers ... map containing the dof type and the current number
void NuTo::NodeBase::SetGlobalDofsNumbers(std::map<Node::eDof, int>& rDofNumbers)
{
    throw MechanicsException("[NuTo::NodeBase::SetGlobalDofsNumbers] Node of type " + GetNodeTypeStr() + " dofs.");
}


//! @brief write dof values to the node (based on global dof number)
//! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeBase::SetGlobalDofValues(int rTimeDerivative, const FullVector<double,Eigen::Dynamic>& rActiveDofValues, const FullVector<double,Eigen::Dynamic>& rDependentDofValues)
{
	throw MechanicsException("[NuTo::NodeBase::GetGlobalDofValues] Node of type " + GetNodeTypeStr() + " dofs.");
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeBase::GetGlobalDofValues(int rTimeDerivative, FullVector<double,Eigen::Dynamic>& rActiveDofValues, FullVector<double,Eigen::Dynamic>& rDependentDofValues) const
{
	throw MechanicsException("[NuTo::NodeBase::GetGlobalDofValues] Node of type " + GetNodeTypeStr() + " dofs.");
}

//! @brief write dof values to the node (based on global dof number)
//! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
//! @param rDofType ... specific dof type
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeBase::SetGlobalDofValues(int rTimeDerivative, Node::eDof rDofType, const FullVector<double,Eigen::Dynamic>& rActiveDofValues, const FullVector<double,Eigen::Dynamic>& rDependentDofValues)
{
    throw MechanicsException("[NuTo::NodeBase::GetGlobalDofValues] Node of type " + GetNodeTypeStr() + " dofs.");
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
//! @param rDofType ... specific dof type
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::NodeBase::GetGlobalDofValues(int rTimeDerivative, Node::eDof rDofType, FullVector<double,Eigen::Dynamic>& rActiveDofValues, FullVector<double,Eigen::Dynamic>& rDependentDofValues) const
{
    throw MechanicsException("[NuTo::NodeBase::GetGlobalDofValues] Node of type " + GetNodeTypeStr() + " dofs.");
}

//! @brief extract all dof numbers from the node (based on global dof number)
//int* NuTo::NodeBase::GetGlobalDofs()
//{
//	throw MechanicsException("[NuTo::NodeBase::GetGlobalDofs] Node of type " + GetNodeTypeStr() + " dofs.");
//}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeBase::RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
	throw MechanicsException("[NuTo::NodeBase::RenumberGlobalDofs] Node of type " + GetNodeTypeStr() + " dofs.");
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rDofType ... specific dof type
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::NodeBase::RenumberGlobalDofs(Node::eDof rDofType, std::vector<int>& rMappingInitialToNewOrdering)
{
    throw MechanicsException("[NuTo::NodeBase::RenumberGlobalDofs] Node of type " + GetNodeTypeStr() + " dofs.");
}

//! @brief returns the number of time derivatives stored at the node
//! @return number of derivatives
int NuTo::NodeBase::GetNumTimeDerivatives()const
{
    throw MechanicsException("[NuTo::NodeBase::GetNumTimeDerivatives] Node of type " + GetNodeTypeStr() + " has no dofs");
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NodeBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeBase::serialize(Archive & ar, const unsigned int version)
{

#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::NodeBase)
BOOST_CLASS_TRACKING(NuTo::NodeBase, track_always)
#endif // ENABLE_SERIALIZATION

//*************************************************
//************     COORDINATES      ***************
//*************************************************
int NuTo::NodeBase::GetNumCoordinates() const
{
    return 0;
}
double NuTo::NodeBase::GetCoordinate(short rComponent) const
{throw MechanicsException("[NuTo::NodeBase::GetCoordinate] Node of type " + GetNodeTypeStr() + " has no coordinates.");}
const Eigen::Matrix<double, 1, 1>& NuTo::NodeBase::GetCoordinates1D() const
{throw MechanicsException("[NuTo::NodeBase::GetCoordinates1D] Node of type " + GetNodeTypeStr() + " has no coordinates.");}
const Eigen::Matrix<double, 2, 1>& NuTo::NodeBase::GetCoordinates2D() const
{throw MechanicsException("[NuTo::NodeBase::GetCoordinates2D] Node of type " + GetNodeTypeStr() + " has no coordinates.");}
const Eigen::Matrix<double, 3, 1>& NuTo::NodeBase::GetCoordinates3D() const
{throw MechanicsException("[NuTo::NodeBase::GetCoordinates3D] Node of type " + GetNodeTypeStr() + " has no coordinates.");}
const Eigen::Matrix<double, Eigen::Dynamic, 1> NuTo::NodeBase::GetCoordinates() const
{throw MechanicsException("[NuTo::NodeBase::GetCoordinates] Node of type " + GetNodeTypeStr() + " has no coordinates.");}

void NuTo::NodeBase::SetCoordinates1D(const Eigen::Matrix<double, 1, 1>& rCoordinates)
{throw MechanicsException("[NuTo::NodeBase::SetCoordinates1D] Node of type " + GetNodeTypeStr() + " has no coordinates.");}
void NuTo::NodeBase::SetCoordinates2D(const Eigen::Matrix<double, 2, 1>& rCoordinates)
{throw MechanicsException("[NuTo::NodeBase::SetCoordinates2D] Node of type " + GetNodeTypeStr() + " has no coordinates.");}
void NuTo::NodeBase::SetCoordinates3D(const Eigen::Matrix<double, 3, 1>& rCoordinates)
{throw MechanicsException("[NuTo::NodeBase::SetCoordinates3D] Node of type " + GetNodeTypeStr() + " has no coordinates.");}
void NuTo::NodeBase::SetCoordinates  (const Eigen::Matrix<double, Eigen::Dynamic, 1>& rCoordinates)
{throw MechanicsException("[NuTo::NodeBase::SetCoordinates] Node of type " + GetNodeTypeStr() + " has no coordinates.");}


//*************************************************
//************    DISPLACEMENTS     ***************
//*************************************************
int NuTo::NodeBase::GetNumDisplacements() const
{
    return 0;
}
int NuTo::NodeBase::GetDofDisplacement(int rComponent) const
{throw MechanicsException("[NuTo::NodeBase::GetDofDisplacement] Node of type " + GetNodeTypeStr() + " has no displacements.");}
double NuTo::NodeBase::GetDisplacement(short rIndex) const
{throw MechanicsException("[NuTo::NodeBase::GetDisplacement] Node of type " + GetNodeTypeStr() + " has no displacements.");}

const Eigen::Matrix<double, 1, 1>& NuTo::NodeBase::GetDisplacements1D() const
{
    return GetDisplacements1D(0);
}
const Eigen::Matrix<double, 2, 1>& NuTo::NodeBase::GetDisplacements2D() const
{
    return GetDisplacements2D(0);
}
const Eigen::Matrix<double, 3, 1>& NuTo::NodeBase::GetDisplacements3D() const
{
    return GetDisplacements3D(0);
}
const Eigen::Matrix<double, Eigen::Dynamic, 1> NuTo::NodeBase::GetDisplacements() const
{
    return GetDisplacements(0);
}

const Eigen::Matrix<double, 1, 1>& NuTo::NodeBase::GetDisplacements1D(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetDisplacements1D] Node of type " + GetNodeTypeStr() + " has no displacements.");}
const Eigen::Matrix<double, 2, 1>& NuTo::NodeBase::GetDisplacements2D(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetDisplacements2D] Node of type " + GetNodeTypeStr() + " has no displacements.");}
const Eigen::Matrix<double, 3, 1>& NuTo::NodeBase::GetDisplacements3D(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetDisplacements3D] Node of type " + GetNodeTypeStr() + " has no displacements.");}
const Eigen::Matrix<double, Eigen::Dynamic, 1> NuTo::NodeBase::GetDisplacements(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetDisplacements] Node of type " + GetNodeTypeStr() + " has no displacements.");}

void NuTo::NodeBase::SetDisplacements1D(const Eigen::Matrix<double, 1, 1>& rDisplacements)
{
    SetDisplacements1D(0, rDisplacements);
}
void NuTo::NodeBase::SetDisplacements2D(const Eigen::Matrix<double, 2, 1>& rDisplacements)
{
    SetDisplacements2D(0, rDisplacements);
}
void NuTo::NodeBase::SetDisplacements3D(const Eigen::Matrix<double, 3, 1>& rDisplacements)
{
    SetDisplacements3D(0, rDisplacements);
}
void NuTo::NodeBase::SetDisplacements(const Eigen::Matrix<double, Eigen::Dynamic, 1>& rDisplacements)
{
    SetDisplacements(0, rDisplacements);
}

void NuTo::NodeBase::SetDisplacements1D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rDisplacements)
{throw MechanicsException("[NuTo::NodeBase::SetDisplacements1D] Node of type " + GetNodeTypeStr() + " has no displacements.");}
void NuTo::NodeBase::SetDisplacements2D(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rDisplacements)
{throw MechanicsException("[NuTo::NodeBase::SetDisplacements2D] Node of type " + GetNodeTypeStr() + " has no displacements.");}
void NuTo::NodeBase::SetDisplacements3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rDisplacements)
{throw MechanicsException("[NuTo::NodeBase::SetDisplacements3D] Node of type " + GetNodeTypeStr() + " has no displacements.");}
void NuTo::NodeBase::SetDisplacements  (int rTimeDerivative, const Eigen::Matrix<double, Eigen::Dynamic, 1>& rDisplacements)
{throw MechanicsException("[NuTo::NodeBase::SetDisplacements] Node of type " + GetNodeTypeStr() + " has no displacements.");}
//*************************************************
//************       ROTATIONS      ***************
//*************************************************

int NuTo::NodeBase::GetNumRotations()const
{
    return 0;
}
int NuTo::NodeBase::GetDofRotation(int rComponent)const
{throw MechanicsException("[NuTo::NodeBase::GetDofRotation] Node of type " + GetNodeTypeStr() + " has no rotations.");}
double NuTo::NodeBase::GetRotation(short rIndex)const
{throw MechanicsException("[NuTo::NodeBase::GetRotation] Node of type " + GetNodeTypeStr() + " has no rotations.");}

const Eigen::Matrix<double, 1, 1>& NuTo::NodeBase::GetRotations2D() const
{
    return GetRotations2D(0);
}
const Eigen::Matrix<double, 3, 1>& NuTo::NodeBase::GetRotations3D() const
{
    return GetRotations3D(0);
}
const Eigen::Matrix<double, Eigen::Dynamic, 1> NuTo::NodeBase::GetRotations() const
{
    return GetRotations(0);
}

const Eigen::Matrix<double, 1, 1>& NuTo::NodeBase::GetRotations2D(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetRotations2D] Node of type " + GetNodeTypeStr() + " has no rotations.");}
const Eigen::Matrix<double, 3, 1>& NuTo::NodeBase::GetRotations3D(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetRotations3D] Node of type " + GetNodeTypeStr() + " has no rotations.");}
const Eigen::Matrix<double, Eigen::Dynamic, 1> NuTo::NodeBase::GetRotations(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetRotations] Node of type " + GetNodeTypeStr() + " has no rotations.");}

void NuTo::NodeBase::SetRotations2D(const Eigen::Matrix<double, 1, 1>& rRotations)
{
    SetRotations2D(0, rRotations);
}
void NuTo::NodeBase::SetRotations3D(const Eigen::Matrix<double, 3, 1>& rRotations)
{
    SetRotations3D(0, rRotations);
}
void NuTo::NodeBase::SetRotations(const Eigen::Matrix<double, Eigen::Dynamic, 1>& rRotations)
{
    SetRotations(0, rRotations);
}

void NuTo::NodeBase::SetRotations2D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rRotations)
{throw MechanicsException("[NuTo::NodeBase::SetRotations2D] Node of type " + GetNodeTypeStr() + " has no rotations.");}
void NuTo::NodeBase::SetRotations3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rRotations)
{throw MechanicsException("[NuTo::NodeBase::SetRotations3D] Node of type " + GetNodeTypeStr() + " has no rotations.");}
void NuTo::NodeBase::SetRotations(int rTimeDerivative, const Eigen::Matrix<double, Eigen::Dynamic, 1>& rRotations)
{throw MechanicsException("[NuTo::NodeBase::SetRotations] Node of type " + GetNodeTypeStr() + " has no rotations.");}



//*************************************************
//************      TEMPERATURE     ***************
//*************************************************

int NuTo::NodeBase::GetNumTemperature() const
{
    return 0;
}
int NuTo::NodeBase::GetDofTemperature() const
{throw MechanicsException("[NuTo::NodeBase::GetDofTemperature] Node of type " + GetNodeTypeStr() + " has no temperature.");}

double NuTo::NodeBase::GetTemperature() const
{
    return GetTemperature(0);
}
double NuTo::NodeBase::GetTemperature(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetTemperature] Node of type " + GetNodeTypeStr() + " has no temperature.");}

void NuTo::NodeBase::SetTemperature(double rTemperature)
{
    SetTemperature(0, rTemperature);
}
void NuTo::NodeBase::SetTemperature(int rTimeDerivative, double rTemperature)
{throw MechanicsException("[NuTo::NodeBase::SetTemperature] Node of type " + GetNodeTypeStr() + " has no temperature.");}

//*************************************************
//********  NONLOCAL EQ PLASTIC STRAIN  ***********
//*************************************************

int NuTo::NodeBase::GetNumNonlocalEqPlasticStrain() const
{
    return 0;
}
int NuTo::NodeBase::GetDofNonlocalEqPlasticStrain(int rComponent) const
{throw MechanicsException("[NuTo::NodeBase::GetDofNonlocalEqPlasticStrain] Node of type " + GetNodeTypeStr() + " has no nonlocal eq plastic strain.");}

const Eigen::Matrix<double, 2, 1>& NuTo::NodeBase::GetNonlocalEqPlasticStrain() const
{
    return GetNonlocalEqPlasticStrain(0);
}
const Eigen::Matrix<double, 2, 1>& NuTo::NodeBase::GetNonlocalEqPlasticStrain(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetNonlocalEqPlasticStrain] Node of type " + GetNodeTypeStr() + " has no nonlocal eq plastic strain.");}

void NuTo::NodeBase::SetNonlocalEqPlasticStrain(const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain)
{
    SetNonlocalEqPlasticStrain(0, rNonlocalEqPlasticStrain);
}
void NuTo::NodeBase::SetNonlocalEqPlasticStrain(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain)
{throw MechanicsException("[NuTo::NodeBase::SetNonlocalEqPlasticStrain] Node of type " + GetNodeTypeStr() + " has no nonlocal eq plastic strain.");}

//*************************************************
//********    NONLOCAL TOTAL STRAIN     ***********
//*************************************************

int NuTo::NodeBase::GetNumNonlocalTotalStrain()const
{
    return 0;
}
int NuTo::NodeBase::GetDofNonlocalTotalStrain(int rComponent)const
{throw MechanicsException("[NuTo::NodeBase::GetDofNonlocalTotalStrain] Node of type " + GetNodeTypeStr() + " has no nonlocal total strain.");}
double NuTo::NodeBase::GetNonlocalTotalStrain(short rIndex)const
{throw MechanicsException("[NuTo::NodeBase::GetNonlocalTotalStrain] Node of type " + GetNodeTypeStr() + " has no nonlocal total strain.");}

const Eigen::Matrix<double, 1, 1>& NuTo::NodeBase::GetNonlocalTotalStrain1D() const
{
    return GetNonlocalTotalStrain1D(0);
}
const Eigen::Matrix<double, 3, 1>& NuTo::NodeBase::GetNonlocalTotalStrain2D() const
{
    return GetNonlocalTotalStrain2D(0);
}
const Eigen::Matrix<double, 6, 1>& NuTo::NodeBase::GetNonlocalTotalStrain3D() const
{
    return GetNonlocalTotalStrain3D(0);
}
const Eigen::Matrix<double, Eigen::Dynamic, 1> NuTo::NodeBase::GetNonlocalTotalStrains() const
{
    return GetNonlocalTotalStrains(0);
}

const Eigen::Matrix<double, 1, 1>& NuTo::NodeBase::GetNonlocalTotalStrain1D(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetNonlocalTotalStrain1D] Node of type " + GetNodeTypeStr() + " has no nonlocal total strain.");}
const Eigen::Matrix<double, 3, 1>& NuTo::NodeBase::GetNonlocalTotalStrain2D(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetNonlocalTotalStrain2D] Node of type " + GetNodeTypeStr() + " has no nonlocal total strain.");}
const Eigen::Matrix<double, 6, 1>& NuTo::NodeBase::GetNonlocalTotalStrain3D(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetNonlocalTotalStrain3D] Node of type " + GetNodeTypeStr() + " has no nonlocal total strain.");}
const Eigen::Matrix<double, Eigen::Dynamic, 1> NuTo::NodeBase::GetNonlocalTotalStrains(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetNonlocalTotalStrains] Node of type " + GetNodeTypeStr() + " has no nonlocal total strain.");}

void NuTo::NodeBase::SetNonlocalTotalStrain1D(const Eigen::Matrix<double, 1, 1>& rNonlocalTotalStrain)
{
    SetNonlocalTotalStrain1D(0, rNonlocalTotalStrain);
}
void NuTo::NodeBase::SetNonlocalTotalStrain2D(const Eigen::Matrix<double, 3, 1>& rNonlocalTotalStrain)
{
    SetNonlocalTotalStrain2D(0, rNonlocalTotalStrain);
}
void NuTo::NodeBase::SetNonlocalTotalStrain3D(const Eigen::Matrix<double, 6, 1>& rNonlocalTotalStrain)
{
    SetNonlocalTotalStrain3D(0, rNonlocalTotalStrain);
}
void NuTo::NodeBase::SetNonlocalTotalStrain  (const Eigen::Matrix<double, Eigen::Dynamic, 1>& rNonlocalTotalStrain)
{
    SetNonlocalTotalStrain(0, rNonlocalTotalStrain);
}


void NuTo::NodeBase::SetNonlocalTotalStrain1D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rNonlocalTotalStrain)
{throw MechanicsException("[NuTo::NodeBase::SetNonlocalTotalStrain1D] Node of type " + GetNodeTypeStr() + " has no nonlocal total strain.");}
void NuTo::NodeBase::SetNonlocalTotalStrain2D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rNonlocalTotalStrain)
{throw MechanicsException("[NuTo::NodeBase::SetNonlocalTotalStrain2D] Node of type " + GetNodeTypeStr() + " has no nonlocal total strain.");}
void NuTo::NodeBase::SetNonlocalTotalStrain3D(int rTimeDerivative, const Eigen::Matrix<double, 6, 1>& rNonlocalTotalStrain)
{throw MechanicsException("[NuTo::NodeBase::SetNonlocalTotalStrain3D] Node of type " + GetNodeTypeStr() + " has no nonlocal total strain.");}
void NuTo::NodeBase::SetNonlocalTotalStrain  (int rTimeDerivative, const Eigen::Matrix<double, Eigen::Dynamic, 1>& rNonlocalTotalStrain)
{throw MechanicsException("[NuTo::NodeBase::SetNonlocalTotalStrain] Node of type " + GetNodeTypeStr() + " has no nonlocal total strain.");}

//*************************************************
//*******   NONLOCAL EQUIVALENT STRAIN   **********
//*************************************************

int NuTo::NodeBase::GetNumNonlocalEqStrain() const
{
    return 0;
}
int NuTo::NodeBase::GetDofNonlocalEqStrain() const
{throw MechanicsException("[NuTo::NodeBase::GetDofNonlocalEqStrain] Node of type " + GetNodeTypeStr() + " has no nonlocal eq strain.");}

double NuTo::NodeBase::GetNonlocalEqStrain() const
{
    return GetNonlocalEqStrain(0);
}
double NuTo::NodeBase::GetNonlocalEqStrain(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetNonlocalEqStrain] Node of type " + GetNodeTypeStr() + " has no nonlocal eq strain.");}

void NuTo::NodeBase::SetNonlocalEqStrain(double rNonlocalEqStrain)
{
    SetNonlocalEqStrain(0, rNonlocalEqStrain);
}
void NuTo::NodeBase::SetNonlocalEqStrain(int rTimeDerivative, double rNonlocalEqStrain)
{throw MechanicsException("[NuTo::NodeBase::SetNonlocalEqStrain] Node of type " + GetNodeTypeStr() + " has no nonlocal eq strain.");}

//**************************************************
//*******      WATER VOLUME FRACTION      **********
//**************************************************

int NuTo::NodeBase::GetNumWaterVolumeFraction() const
{
    return 0;
}
int NuTo::NodeBase::GetDofWaterVolumeFraction() const
{throw MechanicsException("[NuTo::NodeBase::GetDofWaterVolumeFraction] Node of type " + GetNodeTypeStr() + " has no water volume fraction.");}

double NuTo::NodeBase::GetWaterVolumeFraction() const
{
    return GetWaterVolumeFraction(0);
}
double NuTo::NodeBase::GetWaterVolumeFraction(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetWaterVolumeFraction] Node of type " + GetNodeTypeStr() + " has no water volume fraction.");}

void NuTo::NodeBase::SetWaterVolumeFraction(double rWaterVolumeFraction)
{
    SetWaterVolumeFraction(0, rWaterVolumeFraction);
}
void NuTo::NodeBase::SetWaterVolumeFraction(int rTimeDerivative, double rWaterVolumeFraction)
{throw MechanicsException("[NuTo::NodeBase::SetWaterVolumeFraction] Node of type " + GetNodeTypeStr() + " has no water volume fraction.");}

//*************************************************
//*******       RELATIVE HUMIDITY        **********
//*************************************************

int NuTo::NodeBase::GetNumRelativeHumidity() const
{
    return 0;
}
int NuTo::NodeBase::GetDofRelativeHumidity() const
{throw MechanicsException("[NuTo::NodeBase::GetDofRelativeHumidity] Node of type " + GetNodeTypeStr() + " has no relative humidity.");}

double NuTo::NodeBase::GetRelativeHumidity() const
{
    return GetRelativeHumidity(0);
}
double NuTo::NodeBase::GetRelativeHumidity(int rTimeDerivative) const
{throw MechanicsException("[NuTo::NodeBase::GetRelativeHumidity] Node of type " + GetNodeTypeStr() + " has no relative humidity.");}

void NuTo::NodeBase::SetRelativeHumidity(double rRelativeHumidity)
{
    SetRelativeHumidity(0, rRelativeHumidity);
}
void NuTo::NodeBase::SetRelativeHumidity(int rTimeDerivative, double rRelativeHumidity)
{throw MechanicsException("[NuTo::NodeBase::SetRelativeHumidity] Node of type " + GetNodeTypeStr() + " has no relative humidity.");}


#ifdef ENABLE_VISUALIZE
void NuTo::NodeBase::Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList) const
{
    Eigen::Matrix<double, 3, 1> coordinates = Eigen::Matrix<double, 3, 1>::Zero();
	switch (this->GetNumCoordinates())
	{
	case 1:
	    coordinates.block<1,1>(0,0) = this->GetCoordinates1D();
		break;
	case 2:
        coordinates.block<2,1>(0,0) = this->GetCoordinates2D();
		break;
	case 3:
        coordinates.block<3,1>(0,0) = this->GetCoordinates3D();
		break;
	default:
		throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither coordinates in 1D, 2D or 3D.");
	}
	unsigned int PointId = rVisualize.AddPoint(coordinates.data());
//	std::cout << "add point " << PointId << std::endl;

    // store data
    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
			case NuTo::VisualizeBase::DISPLACEMENTS:
			{
			    if (GetNumDisplacements() == 0)
			        break;
			    Eigen::Matrix<double, 3, 1> displacements = Eigen::Matrix<double, 3, 1>::Zero();
				switch (this->GetNumDisplacements())
				{
				case 1:
				    displacements.block<1,1>(0,0) = this->GetDisplacements1D();
					break;
				case 2:
                    displacements.block<2,1>(0,0) = this->GetDisplacements2D();
					break;
				case 3:
                    displacements.block<3,1>(0,0) = this->GetDisplacements3D();
					break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither displacements in 1D, 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), displacements.data());
			}
				break;
			case NuTo::VisualizeBase::ROTATION:
			{
                Eigen::Matrix<double, 3, 1> rotations = Eigen::Matrix<double, 3, 1>::Zero();
				switch (this->GetNumRotations())
				{
				case 1:
				    rotations.block<1,1>(0,0) = this->GetRotations2D();
					break;
				case 3:
                    rotations.block<3,1>(0,0) = this->GetRotations3D();
					break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither rotations in 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), rotations.data());
			}
				break;
			case NuTo::VisualizeBase::VELOCITY:
			{
                Eigen::Matrix<double, 3, 1> velocities = Eigen::Matrix<double, 3, 1>::Zero();
                switch (this->GetNumDisplacements())
                {
                case 1:
                    velocities.block<1,1>(0,0) = this->GetDisplacements1D(1);
                    break;
                case 2:
                    velocities.block<2,1>(0,0) = this->GetDisplacements2D(1);
                    break;
                case 3:
                    velocities.block<3,1>(0,0) = this->GetDisplacements3D(1);
					break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither velocities in 1D, 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), velocities.data());
			}
				break;
			case NuTo::VisualizeBase::ANGULAR_VELOCITY:
			{
                Eigen::Matrix<double, 3, 1> angularVelocities = Eigen::Matrix<double, 3, 1>::Zero();
                switch (this->GetNumRotations())
                {
                case 1:
                    angularVelocities.block<1,1>(0,0) = this->GetRotations2D(1);
                    break;
                case 3:
                    angularVelocities.block<3,1>(0,0) = this->GetRotations3D(1);
                    break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither angular velocities in 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), angularVelocities.data());
			}
				break;
			case NuTo::VisualizeBase::ACCELERATION:
			{
                Eigen::Matrix<double, 3, 1> accelerations = Eigen::Matrix<double, 3, 1>::Zero();
                switch (this->GetNumDisplacements())
                {
                case 1:
                    accelerations.block<1,1>(0,0) = this->GetDisplacements1D(2);
                    break;
                case 2:
                    accelerations.block<2,1>(0,0) = this->GetDisplacements2D(2);
                    break;
                case 3:
                    accelerations.block<3,1>(0,0) = this->GetDisplacements3D(2);
                    break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither accelerations in 1D, 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), accelerations.data());
			}
				break;
			case NuTo::VisualizeBase::ANGULAR_ACCELERATION:
			{
                Eigen::Matrix<double, 3, 1> angularAccelerations = Eigen::Matrix<double, 3, 1>::Zero();
                switch (this->GetNumRotations())
                {
                case 1:
                    angularAccelerations.block<1,1>(0,0) = this->GetRotations2D(2);
                    break;
                case 3:
                    angularAccelerations.block<3,1>(0,0) = this->GetRotations3D(2);
                    break;
				default:
					throw MechanicsException("[NuTo::NodeBase::Visualize] node has neither angular accelerations in 2D or 3D.");
				}
					rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), angularAccelerations.data());
			}
				break;
			default:
				break;
        }
    }
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
namespace boost{
template<NODE_DOF_TEMPLATE_PARAMETERS>
    struct is_virtual_base_of<NuTo::NodeCoordinates<TNumCoordinates>, NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>>: public mpl::true_ {};

template<NODE_DOF_TEMPLATE_PARAMETERS>
    struct is_virtual_base_of<NuTo::NodeDisplacements<TNumDisplacements, TNumTimeDerivatives>, NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>>: public mpl::true_ {};

template<NODE_DOF_TEMPLATE_PARAMETERS>
    struct is_virtual_base_of<NuTo::NodeRotations<TNumRotations, TNumTimeDerivatives>, NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>>: public mpl::true_ {};

template<NODE_DOF_TEMPLATE_PARAMETERS>
    struct is_virtual_base_of<NuTo::NodeNonlocalEqPlasticStrain<TNumNonlocalEqPlasticStrain, TNumTimeDerivatives>, NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>>: public mpl::true_ {};

template<NODE_DOF_TEMPLATE_PARAMETERS>
    struct is_virtual_base_of<NuTo::NodeNonlocalEqStrain<TNumNonlocalEqStrain, TNumTimeDerivatives>, NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>>: public mpl::true_ {};

template<NODE_DOF_TEMPLATE_PARAMETERS>
    struct is_virtual_base_of<NuTo::NodeRelativeHumidity<TNumRelativeHumidity, TNumTimeDerivatives>, NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>>: public mpl::true_ {};

template<NODE_DOF_TEMPLATE_PARAMETERS>
    struct is_virtual_base_of<NuTo::NodeWaterVolumeFraction<TNumWaterVolumeFraction,TNumTimeDerivatives>, NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>>: public mpl::true_ {};
}

BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 1, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 1, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 1, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 1, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 2, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 2, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 3, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 3, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 1, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 1, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 0, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 0, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 1, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 1, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 3, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 3, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 1, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 1, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 3, 0, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 3, 0, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 1, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 1, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 1, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 1, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 1, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 1, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 0, 0, 1, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 0, 0, 1, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 0, 0, 1, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 0, 0, 1, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 0, 0, 1, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 0, 0, 1, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 0, 0, 1, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 0, 0, 1, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 0, 0, 1, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 0, 0, 1, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 0, 0, 1, 0, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 0, 0, 1, 0, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 1, 0, 0, 2, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 1, 0, 0, 2, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 0, 0, 2, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 0, 0, 2, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 0, 0, 2, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 0, 0, 2, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 1, 0, 0, 2, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 1, 0, 0, 2, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 0, 0, 2, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 0, 0, 2, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 0, 0, 2, 0, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 0, 0, 2, 0, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 1, 0, 0, 0, 1, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 1, 0, 0, 0, 1, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 0, 0, 0, 3, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 0, 0, 0, 3, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 0, 0, 0, 6, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 0, 0, 0, 6, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 1, 0, 0, 0, 1, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 1, 0, 0, 0, 1, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 0, 0, 0, 3, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 0, 0, 0, 3, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 0, 0, 0, 6, 0, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 0, 0, 0, 6, 0, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 1, 0, 0, 0, 0, 1, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 1, 0, 0, 0, 0, 1, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 0, 0, 0, 0, 1, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 0, 0, 0, 0, 1, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 0, 0, 0, 0, 1, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 0, 0, 0, 0, 1, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 1, 0, 0, 0, 0, 1, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 1, 0, 0, 0, 0, 1, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 0, 0, 0, 0, 1, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 0, 0, 0, 0, 1, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 0, 0, 0, 0, 1, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 0, 0, 0, 0, 1, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 1, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 1, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 1, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 1, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 1, 0, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 1, 0, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 1, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 1, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 2, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 3, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 1, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 1, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 2, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 2, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 3, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 3, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 1, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 1, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 2, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 3, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 0, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 0, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 0, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 0, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 0, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 0, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 0, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 0, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 0, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 0, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 0, 0, 0, 0, 0, 0, 1, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 0, 0, 0, 0, 0, 0, 1, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 0, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 0, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 0, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 0, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 0, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 0, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 0, 0, 0, 0, 0, 0, 0, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 0, 0, 0, 0, 0, 0, 0, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 0, 0, 0, 0, 0, 0, 0, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 0, 0, 0, 0, 0, 0, 0, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 0, 0, 0, 0, 0, 0, 0, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 0, 0, 0, 0, 0, 0, 0, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 0, 0, 0, 0, 0, 0, 0, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 0, 0, 0, 0, 0, 0, 0, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 0, 0, 0, 0, 0, 0, 0, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 0, 0, 0, 0, 0, 0, 0, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 0, 0, 0, 0, 0, 0, 0, 1>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 0, 0, 0, 0, 0, 0, 0, 1>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 1, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 1, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 1, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 0, 0, 0, 0, 0, 0, 0, 1, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 1, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 0, 0, 0, 0, 0, 0, 0, 1, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 0, 0, 0, 0, 0, 0, 1, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 1, 0, 0, 0, 0, 0, 0, 1, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 0, 0, 0, 0, 0, 0, 1, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 1, 0, 0, 0, 0, 0, 0, 1, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 0, 0, 0, 0, 0, 0, 1, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 1, 0, 0, 0, 0, 0, 0, 1, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 0, 0, 0, 0, 0, 0, 1, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1, 2, 0, 0, 0, 0, 0, 0, 1, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 0, 0, 0, 0, 0, 0, 1, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2, 2, 0, 0, 0, 0, 0, 0, 1, 0>))))
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 0, 0, 0, 0, 0, 0, 1, 0>)), BOOST_PP_STRINGIZE(BOOST_IDENTITY_TYPE((NuTo::NodeDof<3, 2, 0, 0, 0, 0, 0, 0, 1, 0>))))
#endif
