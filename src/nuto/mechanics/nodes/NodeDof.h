#pragma once

#include "nuto/mechanics/nodes/NodeDof_Def.h"
#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/MechanicsException.h"

//! @brief ... constructor
template<NODE_DOF_TEMPLATE_PARAMETERS>
NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
NodeDof() //: NuTo::NodeBase()
{

}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
template<NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetGlobalDofs(int& rDOF)
{
    for (int count = 0; count < TNumDisplacements; count++)
    {
        this->mDofDisplacements[count] = rDOF++;
    }
    for (int count = 0; count < TNumRotations; count++)
    {
        this->mDofRotations[count] = rDOF++;
    }
    for (int count=0; count<TNumTemperature; count++)
    {
        this->mDofTemperature[count] = rDOF++;
    }
    for (int count = 0; count < TNumNonlocalEqPlasticStrain; count++)
    {
        this->mDofNonlocalEqPlasticStrain[count] = rDOF++;
    }
//    for (int count=0; count<TNumNonlocalTotalStrain; count++)
//    {
//    	mDofNonlocalTotalStrain[count]=rDOF++;
//    }
    for (int count = 0; count < TNumNonlocalEqStrain; count++)
    {
        this->mDofNonlocalEqStrain[count] = rDOF++;
    }
    for (int count = 0; count < TNumWaterVolumeFraction; count++)
    {
        this->mDofWaterVolumeFraction[count] = rDOF++;
    }
    for (int count = 0; count < TNumRelativeHumidity; count++)
    {
        this->mDofRelativeHumidity[count] = rDOF++;
    }

}

//! @brief sets the global dofs numbers for each dof type
//! @param rDofNumbers ... map containing the dof type and the current number

template<NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetGlobalDofsNumbers(std::map<Node::eDof, int>& rDofNumbers)
{
    for (int i = 0; i < TNumDisplacements; ++i)
        this->mDofDisplacements[i] = rDofNumbers[Node::DISPLACEMENTS]++;

    for (int i = 0; i < TNumRotations; ++i)
        this->mDofRotations[i] = rDofNumbers[Node::ROTATIONS]++;

    for (int i = 0; i < TNumTemperature; ++i)
        this->mDofTemperature[i] = rDofNumbers[Node::TEMPERATURE]++;

    for (int i = 0; i < TNumNonlocalEqPlasticStrain; ++i)
        this->mDofNonlocalEqPlasticStrain[i] = rDofNumbers[Node::NONLOCALEQPLASTICSTRAIN]++;

//  for (int i = 0; i < TNumNonlocalTotalStrain; ++i)
//      this->mDofNonlocalTotalStrain[i]        = rDofNumbers[Node::NONLOCALEQTOTALSTRAIN]++;

    for (int i = 0; i < TNumNonlocalEqStrain; ++i)
        this->mDofNonlocalEqStrain[i] = rDofNumbers[Node::NONLOCALEQSTRAIN]++;

    for (int i = 0; i < TNumWaterVolumeFraction; ++i)
        this->mDofWaterVolumeFraction[i] = rDofNumbers[Node::WATERVOLUMEFRACTION]++;

    for (int i = 0; i < TNumRelativeHumidity; ++i)
        this->mDofRelativeHumidity[i] = rDofNumbers[Node::RELATIVEHUMIDITY]++;

}

//! @brief write dof values to the node (based on global dof number)
//! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
template<NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetGlobalDofValues(int rTimeDerivative, const NuTo::FullVector<double, Eigen::Dynamic>& rActiveDofValues, const NuTo::FullVector<double, Eigen::Dynamic>& rDependentDofValues)
{
    assert(TNumTimeDerivatives >= rTimeDerivative);

    for (int i = 0; i < TNumDisplacements; ++i)
        this->mDisplacements[rTimeDerivative][i] = GetDofValueFromVector(this->mDofDisplacements[i], rActiveDofValues, rDependentDofValues);

    for (int i = 0; i < TNumRotations; ++i)
        this->mRotations[rTimeDerivative][i] = GetDofValueFromVector(this->mDofRotations[i], rActiveDofValues, rDependentDofValues);


    for (int i = 0; i < TNumTemperature; ++i)
        this->mTemperature[rTimeDerivative][i] = GetDofValueFromVector(this->mDofTemperature[i], rActiveDofValues, rDependentDofValues);

    for (int i = 0; i < TNumNonlocalEqPlasticStrain; ++i)
    {
        double value = GetDofValueFromVector(this->mDofNonlocalEqPlasticStrain[i], rActiveDofValues, rDependentDofValues);
        if (value < 0)
            std::cout << "value is less then zero" << value << std::endl; // no further consequences... lucky you.
    }

//    for (int i=0; i<TNumNonlocalTotalStrain; ++i)
//        this->mNonlocalTotalStrain[rTimeDerivative][i] = GetDofValueFromVector(this->mDofNonlocalTotalStrain[i], rActiveDofValues, rDependentDofValues);

    for (int i = 0; i < TNumNonlocalEqStrain; ++i)
        this->mNonlocalEqStrain[rTimeDerivative][i] = GetDofValueFromVector(this->mDofNonlocalEqStrain[i], rActiveDofValues, rDependentDofValues);

    for (int i = 0; i < TNumWaterVolumeFraction; ++i)
        this->mWaterVolumeFraction[rTimeDerivative][i] = GetDofValueFromVector(this->mDofWaterVolumeFraction[i], rActiveDofValues, rDependentDofValues);

    for (int i = 0; i < TNumRelativeHumidity; ++i)
        this->mRelativeHumidity[rTimeDerivative][i] = GetDofValueFromVector(this->mDofRelativeHumidity[i], rActiveDofValues, rDependentDofValues);

}

//! @brief extract dof values from the node (based on global dof number)
//! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
template<NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetGlobalDofValues(int rTimeDerivative, NuTo::FullVector<double, Eigen::Dynamic>& rActiveDofValues, NuTo::FullVector<double, Eigen::Dynamic>& rDependentDofValues) const
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    assert(TNumTimeDerivatives >= rTimeDerivative);

    for (int i = 0; i < TNumDisplacements; ++i)
        WriteNodeValueToVector(this->mDofDisplacements[i], this->mDisplacements[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);

    for (int i = 0; i < TNumRotations; ++i)
        WriteNodeValueToVector(this->mDofRotations[i], this->mRotations[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);

    for (int i = 0; i < TNumTemperature; ++i)
        WriteNodeValueToVector(this->mDofTemperature[i], this->mTemperature[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);

    for (int i = 0; i < TNumNonlocalEqPlasticStrain; ++i)
        WriteNodeValueToVector(this->mDofNonlocalEqPlasticStrain[i], this->mNonlocalEqPlasticStrain[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);

//    for (int i=0; i<TNumNonlocalTotalStrain; ++i)
//        WriteNodeValueToVector(this->mDofNonlocalTotalStrain[i], this->mNonlocalTotalStrain[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);
//
    for (int i = 0; i < TNumNonlocalEqStrain; ++i)
        WriteNodeValueToVector(this->mDofNonlocalEqStrain[i], this->mNonlocalEqStrain[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);

    for (int i = 0; i < TNumWaterVolumeFraction; ++i)
        WriteNodeValueToVector(this->mDofWaterVolumeFraction[i], this->mWaterVolumeFraction[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);

    for (int i = 0; i < TNumRelativeHumidity; ++i)
        WriteNodeValueToVector(this->mDofRelativeHumidity[i], this->mRelativeHumidity[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);
}

//! @brief write dof values to the node (based on global dof number)
//! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
//! @param rDofType ... specific dof type
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
template<NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetGlobalDofValues(
        int rTimeDerivative,
        Node::eDof rDofType,
        const FullVector<double, Eigen::Dynamic>& rActiveDofValues,
        const FullVector<double, Eigen::Dynamic>& rDependentDofValues)
{
    assert(TNumTimeDerivatives >= rTimeDerivative);
    switch (rDofType)
    {
    case Node::DISPLACEMENTS:
        for (int i = 0; i < TNumDisplacements; ++i)
            this->mDisplacements[rTimeDerivative][i] = GetDofValueFromVector(this->mDofDisplacements[i], rActiveDofValues, rDependentDofValues);
        break;

    case Node::ROTATIONS:
        for (int i = 0; i < TNumRotations; ++i)
            this->mRotations[rTimeDerivative][i] = GetDofValueFromVector(this->mDofRotations[i], rActiveDofValues, rDependentDofValues);
        break;

    case Node::TEMPERATURE:
        for (int i = 0; i < TNumTemperature; ++i)
            this->mTemperature[rTimeDerivative][i] = GetDofValueFromVector(this->mDofTemperature[i], rActiveDofValues, rDependentDofValues);
        break;

    case Node::NONLOCALEQPLASTICSTRAIN:
        for (int i = 0; i < TNumNonlocalEqPlasticStrain; ++i)
        {
            double value = GetDofValueFromVector(this->mDofNonlocalEqPlasticStrain[i], rActiveDofValues, rDependentDofValues);
            if (value < 0)
                std::cout << "value is less then zero" << value << std::endl; // no further consequences... lucky you.
        }
        break;

//        case Node::NONLOCALTOTALSTRAIN:
        //    for (int i=0; i<TNumNonlocalTotalStrain; ++i)
        //        this->mNonlocalTotalStrain[rTimeDerivative][i] = GetDofValueFromVector(this->mDofNonlocalTotalStrain[i], rActiveDofValues, rDependentDofValues);
//            break;

    case Node::NONLOCALEQSTRAIN:
        for (int i = 0; i < TNumNonlocalEqStrain; ++i)
            this->mNonlocalEqStrain[rTimeDerivative][i] = GetDofValueFromVector(this->mDofNonlocalEqStrain[i], rActiveDofValues, rDependentDofValues);
        break;

    case Node::WATERVOLUMEFRACTION:
        for (int i = 0; i < TNumWaterVolumeFraction; ++i)
            this->mWaterVolumeFraction[rTimeDerivative][i] = GetDofValueFromVector(this->mDofWaterVolumeFraction[i], rActiveDofValues, rDependentDofValues);
        break;

    case Node::RELATIVEHUMIDITY:
        for (int i = 0; i < TNumRelativeHumidity; ++i)
            this->mRelativeHumidity[rTimeDerivative][i] = GetDofValueFromVector(this->mDofRelativeHumidity[i], rActiveDofValues, rDependentDofValues);
        break;

    default:
        break;
    }
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
//! @param rDofType ... specific dof type
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
template<NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetGlobalDofValues(
        int rTimeDerivative,
        Node::eDof rDofType,
        FullVector<double, Eigen::Dynamic>& rActiveDofValues,
        FullVector<double, Eigen::Dynamic>& rDependentDofValues) const
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    assert(TNumTimeDerivatives >= rTimeDerivative);

    switch (rDofType)
    {
    case Node::DISPLACEMENTS:
        for (int i = 0; i < TNumDisplacements; ++i)
            WriteNodeValueToVector(this->mDofDisplacements[i], this->mDisplacements[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);
        break;

    case Node::ROTATIONS:
        for (int i = 0; i < TNumRotations; ++i)
            WriteNodeValueToVector(this->mDofRotations[i], this->mRotations[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);
        break;

    case Node::TEMPERATURE:
        for (int i = 0; i < TNumTemperature; ++i)
            WriteNodeValueToVector(this->mDofTemperature[i], this->mTemperature[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);
        break;

    case Node::NONLOCALEQPLASTICSTRAIN:
        for (int i = 0; i < TNumNonlocalEqPlasticStrain; ++i)
            WriteNodeValueToVector(this->mDofNonlocalEqPlasticStrain[i], this->mNonlocalEqPlasticStrain[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);
        break;

//        case Node::NONLOCALTOTALSTRAIN:
        //    for (int i=0; i<TNumNonlocalTotalStrain; ++i)
        //        WriteNodeValueToVector(this->mDofNonlocalTotalStrain[i], this->mNonlocalTotalStrain[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);
//            break;

    case Node::NONLOCALEQSTRAIN:
        for (int i = 0; i < TNumNonlocalEqStrain; ++i)
            WriteNodeValueToVector(this->mDofNonlocalEqStrain[i], this->mNonlocalEqStrain[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);
        break;

    case Node::WATERVOLUMEFRACTION:
        for (int i = 0; i < TNumWaterVolumeFraction; ++i)
            WriteNodeValueToVector(this->mDofWaterVolumeFraction[i], this->mWaterVolumeFraction[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);
        break;

    case Node::RELATIVEHUMIDITY:
        for (int i = 0; i < TNumRelativeHumidity; ++i)
            WriteNodeValueToVector(this->mDofRelativeHumidity[i], this->mRelativeHumidity[rTimeDerivative][i], rActiveDofValues, rDependentDofValues);
        break;

    default:
        break;
    }

}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
template<NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    for (int i = 0; i < TNumDisplacements; ++i)
        this->mDofDisplacements[i] = rMappingInitialToNewOrdering[this->mDofDisplacements[i]];

    for (int i = 0; i < TNumRotations; ++i)
        this->mDofRotations[i] = rMappingInitialToNewOrdering[this->mDofRotations[i]];

    for (int i = 0; i < TNumTemperature; ++i)
        this->mDofTemperature[i] = rMappingInitialToNewOrdering[this->mDofTemperature[i]];

    for (int i = 0; i < TNumNonlocalEqPlasticStrain; ++i)
        this->mDofNonlocalEqPlasticStrain[i] = rMappingInitialToNewOrdering[this->mDofNonlocalEqPlasticStrain[i]];
//
//    for (int i=0; i<TNumNonlocalTotalStrain; ++i)
//    	mDofNonlocalTotalStrain[i]=rMappingInitialToNewOrdering[mDofNonlocalTotalStrain[i]];

    for (int i = 0; i < TNumNonlocalEqStrain; ++i)
        this->mDofNonlocalEqStrain[i] = rMappingInitialToNewOrdering[this->mDofNonlocalEqStrain[i]];

    for (int i = 0; i < TNumWaterVolumeFraction; ++i)
        this->mDofWaterVolumeFraction[i] = rMappingInitialToNewOrdering[this->mDofWaterVolumeFraction[i]];

    for (int i = 0; i < TNumRelativeHumidity; ++i)
        this->mDofRelativeHumidity[i] = rMappingInitialToNewOrdering[this->mDofRelativeHumidity[i]];

}

//! @brief renumber the global dofs according to predefined ordering
//! @param rDofType ... specific dof type
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
template<NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
RenumberGlobalDofs(Node::eDof rDofType, std::vector<int>& rMappingInitialToNewOrdering)
{
    switch (rDofType)
    {
    case Node::DISPLACEMENTS:
        for (int i = 0; i < TNumDisplacements; ++i)
            this->mDofDisplacements[i] = rMappingInitialToNewOrdering[this->mDofDisplacements[i]];
        break;

    case Node::ROTATIONS:
        for (int i = 0; i < TNumRotations; ++i)
            this->mDofRotations[i] = rMappingInitialToNewOrdering[this->mDofRotations[i]];
        break;

    case Node::TEMPERATURE:
        for (int i = 0; i < TNumTemperature; ++i)
            this->mDofTemperature[i]=rMappingInitialToNewOrdering[this->mDofTemperature[i]];
        break;

    case Node::NONLOCALEQPLASTICSTRAIN:
        for (int i = 0; i < TNumNonlocalEqPlasticStrain; ++i)
            this->mDofNonlocalEqPlasticStrain[i] = rMappingInitialToNewOrdering[this->mDofNonlocalEqPlasticStrain[i]];
        break;

//        case Node::NONLOCALTOTALSTRAIN:
//            for (int i=0; i<TNumNonlocalTotalStrain; ++i)
//                mDofNonlocalTotalStrain[i]=rMappingInitialToNewOrdering[mDofNonlocalTotalStrain[i]];
//            break;

    case Node::NONLOCALEQSTRAIN:
        for (int i = 0; i < TNumNonlocalEqStrain; ++i)
            this->mDofNonlocalEqStrain[i] = rMappingInitialToNewOrdering[this->mDofNonlocalEqStrain[i]];
        break;

    case Node::WATERVOLUMEFRACTION:
        for (int i = 0; i < TNumWaterVolumeFraction; ++i)
            this->mDofWaterVolumeFraction[i] = rMappingInitialToNewOrdering[this->mDofWaterVolumeFraction[i]];
        break;

    case Node::RELATIVEHUMIDITY:
        for (int i = 0; i < TNumRelativeHumidity; ++i)
            this->mDofRelativeHumidity[i] = rMappingInitialToNewOrdering[this->mDofRelativeHumidity[i]];
        break;

    default:
        break;
    }

}

//! @brief returns the number of displacements of the node
//! @return number of displacements
template<NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNumTimeDerivatives() const
{
    return TNumTimeDerivatives;
}

//! @brief returns the type of node as a string (all the data stored at the node)
//! @return string

template <NODE_DOF_TEMPLATE_PARAMETERS>
std::string NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::GetNodeTypeStr() const
{
    std::stringstream NodeDofype;

    if (TNumTimeDerivatives > 0)
        NodeDofype << "TimeDerivatives:" << TNumTimeDerivatives << "\n";

    if (TNumCoordinates > 0)
        NodeDofype << "Coordinates:" << TNumCoordinates << "\n";

    if (TNumDisplacements > 0)
        NodeDofype << "Displacements:" << TNumDisplacements << "\n";

    if (TNumRotations > 0)
        NodeDofype << "Rotations:" << TNumRotations << "\n";

    if (TNumTemperature > 0)
        NodeDofype << "Temperature:" << TNumTemperature << "\n";

    if (TNumNonlocalEqPlasticStrain > 0)
        NodeDofype << "NonlocalEqPlasticStrain:" << TNumNonlocalEqPlasticStrain << "\n";

    if (TNumNonlocalTotalStrain > 0)
        NodeDofype << "NonlocalTotalStrain:" << TNumNonlocalTotalStrain << "\n";

    if (TNumNonlocalEqStrain > 0)
        NodeDofype << "NonlocalEqStrain:" << TNumNonlocalEqStrain << "\n";

    if (TNumWaterVolumeFraction > 0)
        NodeDofype << "WaterVolumeFraction:" << TNumWaterVolumeFraction << "\n";

    if (TNumRelativeHumidity > 0)
        NodeDofype << "RelativeHumidity:" << TNumRelativeHumidity << "\n";

    return NodeDofype.str();

}

//! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
template<NODE_DOF_TEMPLATE_PARAMETERS>
NuTo::NodeBase* NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
Clone() const
{
    return new NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>(*this);
}

//! @brief extracts the appropriate dof value from the active dof vector or the dependent dof vector
//! @param rDofNumber ... dof number
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
//! @return dof value that corresponds to the rDofNumber
template<NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDofValueFromVector(
        int rDofNumber,
        const FullVector<double, Eigen::Dynamic>& rActiveDofValues,
        const FullVector<double, Eigen::Dynamic>& rDependentDofValues) const
{
    if (rDofNumber < rActiveDofValues.GetNumRows())
    {
        return rActiveDofValues(rDofNumber);
    }
    else
    {
        rDofNumber -= rActiveDofValues.GetNumRows();
        assert(rDofNumber < rDependentDofValues.GetNumRows());
        return rDependentDofValues(rDofNumber);
    }
}

//! @brief writes the appropriate dof value to the active dof vector or the dependent dof vector
//! @param rDofNumber ... dof number
//! @param rDofValue ... dof value
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
template<NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
WriteNodeValueToVector(
        int rDofNumber,
        double rDofValue,
        FullVector<double, Eigen::Dynamic>& rActiveDofValues,
        FullVector<double, Eigen::Dynamic>& rDependentDofValues) const
{
    if (rDofNumber < rActiveDofValues.GetNumRows())
    {
        rActiveDofValues(rDofNumber) = rDofValue;
    }
    else
    {
        rDofNumber -= rActiveDofValues.GetNumRows();
        assert(rDofNumber < rDependentDofValues.GetNumRows());
        rDependentDofValues(rDofNumber) = rDofValue;
    }
}


