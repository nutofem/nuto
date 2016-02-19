// $Id$

#ifndef NODEDOF_H
#define NODEDOF_H
#include "nuto/mechanics/nodes/NodeDof_Def.h"
#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/MechanicsException.h"




//! @brief ... constructor
template <NODE_DOF_TEMPLATE_PARAMETERS>
NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::NodeDof() //: NuTo::NodeBase()
{

}

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetGlobalDofs(int& rDOF)
{
    for (int count=0; count<TNumDisplacements; count++)
    {
    	this->mDofDisplacements[count]=rDOF++;
    }
    for (int count=0; count<TNumRotations; count++)
    {
    	this->mDofRotations[count]=rDOF++;
    }
//    for (int count=0; count<TNumTemperatures; count++)
//    {
//    	mDofTemperatures[count]=rDOF++;
//    }
    for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
    {
        this->mDofNonlocalEqPlasticStrain[count]=rDOF++;
    }
//    for (int count=0; count<TNumNonlocalTotalStrain; count++)
//    {
//    	mDofNonlocalTotalStrain[count]=rDOF++;
//    }
    for (int count=0; count<TNumNonlocalEqStrain; count++)
    {
        this->mDofNonlocalEqStrain[count]=rDOF++;
    }
    for (int count=0; count<TNumWaterVolumeFraction; count++)
    {
        this->mDofWaterVolumeFraction[count]=rDOF++;
    }
    for (int count=0; count<TNumRelativeHumidity; count++)
    {
        this->mDofRelativeHumidity[count]=rDOF++;
    }

}

//! @brief write dof values to the node (based on global dof number)
//! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetGlobalDofValues(int rTimeDerivative, const NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, const NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues)
{
    assert(TNumTimeDerivatives>=rTimeDerivative);

    for (int count=0; count<TNumDisplacements; count++)
    {
        int dof = this->mDofDisplacements[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof);
        }
        else
        {
            value = rActiveDofValues(dof);
        }
        this->mDisplacements[rTimeDerivative][count] = value;
    }

    for (int count=0; count<TNumRotations; count++)
    {
        int dof = this->mDofRotations[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof);
        }
        else
        {
            value = rActiveDofValues(dof);
        }
        this->mRotations[rTimeDerivative][count] = value;
    }
//
//    for (int count=0; count<TNumTemperatures; count++)
//    {
//        int dof = this->mDofTemperatures[count];
//        double value;
//        if (dof >= rActiveDofValues.GetNumRows())
//        {
//            dof -= rActiveDofValues.GetNumRows();
//            assert(dof < rDependentDofValues.GetNumRows());
//            value = rDependentDofValues(dof);
//        }
//        else
//        {
//            value = rActiveDofValues(dof);
//        }
//        this->mTemperatures[rTimeDerivative][count] = value;
//    }
//
    for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
    {
        int dof = this->mDofNonlocalEqPlasticStrain[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof);
        }
        else
        {
            value = rActiveDofValues(dof);
        }
        if (value<0)
        {
        	std::cout << "value is less then zero" << value << std::endl;
        }
        this->mNonlocalEqPlasticStrain[rTimeDerivative][count] = value;
    }
//
//    for (int count=0; count<TNumNonlocalTotalStrain; count++)
//    {
//        int dof = this->mDofNonlocalTotalStrain[count];
//        double value;
//        if (dof >= rActiveDofValues.GetNumRows())
//        {
//            dof -= rActiveDofValues.GetNumRows();
//            assert(dof < rDependentDofValues.GetNumRows());
//            value = rDependentDofValues(dof);
//        }
//        else
//        {
//            value = rActiveDofValues(dof);
//        }
//        this->mNonlocalTotalStrain[rTimeDerivative][count] = value;
//    }
//
    for (int count=0; count<TNumNonlocalEqStrain; count++)
    {
        int dof = this->mDofNonlocalEqStrain[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof);
        }
        else
        {
            value = rActiveDofValues(dof);
        }
        this->mNonlocalEqStrain[rTimeDerivative][count] = value;
    }

    for (int count=0; count<TNumWaterVolumeFraction; count++)
    {
        int dof = this->mDofWaterVolumeFraction[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof);
        }
        else
        {
            value = rActiveDofValues(dof);
        }
        this->mWaterVolumeFraction[rTimeDerivative][count] = value;
    }

    for (int count=0; count<TNumRelativeHumidity; count++)
    {
        int dof = this->mDofRelativeHumidity[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof);
        }
        else
        {
            value = rActiveDofValues(dof);
        }
        this->mRelativeHumidity[rTimeDerivative][count] = value;
    }

}


//! @brief extract dof values from the node (based on global dof number)
//! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetGlobalDofValues(int rTimeDerivative, NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues) const
{
	assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);
    assert(TNumTimeDerivatives>=rTimeDerivative);

    for (int count=0; count<TNumDisplacements; count++)
    {
        int dof = this->mDofDisplacements[count];
        double value = this->mDisplacements[rTimeDerivative][count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof) = value;
        }
        else
        {
            rActiveDofValues(dof) = value;
        }
    }

    for (int count=0; count<TNumRotations; count++)
    {
        int dof = this->mDofRotations[count];
        double value = this->mRotations[rTimeDerivative][count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof) = value;
        }
        else
        {
            rActiveDofValues(dof) = value;
        }
    }
//
//    for (int count=0; count<TNumTemperatures; count++)
//    {
//        int dof = this->mDofTemperatures[count];
//        double value = this->mTemperatures[rTimeDerivative][count];
//        if (dof >= rActiveDofValues.GetNumRows())
//        {
//            dof -= rActiveDofValues.GetNumRows();
//            assert(dof < rDependentDofValues.GetNumRows());
//            rDependentDofValues(dof) = value;
//        }
//        else
//        {
//            rActiveDofValues(dof) = value;
//        }
//    }
//
    for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
    {
        int dof = this->mDofNonlocalEqPlasticStrain[count];
        double value = this->mNonlocalEqPlasticStrain[rTimeDerivative][count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof) = value;
        }
        else
        {
            rActiveDofValues(dof) = value;
        }
    }
//
//    for (int count=0; count<TNumNonlocalTotalStrain; count++)
//    {
//        int dof = this->mDofNonlocalTotalStrain[count];
//        double value = this->mNonlocalTotalStrain[rTimeDerivative][count];
//        if (dof >= rActiveDofValues.GetNumRows())
//        {
//            dof -= rActiveDofValues.GetNumRows();
//            assert(dof < rDependentDofValues.GetNumRows());
//            rDependentDofValues(dof) = value;
//        }
//        else
//        {
//            rActiveDofValues(dof) = value;
//        }
//    }
//
    for (int count=0; count<TNumNonlocalEqStrain; count++)
    {
        int dof = this->mDofNonlocalEqStrain[count];
        double value = this->mNonlocalEqStrain[rTimeDerivative][count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof) = value;
        }
        else
        {
            rActiveDofValues(dof) = value;
        }
    }

    for (int count=0; count<TNumWaterVolumeFraction; count++)
    {
        int dof = this->mDofWaterVolumeFraction[count];
        double value = this->mWaterVolumeFraction[rTimeDerivative][count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof) = value;
        }
        else
        {
            rActiveDofValues(dof) = value;
        }
    }

    for (int count=0; count<TNumRelativeHumidity; count++)
    {
        int dof = this->mDofRelativeHumidity[count];
        double value = this->mRelativeHumidity[rTimeDerivative][count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof) = value;
        }
        else
        {
            rActiveDofValues(dof) = value;
        }
    }
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    for (int count=0; count<TNumDisplacements; count++)
    {
    	this->mDofDisplacements[count]=rMappingInitialToNewOrdering[this->mDofDisplacements[count]];
    }

    for (int count=0; count<TNumRotations; count++)
    {
    	this->mDofRotations[count]=rMappingInitialToNewOrdering[this->mDofRotations[count]];
    }

//    for (int count=0; count<TNumTemperatures; count++)
//    {
//    	this->mDofTemperatures[count]=rMappingInitialToNewOrdering[this->mDofTemperatures[count]];
//    }
//
    for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
    {
        this->mDofNonlocalEqPlasticStrain[count]=rMappingInitialToNewOrdering[this->mDofNonlocalEqPlasticStrain[count]];
    }
//
//    for (int count=0; count<TNumNonlocalTotalStrain; count++)
//    {
//    	mDofNonlocalTotalStrain[count]=rMappingInitialToNewOrdering[mDofNonlocalTotalStrain[count]];
//    }
    for (int count=0; count<TNumNonlocalEqStrain; count++)
    {
        this->mDofNonlocalEqStrain[count]=rMappingInitialToNewOrdering[this->mDofNonlocalEqStrain[count]];
    }
    for (int count=0; count<TNumWaterVolumeFraction; count++)
    {
        this->mDofWaterVolumeFraction[count]=rMappingInitialToNewOrdering[this->mDofWaterVolumeFraction[count]];
    }
    for (int count=0; count<TNumRelativeHumidity; count++)
    {
        this->mDofRelativeHumidity[count]=rMappingInitialToNewOrdering[this->mDofRelativeHumidity[count]];
    }
}

//! @brief returns the number of displacements of the node
//! @return number of displacements
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNumTimeDerivatives()const
{
	return TNumTimeDerivatives;
}

//
//
////*************************************************
////************      TEMPERATURE     ***************
////*************************************************
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNumTemperatures()const
//{
//	return TNumTemperatures;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetDofTemperature()const
//{
//    assert(TNumTemperatures==1);
//    return mDofTemperatures[0];
//}
//
////*************************************************
////************      TEMPERATURE  GET   ************
////*************************************************
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetTemperature()const
//{
//	assert(TNumTemperatures==1);
//	return mTemperatures[0][0];
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetTemperature(int rTimeDerivative)const
//{
//	assert(TNumTemperatures==1);
//	assert(rTimeDerivative>=0);
//	assert(rTimeDerivative<=TNumTimeDerivatives);
//	return mTemperatures[rTimeDerivative][0];
//}
//
////*************************************************
////************      TEMPERATURE  SET   ************
////*************************************************
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetTemperature(double rTemperatures)
//{
//	assert(TNumTemperatures==1);
//	mTemperatures[0][0] = rTemperatures;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetTemperature(int rTimeDerivative, double rTemperatures)
//{
//	assert(TNumTemperatures==1);
//	assert(rTimeDerivative>=0);
//	assert(rTimeDerivative<=TNumTimeDerivatives);
//	mTemperatures[rTimeDerivative][0] = rTemperatures;
//}
//
//

//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetNonlocalEqPlasticStrain(const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain)
//{
//	assert(TNumNonlocalEqPlasticStrain==2);
//	mNonlocalEqPlasticStrain[0] = rNonlocalEqPlasticStrain;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetNonlocalEqPlasticStrain(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain)
//{
//	assert(TNumNonlocalEqPlasticStrain==2);
//	assert(rTimeDerivative>=0);
//	assert(rTimeDerivative<=TNumTimeDerivatives);
//    mNonlocalEqPlasticStrain[rTimeDerivative] = rNonlocalEqPlasticStrain;
//}
//
////*************************************************
////********    NONLOCAL TOTAL STRAIN     ***********
////*************************************************
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNumNonlocalTotalStrain()const
//{
//	return TNumNonlocalTotalStrain;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetDofNonlocalTotalStrain(int rComponent)const
//{
//    assert(rComponent<TNumNonlocalTotalStrain);
//    return mDofNonlocalTotalStrain[rComponent];
//
//}
//
////*************************************************
////*******    NONLOCAL TOTAL STRAIN  GET  **********
////*************************************************
//
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNonlocalTotalStrain(short rIndex)const
//{
//    assert(rIndex<TNumNonlocalTotalStrain);
//    return mNonlocalTotalStrain[0][rIndex];
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//const Eigen::Matrix<double, 1, 1>& NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNonlocalTotalStrain1D() const
//{
//    assert(TNumNonlocalTotalStrain == 1);
//    return mNonlocalTotalStrain[0];
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//const Eigen::Matrix<double, 3, 1>& NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNonlocalTotalStrain2D() const
//{
//    assert(TNumNonlocalTotalStrain == 2);
//    return mNonlocalTotalStrain[0];
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//const Eigen::Matrix<double, 6, 1>& NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNonlocalTotalStrain3D() const
//{
//    assert(TNumNonlocalTotalStrain == 3);
//    return mNonlocalTotalStrain[0];
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//const Eigen::Matrix<double, Eigen::Dynamic, 1>& NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNonlocalTotalStrains() const
//{
//    return mNonlocalTotalStrain[0];
//}
//
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//const Eigen::Matrix<double, 1, 1>& NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNonlocalTotalStrain1D(int rTimeDerivative) const
//{
//    assert(TNumNonlocalTotalStrain==1);
//    assert(rTimeDerivative>=0);
//    assert(rTimeDerivative<=TNumTimeDerivatives);
//    return mNonlocalTotalStrain[rTimeDerivative];
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//const Eigen::Matrix<double, 3, 1>& NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNonlocalTotalStrain2D(int rTimeDerivative) const
//{
//    assert(TNumNonlocalTotalStrain==2);
//    assert(rTimeDerivative>=0);
//    assert(rTimeDerivative<=TNumTimeDerivatives);
//    return mNonlocalTotalStrain[rTimeDerivative];
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//const Eigen::Matrix<double, 6, 1>& NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNonlocalTotalStrain3D(int rTimeDerivative) const
//{
//    assert(TNumNonlocalTotalStrain==3);
//    assert(rTimeDerivative>=0);
//    assert(rTimeDerivative<=TNumTimeDerivatives);
//    return mNonlocalTotalStrain[rTimeDerivative];
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//const Eigen::Matrix<double, Eigen::Dynamic, 1>& NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNonlocalTotalStrains(int rTimeDerivative) const
//{
//    assert(rTimeDerivative<=TNumTimeDerivatives);
//    return mNonlocalTotalStrain[rTimeDerivative];
//}
//
////*************************************************
////*******    NONLOCAL TOTAL STRAIN  SET  **********
////*************************************************
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetNonlocalTotalStrain1D(const Eigen::Matrix<double, 1, 1>& rNonlocalTotalStrain)
//{
//    assert(TNumNonlocalTotalStrain == 1);
//    mNonlocalTotalStrain[0] = rNonlocalTotalStrain;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetNonlocalTotalStrain2D(const Eigen::Matrix<double, 3, 1>& rNonlocalTotalStrain)
//{
//    assert(TNumNonlocalTotalStrain == 2);
//    mNonlocalTotalStrain[0] = rNonlocalTotalStrain;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetNonlocalTotalStrain3D(const Eigen::Matrix<double, 6, 1>& rNonlocalTotalStrain)
//{
//    assert(TNumNonlocalTotalStrain == 3);
//    mNonlocalTotalStrain[0] = rNonlocalTotalStrain;
//}
//
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetNonlocalTotalStrain1D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rNonlocalTotalStrain)
//{
//    assert(TNumNonlocalTotalStrain == 1);
//    assert(rTimeDerivative>=0);
//    assert(rTimeDerivative<=TNumTimeDerivatives);
//    mNonlocalTotalStrain[rTimeDerivative] = rNonlocalTotalStrain;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetNonlocalTotalStrain2D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rNonlocalTotalStrain)
//{
//    assert(TNumNonlocalTotalStrain == 2);
//    assert(rTimeDerivative>=0);
//    assert(rTimeDerivative<=TNumTimeDerivatives);
//    mNonlocalTotalStrain[rTimeDerivative] = rNonlocalTotalStrain;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetNonlocalTotalStrain3D(int rTimeDerivative, const Eigen::Matrix<double, 6, 1>& rNonlocalTotalStrain)
//{
//    assert(TNumNonlocalTotalStrain == 3);
//    assert(rTimeDerivative>=0);
//    assert(rTimeDerivative<=TNumTimeDerivatives);
//    mNonlocalTotalStrain[rTimeDerivative] = rNonlocalTotalStrain;
//}
//

//
////**************************************************
////*******      WATER VOLUME FRACTION      **********
////**************************************************
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNumWaterVolumeFraction()const
//{
//    return TNumWaterVolumeFraction;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetDofWaterVolumeFraction()const
//{
//    assert(TNumWaterVolumeFraction==1);
//    return mDofWaterVolumeFraction[0];
//}
//
////**************************************************
////*******      WATER VOLUME FRACTION  GET  *********
////**************************************************
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetWaterVolumeFraction()const
//{
//    assert(TNumWaterVolumeFraction==1);
//    return mWaterVolumeFraction[0][0];
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetWaterVolumeFraction(int rTimeDerivative)const
//{
//    assert(TNumWaterVolumeFraction==1);
//    assert(rTimeDerivative>=0 && rTimeDerivative<=TNumTimeDerivatives);
//    return mWaterVolumeFraction[rTimeDerivative][0];
//}
//
////**************************************************
////*******      WATER VOLUME FRACTION  SET  *********
////**************************************************
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetWaterVolumeFraction(double rWaterVolumeFraction)
//{
//    assert(TNumWaterVolumeFraction==1);
//    mWaterVolumeFraction[0][0] = rWaterVolumeFraction;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetWaterVolumeFraction(int rTimeDerivative, double rWaterVolumeFraction)
//{
//    assert(TNumWaterVolumeFraction==1);
//    assert(rTimeDerivative>=0 && rTimeDerivative<=TNumTimeDerivatives);
//    mWaterVolumeFraction[0][0] = rWaterVolumeFraction;
//}
//
////*************************************************
////*******       RELATIVE HUMIDITY        **********
////*************************************************
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetNumRelativeHumidity()const
//{
//    return TNumRelativeHumidity;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetDofRelativeHumidity()const
//{
//    assert(TNumRelativeHumidity==1);
//    return mDofRelativeHumidity[0];
//}
//
////*************************************************
////*******       RELATIVE HUMIDITY  GET   **********
////*************************************************
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetRelativeHumidity()const
//{
//    assert(TNumRelativeHumidity==1);
//    return mRelativeHumidity[0][0];
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//GetRelativeHumidity(int rTimeDerivative)const
//{
//    assert(TNumRelativeHumidity==1);
//    assert(rTimeDerivative>=0 && rTimeDerivative<=TNumTimeDerivatives);
//    return mRelativeHumidity[rTimeDerivative][0];
//}
//
////*************************************************
////*******       RELATIVE HUMIDITY  SET   **********
////*************************************************
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetRelativeHumidity(double rRelativeHumidity)
//{
//    assert(TNumRelativeHumidity==1);
//    mRelativeHumidity[0][0] = rRelativeHumidity;
//}
//
//template <NODE_DOF_TEMPLATE_PARAMETERS>
//void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
//SetRelativeHumidity(int rTimeDerivative, double rRelativeHumidity)
//{
//    assert(TNumRelativeHumidity==1);
//    assert(rTimeDerivative>=0 && rTimeDerivative<=TNumTimeDerivatives);
//    mRelativeHumidity[0][0] = rRelativeHumidity;
//}
//



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

    if (TNumTemperatures > 0)
        NodeDofype << "Temperatures:" << TNumTemperatures << "\n";

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
template <NODE_DOF_TEMPLATE_PARAMETERS>
NuTo::NodeBase* NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
Clone()const
{
    return new NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>(*this);
}

#endif //NodeDof_H


