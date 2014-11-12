// $Id$

#ifndef NODE_DOF_H
#define NODE_DOF_H
#include "nuto/mechanics/nodes/NodeDof_Def.h"
#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/MechanicsException.h"




//! @brief ... constructor
template <NODE_DOF_TEMPLATE_PARAMETERS>
NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::NodeDof() //: NuTo::NodeBase()
{
	for (int countDerivatives=0; countDerivatives<TNumTimeDerivatives+1; countDerivatives++)
	{
		for (int count=0; count<TNumDisplacements; count++)
		{
			mDisplacements[countDerivatives][count]=0.;
		}
		for (int count=0; count<TNumRotations; count++)
		{
			mRotations[countDerivatives][count]=0.;
		}
		for (int count=0; count<TNumTemperatures; count++)
		{
			mTemperatures[countDerivatives][count]=0.;
		}
		for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
		{
			mNonlocalEqPlasticStrain[countDerivatives][count]=0.;
		}
		for (int count=0; count<TNumNonlocalTotalStrain; count++)
		{
			mNonlocalTotalStrain[countDerivatives][count]=0.;
		}
        for (int count=0; count<TNumNonlocalEqStrain; count++)
        {
            mNonlocalEqStrain[countDerivatives][count]=0.;
        }
	}
	for (int count=0; count<TNumDisplacements; count++)
	{
		mDofDisplacements[count]=-1;
	}
	for (int count=0; count<TNumRotations; count++)
	{
		mDofRotations[count]=-1;
	}
	for (int count=0; count<TNumTemperatures; count++)
	{
		mDofTemperatures[count]=-1;
	}
	for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
	{
		mDofNonlocalEqPlasticStrain[count]=-1;
	}
	for (int count=0; count<TNumNonlocalTotalStrain; count++)
	{
		mDofNonlocalTotalStrain[count]=-1;
	}
    for (int count=0; count<TNumNonlocalEqStrain; count++)
    {
        mDofNonlocalEqStrain[count]=-1;
    }


}

//! @brief ... destructor
template <NODE_DOF_TEMPLATE_PARAMETERS>
NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::~NodeDof() //: NuTo::~NodeBase()
{
}

//! @brief assignment operator
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
     operator=(NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION> const& rOther)
{
    mDisplacements = rOther.mDisplacements;
    mDofDisplacements = rOther.mDofDisplacements;

    mRotations = rOther.mRotations;
    mDofRotations = rOther.mDofRotations;

    mTemperatures = rOther.mTemperatures;
    mDofTemperatures = rOther.mDofTemperatures;

    mNonlocalEqPlasticStrain = rOther.mNonlocalEqPlasticStrain;
    mDofNonlocalEqPlasticStrain = rOther.mDofNonlocalEqPlasticStrain;

    mNonlocalTotalStrain = rOther.mNonlocalTotalStrain;
    mDofNonlocalTotalStrain = rOther.mDofNonlocalTotalStrain;

    mNonlocalEqStrain = rOther.mNonlocalEqStrain;
    mDofNonlocalEqStrain = rOther.mDofNonlocalEqStrain;

    mWaterPhaseFraction = rOther.mWaterPhaseFraction;
    mDofWaterPhaseFraction = rOther.mDofWaterPhaseFraction;

    mRelativeHumidity = rOther.mRelativeHumidity;
    mDofRelativeHumidity=rOther.mDofRelativeHumidity;

}

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetGlobalDofs(int& rDOF)
{
    for (int count=0; count<TNumDisplacements; count++)
    {
    	mDofDisplacements[count]=rDOF++;
    }
    for (int count=0; count<TNumRotations; count++)
    {
    	mDofRotations[count]=rDOF++;
    }
    for (int count=0; count<TNumTemperatures; count++)
    {
    	mDofTemperatures[count]=rDOF++;
    }
    for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
    {
    	mDofNonlocalEqPlasticStrain[count]=rDOF++;
    }
    for (int count=0; count<TNumNonlocalTotalStrain; count++)
    {
    	mDofNonlocalTotalStrain[count]=rDOF++;
    }
    for (int count=0; count<TNumNonlocalEqStrain; count++)
    {
        mDofNonlocalEqStrain[count]=rDOF++;
    }
    for (int count=0; count<TNumWaterPhaseFraction; count++)
    {
        mDofWaterPhaseFraction[count]=rDOF++;
    }
    for (int count=0; count<TNumRelativeHumidity; count++)
    {
        mDofRelativeHumidity[count]=rDOF++;
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

    for (int count=0; count<TNumTemperatures; count++)
    {
        int dof = this->mDofTemperatures[count];
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
        this->mTemperatures[rTimeDerivative][count] = value;
    }

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

    for (int count=0; count<TNumNonlocalTotalStrain; count++)
    {
        int dof = this->mDofNonlocalTotalStrain[count];
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
        this->mNonlocalTotalStrain[rTimeDerivative][count] = value;
    }

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

    for (int count=0; count<TNumWaterPhaseFraction; count++)
    {
        int dof = this->mDofWaterPhaseFraction[count];
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
        this->mWaterPhaseFraction[rTimeDerivative][count] = value;
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

    for (int count=0; count<TNumTemperatures; count++)
    {
        int dof = this->mDofTemperatures[count];
        double value = this->mTemperatures[rTimeDerivative][count];
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

    for (int count=0; count<TNumNonlocalTotalStrain; count++)
    {
        int dof = this->mDofNonlocalTotalStrain[count];
        double value = this->mNonlocalTotalStrain[rTimeDerivative][count];
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

    for (int count=0; count<TNumWaterPhaseFraction; count++)
    {
        int dof = this->mDofWaterPhaseFraction[count];
        double value = this->mWaterPhaseFraction[rTimeDerivative][count];
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
    	mDofDisplacements[count]=rMappingInitialToNewOrdering[mDofDisplacements[count]];
    }

    for (int count=0; count<TNumRotations; count++)
    {
    	mDofRotations[count]=rMappingInitialToNewOrdering[mDofRotations[count]];
    }

    for (int count=0; count<TNumTemperatures; count++)
    {
    	mDofTemperatures[count]=rMappingInitialToNewOrdering[mDofTemperatures[count]];
    }

    for (int count=0; count<TNumNonlocalEqPlasticStrain; count++)
    {
    	mDofNonlocalEqPlasticStrain[count]=rMappingInitialToNewOrdering[mDofNonlocalEqPlasticStrain[count]];
    }

    for (int count=0; count<TNumNonlocalTotalStrain; count++)
    {
    	mDofNonlocalTotalStrain[count]=rMappingInitialToNewOrdering[mDofNonlocalTotalStrain[count]];
    }
    for (int count=0; count<TNumNonlocalEqStrain; count++)
    {
        mDofNonlocalEqStrain[count]=rMappingInitialToNewOrdering[mDofNonlocalEqStrain[count]];
    }
    for (int count=0; count<TNumNonlocalEqStrain; count++)
    {
        mDofWaterPhaseFraction[count]=rMappingInitialToNewOrdering[mDofWaterPhaseFraction[count]];
    }
    for (int count=0; count<TNumNonlocalEqStrain; count++)
    {
        mDofRelativeHumidity[count]=rMappingInitialToNewOrdering[mDofRelativeHumidity[count]];
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

//! @brief returns the number of displacements of the node
//! @return number of displacements
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNumDisplacements()const
{
	return TNumDisplacements;
}

//! @brief gives the global DOF of a displacement component
//! @param rComponent component
//! @return global DOF
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDofDisplacement(int rComponent)const
{
	assert(rComponent<TNumDisplacements);
	return mDofDisplacements[rComponent];
}

//! @brief returns the displacements of the node
//! @return displacement
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDisplacements1D(double rDisplacements[1])const
{
	assert(TNumDisplacements==1);
	rDisplacements[0] = mDisplacements[0][0];
}

//! @brief returns the displacements of the node
//! @return displacement
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDisplacements1D(int rTimeDerivative, double rDisplacements[1])const
{
	assert(TNumDisplacements==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rDisplacements[0] = mDisplacements[rTimeDerivative][0];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetDisplacements1D(const double rDisplacements[1])
{
	assert(TNumDisplacements==1);
	mDisplacements[0][0] = rDisplacements[0];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetDisplacements1D(int rTimeDerivative, const double rDisplacements[1])
{
	assert(TNumDisplacements==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mDisplacements[rTimeDerivative][0] = rDisplacements[0];
}

//! @brief returns the displacements of the node
//! @return displacement
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDisplacements2D(double rDisplacements[2])const
{
	assert(TNumDisplacements==2);
	rDisplacements[0] = mDisplacements[0][0];
	rDisplacements[1] = mDisplacements[0][1];
}

//! @brief returns the displacements of the node
//! @return displacement
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDisplacements2D(int rTimeDerivative, double rDisplacements[2])const
{
	assert(TNumDisplacements==2);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rDisplacements[0] = mDisplacements[rTimeDerivative][0];
	rDisplacements[1] = mDisplacements[rTimeDerivative][1];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetDisplacements2D(const double rDisplacements[2])
{
	assert(TNumDisplacements==2);
	mDisplacements[0][0] = rDisplacements[0];
	mDisplacements[0][1] = rDisplacements[1];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetDisplacements2D(int rTimeDerivative, const double rDisplacements[2])
{
	assert(TNumDisplacements==2);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mDisplacements[rTimeDerivative][0] = rDisplacements[0];
	mDisplacements[rTimeDerivative][1] = rDisplacements[1];
}

//! @brief returns the displacements of the node
//! @return displacement
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDisplacements3D(double rDisplacements[3])const
{
	assert(TNumDisplacements==3);
	rDisplacements[0] = mDisplacements[0][0];
	rDisplacements[1] = mDisplacements[0][1];
	rDisplacements[2] = mDisplacements[0][2];
}

//! @brief returns the displacements of the node
//! @return displacement
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDisplacements3D(int rTimeDerivative, double rDisplacements[3])const
{
	assert(TNumDisplacements==3);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rDisplacements[0] = mDisplacements[rTimeDerivative][0];
	rDisplacements[1] = mDisplacements[rTimeDerivative][1];
	rDisplacements[2] = mDisplacements[rTimeDerivative][2];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetDisplacements3D(const double rDisplacements[3])
{
	assert(TNumDisplacements==3);
	mDisplacements[0][0] = rDisplacements[0];
	mDisplacements[0][1] = rDisplacements[1];
	mDisplacements[0][2] = rDisplacements[2];
}

//! @brief set the displacements
//! @param rDisplacements  given displacements
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetDisplacements3D(int rTimeDerivative, const double rDisplacements[3])
{
	assert(TNumDisplacements==3);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mDisplacements[rTimeDerivative][0] = rDisplacements[0];
	mDisplacements[rTimeDerivative][1] = rDisplacements[1];
	mDisplacements[rTimeDerivative][2] = rDisplacements[2];
}

//! @brief returns the displacements of the node
//! @return displacement
template <NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDisplacement(short rIndex)const
{
	assert(rIndex<TNumDisplacements);
	return mDisplacements[0][rIndex];
}


//! @brief returns the number of Rotations of the node
//! @return number of Rotations
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNumRotations()const
{
	return TNumRotations;
}

//! @brief gives the global DOF of a Rotation component
//! @param rComponent component
//! @return global DOF
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDofRotation(int rComponent)const
{
	assert(TNumRotations>rComponent);
	return mDofRotations[rComponent];
}

//! @brief returns the Rotations of the node
//! @return Rotation
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetRotations2D(double rRotations[1])const
{
	assert(TNumRotations==1);
	rRotations[0] = mRotations[0][0];
}

//! @return Rotation
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetRotations2D(int rTimeDerivative, double rRotations[1])const
{
	assert(TNumRotations==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rRotations[rTimeDerivative] = mRotations[0][0];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetRotations2D(const double rRotations[1])
{
	assert(TNumRotations==1);
	mRotations[0][0] = rRotations[0];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetRotations2D(int rTimeDerivative, const double rRotations[1])
{
	assert(TNumRotations==1);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mRotations[rTimeDerivative][0] = rRotations[0];
}

//! @brief returns the Rotations of the node
//! @return Rotation
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetRotations3D(double rRotations[3])const
{
	assert(TNumRotations==3);
	rRotations[0] = mRotations[0][0];
	rRotations[1] = mRotations[0][1];
	rRotations[2] = mRotations[0][2];
}

//! @brief returns the Rotations of the node
//! @return Rotation
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetRotations3D(int rTimeDerivative, double rRotations[3])const
{
	assert(TNumRotations==3);
	rRotations[0] = mRotations[rTimeDerivative][0];
	rRotations[1] = mRotations[rTimeDerivative][1];
	rRotations[2] = mRotations[rTimeDerivative][2];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetRotations3D(const double rRotations[3])
{
	assert(TNumRotations==3);
	mRotations[0][0] = rRotations[0];
	mRotations[0][1] = rRotations[1];
	mRotations[0][2] = rRotations[2];
}

//! @brief set the Rotations
//! @param rRotations  given Rotations
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetRotations3D(int rTimeDerivative, const double rRotations[3])
{
	assert(TNumRotations==3);
	mRotations[rTimeDerivative][0] = rRotations[0];
	mRotations[rTimeDerivative][1] = rRotations[1];
	mRotations[rTimeDerivative][2] = rRotations[2];
}


//! @brief returns the Rotations of the node
//! @return Rotation
template <NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetRotation(short rIndex)const
{
	assert(TNumRotations>rIndex);
	return mRotations[0][rIndex];
}

//! @brief returns the number of temperatures of the node
//! @return number of temperatures
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNumTemperatures()const
{
	return TNumTemperatures;
}

//! @brief returns the temperature of the node
//! @return temperature
template <NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetTemperature()const
{
	assert(TNumTemperatures==1);
	return mTemperatures[0][0];
}

//! @brief returns the temperature of the node
//! @return temperature
template <NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetTemperature(int rTimeDerivative)const
{
	assert(TNumTemperatures==1);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	return mTemperatures[rTimeDerivative][0];
}

//! @brief set the temperature of the node
//! @param rTemperature  given temperature
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetTemperature(double rTemperatures)
{
	assert(TNumTemperatures==1);
	mTemperatures[0][0] = rTemperatures;
}

//! @brief set the temperature of the node
//! @param rTemperature  given temperature
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetTemperature(int rTimeDerivative, double rTemperatures)
{
	assert(TNumTemperatures==1);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mTemperatures[rTimeDerivative][0] = rTemperatures;
}

//! @brief gives the global DOF of a temperature component
//! @param rComponent component
//! @return global DOF
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDofTemperature()const
{
	assert(TNumTemperatures==1);
	return mDofTemperatures[0];

}

//! @brief returns the number of Damage of the node
//! @return number of Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNumNonlocalEqPlasticStrain()const
{
	return TNumNonlocalEqPlasticStrain;
}

//! @brief returns the Damage of the node
//! @return Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNonlocalEqPlasticStrain(double* rNonlocalEqPlasticStrain)const
{
	assert(TNumNonlocalEqPlasticStrain==2);
	rNonlocalEqPlasticStrain[0] = mNonlocalEqPlasticStrain[0][0];
	rNonlocalEqPlasticStrain[1] = mNonlocalEqPlasticStrain[0][1];
}

//! @brief returns the Damage of the node
//! @return Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNonlocalEqPlasticStrain(int rTimeDerivative, double* rNonlocalEqPlasticStrain)const
{
	assert(TNumNonlocalEqPlasticStrain==2);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rNonlocalEqPlasticStrain[0] = mNonlocalEqPlasticStrain[rTimeDerivative][0];
	rNonlocalEqPlasticStrain[1] = mNonlocalEqPlasticStrain[rTimeDerivative][1];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetNonlocalEqPlasticStrain(const double* rNonlocalEqPlasticStrain)
{
	assert(TNumNonlocalEqPlasticStrain==2);
	mNonlocalEqPlasticStrain[0][0] = rNonlocalEqPlasticStrain[0];
	mNonlocalEqPlasticStrain[0][1] = rNonlocalEqPlasticStrain[1];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetNonlocalEqPlasticStrain(int rTimeDerivative, const double* rNonlocalEqPlasticStrain)
{
	assert(TNumNonlocalEqPlasticStrain==2);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mNonlocalEqPlasticStrain[rTimeDerivative][0] = rNonlocalEqPlasticStrain[0];
	mNonlocalEqPlasticStrain[rTimeDerivative][1] = rNonlocalEqPlasticStrain[1];
}

//! @brief gives the global DOF of a Damage component
//! @param rComponent component
//! @return global DOF
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDofNonlocalEqPlasticStrain(int rComponent)const
{
	assert(rComponent<TNumNonlocalEqPlasticStrain);
	return mDofNonlocalEqPlasticStrain[rComponent];

}

//! @brief returns the number of Damage of the node
//! @return number of Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNumNonlocalTotalStrain()const
{
	return TNumNonlocalTotalStrain;
}

//! @brief returns the Damage of the node
//! @return Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNonlocalTotalStrain1D(double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==1);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[0][0];
}

//! @brief returns the Damage of the node
//! @return Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNonlocalTotalStrain2D(double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==3);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[0][0];
	rNonlocalTotalStrain[1] = mNonlocalTotalStrain[0][1];
	rNonlocalTotalStrain[2] = mNonlocalTotalStrain[0][2];
}

//! @brief returns the Damage of the node
//! @return Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNonlocalTotalStrain3D(double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==6);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[0][0];
	rNonlocalTotalStrain[1] = mNonlocalTotalStrain[0][1];
	rNonlocalTotalStrain[2] = mNonlocalTotalStrain[0][2];
	rNonlocalTotalStrain[3] = mNonlocalTotalStrain[0][3];
	rNonlocalTotalStrain[4] = mNonlocalTotalStrain[0][4];
	rNonlocalTotalStrain[5] = mNonlocalTotalStrain[0][5];
}

//! @brief returns the Damage of the node
//! @return Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNonlocalTotalStrain1D(int rTimeDerivative, double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==1);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[rTimeDerivative][0];
}

//! @brief returns the Damage of the node
//! @return Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNonlocalTotalStrain2D(int rTimeDerivative, double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==3);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[rTimeDerivative][0];
	rNonlocalTotalStrain[1] = mNonlocalTotalStrain[rTimeDerivative][1];
	rNonlocalTotalStrain[1] = mNonlocalTotalStrain[rTimeDerivative][2];
}

//! @brief returns the Damage of the node
//! @return Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNonlocalTotalStrain3D(int rTimeDerivative, double* rNonlocalTotalStrain)const
{
	assert(TNumNonlocalTotalStrain==6);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	rNonlocalTotalStrain[0] = mNonlocalTotalStrain[rTimeDerivative][0];
	rNonlocalTotalStrain[1] = mNonlocalTotalStrain[rTimeDerivative][1];
	rNonlocalTotalStrain[2] = mNonlocalTotalStrain[rTimeDerivative][2];
	rNonlocalTotalStrain[3] = mNonlocalTotalStrain[rTimeDerivative][3];
	rNonlocalTotalStrain[4] = mNonlocalTotalStrain[rTimeDerivative][4];
	rNonlocalTotalStrain[5] = mNonlocalTotalStrain[rTimeDerivative][5];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetNonlocalTotalStrain1D(const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==1);
	mNonlocalTotalStrain[0][0] = rNonlocalTotalStrain[0];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetNonlocalTotalStrain2D(const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==3);
	mNonlocalTotalStrain[0][0] = rNonlocalTotalStrain[0];
	mNonlocalTotalStrain[0][1] = rNonlocalTotalStrain[1];
	mNonlocalTotalStrain[0][2] = rNonlocalTotalStrain[2];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetNonlocalTotalStrain3D(const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==6);
	mNonlocalTotalStrain[0][0] = rNonlocalTotalStrain[0];
	mNonlocalTotalStrain[0][1] = rNonlocalTotalStrain[1];
	mNonlocalTotalStrain[0][2] = rNonlocalTotalStrain[2];
	mNonlocalTotalStrain[0][3] = rNonlocalTotalStrain[3];
	mNonlocalTotalStrain[0][4] = rNonlocalTotalStrain[4];
	mNonlocalTotalStrain[0][5] = rNonlocalTotalStrain[5];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetNonlocalTotalStrain1D(int rTimeDerivative, const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==1);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mNonlocalTotalStrain[rTimeDerivative][0] = rNonlocalTotalStrain[0];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetNonlocalTotalStrain2D(int rTimeDerivative, const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==3);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mNonlocalTotalStrain[rTimeDerivative][0] = rNonlocalTotalStrain[0];
	mNonlocalTotalStrain[rTimeDerivative][1] = rNonlocalTotalStrain[1];
	mNonlocalTotalStrain[rTimeDerivative][2] = rNonlocalTotalStrain[2];
}

//! @brief set the Damage of the node
//! @param rDamage  given Damage
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetNonlocalTotalStrain3D(int rTimeDerivative, const double* rNonlocalTotalStrain)
{
	assert(TNumNonlocalTotalStrain==6);
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<=TNumTimeDerivatives);
	mNonlocalTotalStrain[rTimeDerivative][0] = rNonlocalTotalStrain[0];
	mNonlocalTotalStrain[rTimeDerivative][1] = rNonlocalTotalStrain[1];
	mNonlocalTotalStrain[rTimeDerivative][2] = rNonlocalTotalStrain[2];
	mNonlocalTotalStrain[rTimeDerivative][3] = rNonlocalTotalStrain[3];
	mNonlocalTotalStrain[rTimeDerivative][4] = rNonlocalTotalStrain[4];
	mNonlocalTotalStrain[rTimeDerivative][5] = rNonlocalTotalStrain[5];
}

//! @brief returns the nonlocal total strain component of the node
//! @return strain component (rTimeDerivative=0)
template <NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNonlocalTotalStrain(short rIndex)const
{
	assert(rIndex<TNumNonlocalTotalStrain);
	return mNonlocalTotalStrain[0][rIndex];
}

//! @brief gives the global DOF of a Damage component
//! @param rComponent component
//! @return global DOF
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDofNonlocalTotalStrain(int rComponent)const
{
	assert(rComponent<TNumNonlocalTotalStrain);
	return mDofNonlocalTotalStrain[rComponent];

}


//! @brief returns the number of nonlocal eq. strain of the node
//! @return number of temperatures
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNumNonlocalEqStrain()const
{
    return TNumNonlocalEqStrain;
}

//! @brief returns the nonlocal eq. strain of the node
//! @return nonlocal eq. strain
template <NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNonlocalEqStrain()const
{
    assert(TNumNonlocalEqStrain==1);
    return mNonlocalEqStrain[0][0];
}

//! @brief returns the nonlocal eq. strain of the node
//! @return nonlocal eq. strain
template <NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNonlocalEqStrain(int rTimeDerivative)const
{
    assert(TNumNonlocalEqStrain==1);
    assert(rTimeDerivative>=0);
    assert(rTimeDerivative<=TNumTimeDerivatives);
    return mNonlocalEqStrain[rTimeDerivative][0];
}

//! @brief set the nonlocal eq. strain of the node
//! @param rNonlocalEqStrain  given nonlocal eq. strain
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetNonlocalEqStrain(double rNonlocalEqStrain)
{
    assert(TNumNonlocalEqStrain==1);
    mNonlocalEqStrain[0][0] = rNonlocalEqStrain;
}

//! @brief set the nonlocal eq. strain of the node
//! @param rNonlocalEqStrain  given nonlocal eq. strain
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetNonlocalEqStrain(int rTimeDerivative, double rNonlocalEqStrain)
{
    assert(TNumNonlocalEqStrain==1);
    assert(rTimeDerivative>=0);
    assert(rTimeDerivative<=TNumNonlocalEqStrain);
    mNonlocalEqStrain[rTimeDerivative][0] = rNonlocalEqStrain;
}

//! @brief gives the global DOF of a nonlocal eq. strain component
//! @param rComponent component
//! @return global DOF
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDofNonlocalEqStrain()const
{
    assert(TNumNonlocalEqStrain==1);
    return mDofNonlocalEqStrain[0];

}

// Moisture Transport --- Begin

// WaterPhaseFraction, int TNumRelativeHumidity

//! @brief returns the number of water phase fraction components of the node
//! @return number of water phase fraction components
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNumWaterPhaseFraction()const
{
    return TNumWaterPhaseFraction;
}

//! @brief returns the water phase fraction of the node
//! @return water phase fraction
template <NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetWaterPhaseFraction()const
{
    assert(TNumWaterPhaseFraction==1);
    return mWaterPhaseFraction[0][0];
}

//! @brief returns the water phase fraction of the node
//! @param rTimeDerivative time derivative
//! @return water phase fraction
template <NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetWaterPhaseFraction(int rTimeDerivative)const
{
    assert(TNumWaterPhaseFraction==1);
    assert(rTimeDerivative>=0 && rTimeDerivative<=TNumTimeDerivatives);
    return mWaterPhaseFraction[rTimeDerivative][0];
}

//! @brief set the water phase fraction of the node
//! @param rTemperature  given temperature
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetWaterPhaseFraction(double rWaterPhaseFraction)
{
    assert(TNumWaterPhaseFraction==1);
    mWaterPhaseFraction[0][0] = rWaterPhaseFraction;
}

//! @brief set the water phase fraction of the node
//! @param rTimeDerivative time derivative
//! @param rTemperature  given temperature
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetWaterPhaseFraction(int rTimeDerivative, double rWaterPhaseFraction)
{
    assert(TNumWaterPhaseFraction==1);
    assert(rTimeDerivative>=0 && rTimeDerivative<=TNumTimeDerivatives);
    mWaterPhaseFraction[0][0] = rWaterPhaseFraction;
}

//! @brief gives the global DOF of a water phase fraction component
//! @param rComponent component
//! @return global DOF
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDofWaterPhaseFraction()const
{
    assert(TNumWaterPhaseFraction==1);
    return mDofWaterPhaseFraction[0];
}

// --- relative Feuchte

//! @brief returns the number of relative humidity components of the node
//! @return number of relative humidity components
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNumRelativeHumidity()const
{
    return TNumRelativeHumidity;
}

//! @brief returns the relative humidity of the node
//! @return relative humidity
template <NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetRelativeHumidity()const
{
    assert(TNumRelativeHumidity==1);
    return mRelativeHumidity[0][0];
}

//! @brief returns the relative humidity of the node
//! @param rTimeDerivative time derivative
//! @return relative humidity
template <NODE_DOF_TEMPLATE_PARAMETERS>
double NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetRelativeHumidity(int rTimeDerivative)const
{
    assert(TNumRelativeHumidity==1);
    assert(rTimeDerivative>=0 && rTimeDerivative<=TNumTimeDerivatives);
    return mRelativeHumidity[rTimeDerivative][0];
}

//! @brief set the relative humidity of the node
//! @param rTemperature  given relative humidity
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetRelativeHumidity(double rRelativeHumidity)
{
    assert(TNumRelativeHumidity==1);
    mRelativeHumidity[0][0] = rRelativeHumidity;
}

//! @brief set the relative humidity of the node
//! @param rTimeDerivative time derivative
//! @param rTemperature  given relative humidity
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
SetRelativeHumidity(int rTimeDerivative, double rRelativeHumidity)
{
    assert(TNumRelativeHumidity==1);
    assert(rTimeDerivative>=0 && rTimeDerivative<=TNumTimeDerivatives);
    mRelativeHumidity[0][0] = rRelativeHumidity;
}

//! @brief gives the global DOF of a relative humidity component
//! @param rComponent component
//! @return global DOF
template <NODE_DOF_TEMPLATE_PARAMETERS>
int NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetDofRelativeHumidity()const
{
    assert(TNumRelativeHumidity==1);
    return mDofRelativeHumidity[0];
}



// Moisture Transport --- End





//! @brief returns the type of node as a string (all the data stored at the node)
//! @return string
template <NODE_DOF_TEMPLATE_PARAMETERS>
std::string NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
GetNodeTypeStr()const
{
	throw MechanicsException("");
    return std::string("[NuTo::NodeDof::GetNodeTypeStr]to be done");
}

/*
//! @brief returns the type of node as an enum (all the data stored at the node)
//! @return enum
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures >
NuTo::Node::eNodeType NuTo::NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>::
GetNodeType()const
{
	return Nodes::
    NodeCoordinates1D,
    NodeCoordinates2D,
    NodeCoordinates3D,
    NodeCoordinatesTimeDerivative_0_Displacements1D,
    NodeCoordinatesTimeDerivative_0_Displacements2D,
    NodeCoordinatesTimeDerivative_0_Displacements3D,
    NodeCoordinatesTimeDerivative_0_DisplacementsMultiscale1D,
    NodeCoordinatesTimeDerivative_0_DisplacementsMultiscale2D,
    NodeCoordinatesTimeDerivative_0_DisplacementsMultiscale3D,
    NodeCoordinatesTimeDerivative_0_DisplacementsNonlocalData2D,
    NodeCoordinatesTimeDerivative_0_DisplacementsNonlocalData3D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRadius2D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRadius3D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRotations2D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRotations3D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRotationsRadius2D,
    NodeCoordinatesTimeDerivative_0_DisplacementsRotationsRadius3D,
    NodeCoordinatesTimeDerivative_2_DisplacementsRotations2D,
    NodeCoordinatesTimeDerivative_2_DisplacementsRotationsRadius2D,
    NodeCoordinatesTimeDerivative_2_Displacements1D,
    NodeCoordinatesTimeDerivative_2_Displacements2D,
    NodeCoordinatesTimeDerivative_2_Displacements3D,
    NodeCoordinatesTimeDerivative_0_Temperature1D,
    NodeCoordinatesTimeDerivative_0_Temperature2D,
    NodeCoordinatesTimeDerivative_0_Temperature3D,

}
*/

#ifdef ENABLE_VISUALIZE
template <NODE_DOF_TEMPLATE_PARAMETERS>
void NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const
{

}
#endif // ENABLE_VISUALIZE

//! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
template <NODE_DOF_TEMPLATE_PARAMETERS>
NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>* NuTo::NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>::
Clone()const
{
    return new NodeDof<NODE_DOF_TEMPLATE_INITIALIZATION>(*this);
}


/*
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,1,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,2,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,3,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,1,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,2,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,3,0,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,2,1,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,3,3,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,2,1,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,3,3,0,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,0,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1,0,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,1,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,2,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,3,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,1,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,2,0,1,0>)))
BOOST_CLASS_EXPORT_IMPLEMENT(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,3,0,1,0>)))
#endif // ENABLE_SERIALIZATION
*/
#endif //NODEDOF_H

