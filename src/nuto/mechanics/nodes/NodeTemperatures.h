// $Id: $
#ifndef NodeTemperature_H
#define NodeTemperature_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeTemperaturesBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for nodes having rotational degrees of freedom
template<int NUMTEMPERATURES>
class NodeTemperatures : public NodeTemperaturesBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeTemperatures() : NodeTemperaturesBase ()
    {
        for (int count=0; count<NUMTEMPERATURES; count++)
            mTemperatures[count]=0;
    }

    //! @brief constructor
    NodeTemperatures (const double rTemperatures[NUMTEMPERATURES])  : NodeTemperaturesBase ()
    {
        memcpy(mTemperatures,rTemperatures,NUMTEMPERATURES*sizeof(double));
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
        & BOOST_SERIALIZATION_NVP(mTemperatures);
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the number of temperatures of the node
    //! @return number of temperatures
    int GetNumTemperatures()const
    {
        return NUMTEMPERATURES;
    }

    //! @brief gives the global DOF of a temperature component
    //! @param rComponent component
    //! @return global DOF
    int GetDofTemperature(int rComponent)const
    {
        assert(rComponent>=0 && rComponent<NUMTEMPERATURES);
        return this->mDOF[rComponent];
    }

    //! @brief set the rotations
    //! @param rRotations  given rotations
    void SetTemperatures(const double rTemperatures[NUMTEMPERATURES])
    {
        memcpy(mTemperatures,rTemperatures,NUMTEMPERATURES*sizeof(double));
    }

    //! @brief writes the temperature of a node to the prescribed pointer
    //! @param rTemperatur temperature
    void GetTemperatur(double rTemperatures[NUMTEMPERATURES])const
    {
        memcpy(rTemperatures,mTemperatures,NUMTEMPERATURES*sizeof(double));
    }

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        for (int count=0; count<NUMTEMPERATURES; count++)
        {
            mDOF[count]=rDOF++;
        }
    }


    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
    {
        assert(rActiveDofValues.GetNumColumns() == 1);
        assert(rDependentDofValues.GetNumColumns() == 1);
        for (int count=0; count<NUMTEMPERATURES; count++)
        {
            int dof = this->mDOF[count];
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
            this->mTemperatures[count] = value;
        }
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        assert(rActiveDofValues.GetNumColumns() == 1);
        assert(rDependentDofValues.GetNumColumns() == 1);
        for (int count=0; count<NUMTEMPERATURES; count++)
        {
            int dof = this->mDOF[count];
            double value = this->mTemperatures[count];
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
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
    {
        for (int count=0; count<NUMTEMPERATURES; count++)
        {
            mDOF[count]=rMappingInitialToNewOrdering[mDOF[count]];
        }
    }

protected:
    double mTemperatures[NUMTEMPERATURES];
    int mDOF[NUMTEMPERATURES];
};

}

#endif //NodeTemperature_H
