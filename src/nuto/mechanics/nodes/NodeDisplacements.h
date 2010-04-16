#ifndef NODE_DISPLACEMENTS_H
#define NODE_DISPLACEMENTS_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeDisplacementsBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for nodes having displacement degrees of freedom
template<int NUMDISPLACEMENTS>
class NodeDisplacements : public NodeDisplacementsBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeDisplacements() : NodeDisplacementsBase ()
    {
        for (int count=0; count<NUMDISPLACEMENTS; count++)
            mDisplacements[count]=0;
    }

    //! @brief constructor
    NodeDisplacements (const double rDisplacements[NUMDISPLACEMENTS])  : NodeDisplacementsBase ()
    {
        memcpy(mDisplacements,rDisplacements,NUMDISPLACEMENTS*sizeof(double));
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
        & BOOST_SERIALIZATION_NVP(mDisplacements);
    }
#endif // ENABLE_SERIALIZATION

    //! @brief returns the number of displacements of the node
    //! @return number of displacements
    int GetNumDisplacements()const
    {
        return NUMDISPLACEMENTS;
    }

    //! @brief gives the global DOF of a displacement component
    //! @param rComponent component
    //! @return global DOF
    int GetDofDisplacement(int rComponent)const
    {
        assert(rComponent>=0 && rComponent<NUMDISPLACEMENTS);
        return mDOF[rComponent];
    }

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    void SetDisplacements(const double rDisplacements[NUMDISPLACEMENTS])
    {
        memcpy(mDisplacements,rDisplacements,NUMDISPLACEMENTS*sizeof(double));
    }

    //! @brief writes the displacements of a node to the prescribed pointer
    //! @param rDisplacements displacements
    void GetDisplacements(double rDisplacements[NUMDISPLACEMENTS])const
    {
        memcpy(rDisplacements,mDisplacements,NUMDISPLACEMENTS*sizeof(double));
    }

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)
    {
        for (int count=0; count<NUMDISPLACEMENTS; count++)
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
        for (int count=0; count<NUMDISPLACEMENTS; count++)
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
            this->mDisplacements[count] = value;
        }
    }

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
    {
        assert(rActiveDofValues.GetNumColumns() == 1);
        assert(rDependentDofValues.GetNumColumns() == 1);
        for (int count=0; count<NUMDISPLACEMENTS; count++)
        {
            int dof = this->mDOF[count];
            double value = this->mDisplacements[count];
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
        for (int count=0; count<NUMDISPLACEMENTS; count++)
        {
            mDOF[count]=rMappingInitialToNewOrdering[mDOF[count]];
        }
    }

protected:
    double mDisplacements[NUMDISPLACEMENTS];
    int mDOF[NUMDISPLACEMENTS];

};
}//namespace NuTo
#endif //NODE_DISPLACEMENTS_H
