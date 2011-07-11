// $Id$
#ifndef NODE_GRID_3D_H
#define NODE_GRID_3D_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for 3D reference nodes
class NodeGrid3D : public virtual NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

private:
    //! @brief constructor
     NodeGrid3D(){}

public:
   //! @brief constructor
    NodeGrid3D(int rNodeGridNum);

    //! @brief destructor
    ~ NodeGrid3D();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the number of coordinates of the node
    //! @return number of coordinates
    virtual int GetNumCoordinates()const;

    //! @brief set the grid node number
    //! @param rNodeGridNum  given grid node number
    virtual void SetNodeGridNum(int rNodeGridNum);

    //! @brief get the grid node number
    //! @return rNodeGridNum  given grid node number
    int GetNodeGridNum()const;

    //! @brief Get the Number of belonging elements
    //! @return Number of belonging elements
    int GetNumElems();

    //!@brief Set Number of belonging elements
    //! @param Number of belonging elements
    void SetNumElems(int rNumElems);


    //! @brief Get Ids of the elements
    //! @return int * Ids of the elements
    int* GetElementIds();

    //! @brief Set Ids of the elements
    //! @param int * Ids of the elements
    void SetElementIds(int * rElementIds);

    //! @brief Get Ids of the nodes
    //! @return int * Ids of the nodes
    int* GetNodeIds();

    //! @brief Set Ids of the nodes
    //! @param int * Ids of the nodes
    void SetNodeIds(int * rNodeIds);

    FullMatrix<double>* GetPartCoefficientMatrix0(int node);
    double* GetPartCoefficient0(int node);

    void SetPartCoefficientMatrix0(int node, NuTo::FullMatrix<double>& rCoefficientMatrix0);
    void SetPartCoefficient0(int node, NuTo::FullMatrix<double>& rCoefficientMatrix0);

    //! @brief writes the coordinates of a node to the prescribed pointer
    //! @param rCoordinates coordinates
    virtual void GetCoordinates3D(double rCoordinates[3])const;

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

protected:
    //! @brief stores the grid number of the node
    //! @param grid number of the node

    int mNodeGridNum;
    int mNumElems;
    int* mElementIds;
    int* mNodeIds;
    std::vector<FullMatrix<double>*> mCoefficientMatrix0;
    double* mCoefficient0;
};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NodeGrid3D)
#endif // ENABLE_SERIALIZATION

#endif //NODE_GRID_3D_H
