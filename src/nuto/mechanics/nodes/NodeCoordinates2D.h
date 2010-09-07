#ifndef NODE_COORDINATES_2D_H
#define NODE_COORDINATES_2D_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeBase.h"
namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for 2D reference nodes
class NodeCoordinates2D : public virtual NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief constructor
    NodeCoordinates2D();

    //! @brief constructor
    NodeCoordinates2D(const double rCoordinates[2]);

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

    //! @brief set the coordinates
    //! @param rCoordinates  given coordinates
    virtual void SetCoordinates2D(const double rCoordinates[2]);

    //! @brief writes the coordinates of a node to the prescribed pointer
    //! @param rCoordinates coordinates
    virtual void GetCoordinates2D(double rCoordinates[2])const;

    //! @brief returns the coordinate of a given direction of the node
    //! @param rIndex index of the direction
    //! @return coordinate
    virtual double GetCoordinate(short rIndex)const;

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

protected:
    double mCoordinates[2];

};

class less_XCoordinate2D : public std::binary_function<NodeCoordinates2D*, NodeCoordinates2D* , bool>
{
public:

    //! @brief sorts the nodes in increasing x-direction with a shift, that is subtracted
    //! and if its less than zero, it is added right to the maximum coordinates
    //! this actually corresponds to sorting the nodes according to their x value,
    //! and then cutting at a certain position (shift) and move all the previous nodes at the end
    //! while keeping their order
    //! if shift is set to zero (and maxcoordinate any value), the standard ordering is obtained
    less_XCoordinate2D(double rShift, double rMaxCoordinate)
    {
        mShift  = rShift;
        //use factor two. in order to make sure that the shiftes nodes have a higher node value than the nonshifted ones
        m2MaxCoordinate = 2.*rMaxCoordinate;
    }

    bool operator()(NodeCoordinates2D* nodePtr1, NodeBase* nodePtr2)
    {
        double coord1[2], coord2[2];
        nodePtr1->GetCoordinates2D(coord1);
        nodePtr2->GetCoordinates2D(coord2);
        coord1[0]-=mShift;
        coord2[0]-=mShift;
        if(coord1[0]<0)
            coord1[0]+=m2MaxCoordinate;
        if(coord2[0]<0)
            coord2[0]+=m2MaxCoordinate;
        return coord1[0] < coord2[0];
    }
protected:
    double mShift, m2MaxCoordinate;
};

class greater_XCoordinate2D : public std::binary_function<NodeCoordinates2D*, NodeCoordinates2D* , bool>
{
public:
    greater_XCoordinate2D(double rShift, double rMaxCoordinate)
    {
        mShift  = rShift;
        //use factor two. in order to make sure that the shiftes nodes have a higher node value than the nonshifted ones
        m2MaxCoordinate = 2.*rMaxCoordinate;
    }

    bool operator()(NodeCoordinates2D* nodePtr1, NodeBase* nodePtr2)
    {
        double coord1[2], coord2[2];
        nodePtr1->GetCoordinates2D(coord1);
        nodePtr2->GetCoordinates2D(coord2);
        coord1[0]-=mShift;
        coord2[0]-=mShift;
        if(coord1[0]<0)
            coord1[0]+=m2MaxCoordinate;
        if(coord2[0]<0)
            coord2[0]+=m2MaxCoordinate;
        return coord1[0] > coord2[0];
    }
protected:
    double mShift, m2MaxCoordinate;
};

class less_YCoordinate2D : public std::binary_function<NodeCoordinates2D*, NodeCoordinates2D* , bool>
{
public:
    less_YCoordinate2D(double rShift, double rMaxCoordinate)
    {
        mShift  = rShift;
        //use factor two. in order to make sure that the shiftes nodes have a higher node value than the nonshifted ones
        m2MaxCoordinate = 2.*rMaxCoordinate;
    }

    bool operator()(NodeCoordinates2D* nodePtr1, NodeBase* nodePtr2)
    {
        double coord1[2], coord2[2];
        nodePtr1->GetCoordinates2D(coord1);
        nodePtr2->GetCoordinates2D(coord2);
        coord1[1]-=mShift;
        coord2[1]-=mShift;
        if(coord1[1]<0)
            coord1[1]+=m2MaxCoordinate;
        if(coord2[1]<0)
            coord2[1]+=m2MaxCoordinate;
        return coord1[1] < coord2[1];
    }
protected:
    double mShift, m2MaxCoordinate;
};

class greater_YCoordinate2D : public std::binary_function<NodeCoordinates2D*, NodeCoordinates2D* , bool>
{
public:
    greater_YCoordinate2D(double rShift, double rMaxCoordinate)
    {
        mShift  = rShift;
        //use factor two. in order to make sure that the shiftes nodes have a higher node value than the nonshifted ones
        m2MaxCoordinate = 2.*rMaxCoordinate;
    }

    bool operator()(NodeCoordinates2D* nodePtr1, NodeBase* nodePtr2)
    {
        double coord1[2], coord2[2];
        nodePtr1->GetCoordinates2D(coord1);
        nodePtr2->GetCoordinates2D(coord2);
        coord1[1]-=mShift;
        coord2[1]-=mShift;
        if(coord1[1]<0)
            coord1[1]+=m2MaxCoordinate;
        if(coord2[1]<0)
            coord2[1]+=m2MaxCoordinate;

        return coord1[1] > coord2[1];
    }
protected:
    double mShift, m2MaxCoordinate;
};

}//namespace NuTo
#endif //NODE_COORDINATES_2D_H
