#ifndef NODEBASE_H
#define NODEBASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif  // ENABLE_SERIALIZATION
#include <vector>

#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
template <class T> class FullMatrix;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all nodes
class NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    NodeBase();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF)=0;

    //! @brief write dof values to the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues) = 0;

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const = 0;

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering) = 0;

    //! @brief returns the number of coordinates of the node
    //! @return number of coordinates
    virtual int GetNumCoordinates()const;

    //! @brief returns the coordinates of the node
    //! @return coordinates
    virtual void GetCoordinates1D(double rCoordinates[1])const;

    //! @brief set the coordinates
    //! @param rCoordinates  given coordinates
    virtual void SetCoordinates1D(const double rCoordinates[1]);

    //! @brief returns the coordinates of the node
    //! @return coordinates
    virtual void GetCoordinates2D(double rCoordinates[2])const;

    //! @brief set the coordinates
    //! @param rCoordinates  given coordinates
    virtual void SetCoordinates2D(const double rCoordinates[2]);

    //! @brief returns the coordinates of the node
    //! @return coordinates
    virtual void GetCoordinates3D(double rCoordinates[3])const;

    //! @brief set the coordinates
    //! @param rCoordinates  given coordinates
    virtual void SetCoordinates3D(const double rCoordinates[3]);

    //! @brief returns the number of coordinates of the node
    //! @return coordinates
    virtual double GetCoordinate(short rIndex)const;

    //! @brief returns the number of displacements of the node
    //! @return number of displacements
    virtual int GetNumDisplacements()const;

    //! @brief gives the global DOF of a displacement component
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofDisplacement(int rComponent)const;

    //! @brief returns the displacements of the node
    //! @return displacement
    virtual void GetDisplacements1D(double rCoordinates[1])const;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    virtual void SetDisplacements1D(const double rDisplacements[1]);

    //! @brief returns the displacements of the node
    //! @return displacement
    virtual void GetDisplacements2D(double rCoordinates[2])const;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    virtual void SetDisplacements2D(const double rDisplacements[2]);

    //! @brief returns the displacements of the node
    //! @return displacement
    virtual void GetDisplacements3D(double rCoordinates[3])const;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    virtual void SetDisplacements3D(const double rDisplacements[3]);

    //! @brief returns the displacements of the node
    //! @return displacement
    virtual double GetDisplacement(short rIndex)const;

    //! @brief returns the number of Rotations of the node
    //! @return number of Rotations
    virtual int GetNumRotations()const;

    //! @brief gives the global DOF of a Rotation component
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofRotation(int rComponent)const;

    //! @brief returns the Rotations of the node
    //! @return Rotation
    virtual void GetRotations2D(double rCoordinates[2])const;

    //! @brief set the Rotations
    //! @param rRotations  given Rotations
    virtual void SetRotations2D(const double rRotations[2]);

    //! @brief returns the Rotations of the node
    //! @return Rotation
    virtual void GetRotations3D(double rCoordinates[3])const;

    //! @brief set the Rotations
    //! @param rRotations  given Rotations
    virtual void SetRotations3D(const double rRotations[3]);

    //! @brief returns the Rotations of the node
    //! @return Rotation
    virtual double GetRotation(short rIndex)const;

    //! @brief returns the number of temperatures of the node
    //! @return number of temperatures
    virtual int GetNumTemperatures()const;

    //! @brief gives the global DOF of a temperature component
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofTemperature(int rComponent)const;

    //! @brief gives the grid number of the node, this is only implemented for the Grid node, since it stores the grid number
    //! @return NodeGridNum
    virtual int GetNodeGridNum()const;

    //! @brief returns the number of temperatures of the node
    //! @return number of temperatures
    virtual std::string GetNodeTypeStr()const=0;

protected:
    //the base class of the nodes must not contain any data

};
}//namespace NuTo
#endif //NODEBASE_H

