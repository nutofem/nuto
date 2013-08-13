// $Id$
#ifndef NODEDOF_DEF_H
#define NODEDOF_DEF_H

#ifdef ENABLE_SERIALIZATION
// serialize
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/utility/identity_type.hpp>
#include <boost/serialization/array.hpp>
#else
#include <boost/array.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/math/FullMatrix_Def.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeComponentBase.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for all nodes
template <int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain>
class NodeDof: public virtual NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    NodeDof();


    //! @brief destructor
    ~NodeDof();

    //! @brief assignment operator
    void operator=(NodeDof const& rOther);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp("NodeDof",boost::serialization::base_object< NodeBase >(*this));

    	ar & BOOST_SERIALIZATION_NVP(mDisplacements);
   	    ar & BOOST_SERIALIZATION_NVP(mDofDisplacements);

    	ar & BOOST_SERIALIZATION_NVP(mRotations);
    	ar & BOOST_SERIALIZATION_NVP(mDofRotations);

    	ar & BOOST_SERIALIZATION_NVP(mTemperatures);
    	ar & BOOST_SERIALIZATION_NVP(mDofTemperatures);

    	ar & BOOST_SERIALIZATION_NVP(mNonlocalEqPlasticStrain);
    	ar & BOOST_SERIALIZATION_NVP(mDofNonlocalEqPlasticStrain);

    	ar & BOOST_SERIALIZATION_NVP(mNonlocalTotalStrain);
    	ar & BOOST_SERIALIZATION_NVP(mDofNonlocalTotalStrain);
    }
#endif // ENABLE_SERIALIZATION

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    void SetGlobalDofs(int& rDOF) override ;

    //! @brief write dof values to the node (based on global dof number)
    //! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
    //! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void SetGlobalDofValues(int rTimeDerivative, const FullVector<double,Eigen::Dynamic>& rActiveDofValues, const FullVector<double,Eigen::Dynamic>& rDependentDofValues) override;

    //! @brief extract dof values from the node (based on global dof number)
    //! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
    //! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void GetGlobalDofValues(int rTimeDerivative, FullVector<double,Eigen::Dynamic>& rActiveDofValues, FullVector<double,Eigen::Dynamic>& rDependentDofValues) const override;

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering) override;

    //! @brief returns the number of time derivatives stored at the node
    //! @return number of derivatives
    int GetNumTimeDerivatives()const override;

    //! @brief returns the number of displacements of the node
    //! @return number of displacements
    int GetNumDisplacements()const override;

    //! @brief gives the global DOF of a displacement component
    //! @param rComponent component
    //! @return global DOF
    int GetDofDisplacement(int rComponent)const override;

    //! @brief returns the displacements of the node
    //! @return displacement
    void GetDisplacements1D(double rCoordinates[1])const override;

    //! @brief returns the displacements of the node
    //! @param rTimeDerivative time derivative
    //! @return displacement
    void GetDisplacements1D(int rTimeDerivative, double rDisplacements[1])const override;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    void SetDisplacements1D(const double rDisplacements[1]) override;

    //! @brief set the displacements
    //! @param rTimeDerivative time derivative
    //! @param rDisplacements  given displacements
    void SetDisplacements1D(int rTimeDerivative, const double rDisplacements[1]) override;

    //! @brief returns the displacements of the node
    //! @return displacement
    void GetDisplacements2D(double rCoordinates[2])const override;

    //! @brief returns the displacements of the node
    //! @param rTimeDerivative time derivative
    //! @return displacement
    void GetDisplacements2D(int rTimeDerivative, double rDisplacements[2])const override;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    void SetDisplacements2D(const double rDisplacements[2]) override;

    //! @brief set the displacements
    //! @param rTimeDerivative time derivative
    //! @param rDisplacements  given displacements
    void SetDisplacements2D(int rTimeDerivative, const double rDisplacements[2]) override;

    //! @brief returns the displacements of the node
    //! @return displacement
    void GetDisplacements3D(double rCoordinates[3])const override;

    //! @brief returns the displacements of the node
    //! @param rTimeDerivative time derivative
    //! @return displacement
    void GetDisplacements3D(int rTimeDerivative, double rDisplacements[3])const override;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    void SetDisplacements3D(const double rDisplacements[3]) override;

    //! @brief set the displacements
    //! @param rTimeDerivative time derivative
    //! @param rDisplacements  given displacements
    void SetDisplacements3D(int rTimeDerivative, const double rDisplacements[3])override;

    //! @brief returns the displacements of the node
    //! @return displacement
    double GetDisplacement(short rIndex)const override;

    //! @brief returns the number of Rotations of the node
    //! @return number of Rotations
    int GetNumRotations()const override;

    //! @brief gives the global DOF of a Rotation component
    //! @param rComponent component
    //! @return global DOF
    int GetDofRotation(int rComponent)const override;

    //! @brief returns the Rotations of the node
    //! @return Rotation
    void GetRotations2D(double rRotations[1])const override;

    //! @brief returns the Rotations of the node
    //! @param rTimeDerivative time derivative
    //! @return Rotation
    void GetRotations2D(int rTimeDerivative, double rRotations[1])const override;

    //! @brief set the Rotations
    //! @param rRotations  given Rotations
    void SetRotations2D(const double rRotations[1]) override;

    //! @brief set the Rotations
    //! @param rTimeDerivative time derivative
    //! @param rRotations  given Rotations
    void SetRotations2D(int rTimeDerivative, const double rRotations[1]) override;

    //! @brief returns the Rotations of the node
    //! @return Rotation
    void GetRotations3D(double rRotations[3])const override;

    //! @brief returns the Rotations of the node
    //! @param rTimeDerivative time derivative
    //! @return Rotation
    void GetRotations3D(int rTimeDerivative, double rRotations[1])const override;

    //! @brief set the Rotations
    //! @param rRotations  given Rotations
    void SetRotations3D(const double rRotations[3]) override;

    //! @brief set the Rotations
    //! @param rTimeDerivative time derivative
    //! @param rRotations  given Rotations
    void SetRotations3D(int rTimeDerivative, const double rRotations[3]) override;

    //! @brief returns the Rotations of the node
    //! @return Rotation
    double GetRotation(short rIndex)const override;

    //! @brief returns the number of temperatures of the node
    //! @return number of temperatures
    int GetNumTemperatures()const override;

    //! @brief returns the temperature of the node
    //! @return temperature
    double GetTemperature()const override;

    //! @brief returns the temperature of the node
    //! @param rTimeDerivative time derivative
    //! @return temperature
    double GetTemperature(int rTimeDerivative)const override;

    //! @brief set the temperature of the node
    //! @param rTemperature  given temperature
    void SetTemperature(double rTemperature) override;

    //! @brief set the temperature of the node
    //! @param rTimeDerivative time derivative
    //! @param rTemperature  given temperature
    void SetTemperature(int rTimeDerivative, double rTemperature) override;

    //! @brief gives the global DOF of a temperature component
    //! @param rComponent component
    //! @return global DOF
    int GetDofTemperature()const override;

    //! @brief returns the number of Damage dofs of the node
    //! @return number of Damages
    int GetNumNonlocalEqPlasticStrain()const override;

    //! @brief returns the Damage of the node
    //! @return Damage
    void GetNonlocalEqPlasticStrain(double* rNonlocalEqPlasticStrain)const override;

    //! @brief returns the Damage of the node
    //! @param rTimeDerivative time derivative
    //! @return Damage
    void GetNonlocalEqPlasticStrain(int rTimeDerivative, double* rNonlocalEqPlasticStrain)const override;

    //! @brief set the Damage of the node
    //! @param rDamage  given Damage
    void SetNonlocalEqPlasticStrain(const double* rNonlocalEqPlasticStrain) override;

    //! @brief set the Damage of the node
    //! @param rTimeDerivative time derivative
    //! @param rDamage  given Damage
    void SetNonlocalEqPlasticStrain(int rTimeDerivative, const double* rNonlocalEqPlasticStrain) override;

    //! @brief gives the global DOF of a Damage component
    //! @param rComponent component
    //! @return global DOF
    int GetDofNonlocalEqPlasticStrain(int rComponent)const override;

    //! @brief returns the number of Damage dofs of the node
    //! @return number of Damages
    int GetNumNonlocalTotalStrain()const override;

    //! @brief returns the Damage of the node
    //! @return Damage
    void GetNonlocalTotalStrain1D(double* rNonlocalTotalStrain)const override;

    //! @brief returns the Damage of the node
    //! @return Damage
    void GetNonlocalTotalStrain2D(double* rNonlocalTotalStrain)const override;

    //! @brief returns the Damage of the node
    //! @return Damage
    void GetNonlocalTotalStrain3D(double* rNonlocalTotalStrain)const override;

    //! @brief returns the Damage of the node
    //! @param rTimeDerivative time derivative
    //! @return Damage
    void GetNonlocalTotalStrain1D(int rTimeDerivative, double* rNonlocalTotalStrain)const override;

    //! @brief returns the Damage of the node
    //! @param rTimeDerivative time derivative
    //! @return Damage
    void GetNonlocalTotalStrain2D(int rTimeDerivative, double* rNonlocalTotalStrain)const override;

    //! @brief returns the Damage of the node
    //! @param rTimeDerivative time derivative
    //! @return Damage
    void GetNonlocalTotalStrain3D(int rTimeDerivative, double* rNonlocalTotalStrain)const override;

    //! @brief set the Damage of the node
    //! @param rDamage  given Damage
    void SetNonlocalTotalStrain1D(const double* rNonlocalTotalStrain) override;

    //! @brief set the Damage of the node
    //! @param rDamage  given Damage
    void SetNonlocalTotalStrain2D(const double* rNonlocalTotalStrain) override;

    //! @brief set the Damage of the node
    //! @param rDamage  given Damage
    void SetNonlocalTotalStrain3D(const double* rNonlocalTotalStrain) override;

    //! @brief set the Damage of the node
    //! @param rTimeDerivative time derivative
    //! @param rDamage  given Damage
    void SetNonlocalTotalStrain1D(int rTimeDerivative, const double* rNonlocalTotalStrain) override;

    //! @brief set the Damage of the node
    //! @param rTimeDerivative time derivative
    //! @param rDamage  given Damage
    void SetNonlocalTotalStrain2D(int rTimeDerivative, const double* rNonlocalTotalStrain) override;

    //! @brief set the Damage of the node
    //! @param rTimeDerivative time derivative
    //! @param rDamage  given Damage
    void SetNonlocalTotalStrain3D(int rTimeDerivative, const double* rNonlocalTotalStrain) override;

    //! @brief returns the nonlocal total strain component of the node
    //! @return strain component (rTimeDerivative=0)
    double GetNonlocalTotalStrain(short rIndex)const;

    //! @brief gives the global DOF of a Damage component
    //! @param rComponent component
    //! @return global DOF
    int GetDofNonlocalTotalStrain(int rComponent)const override;

    //! @brief returns the type of node as a string (all the data stored at the node)
    //! @return string
    std::string GetNodeTypeStr()const override;

    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    //Node::eNodeType GetNodeType()const;

#ifdef ENABLE_VISUALIZE
    void Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const;
#endif // ENABLE_VISUALIZE

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    virtual NodeDof* Clone()const;

protected:
    boost::array<boost::array<double, TNumDisplacements>,TNumTimeDerivatives+1> mDisplacements;
    boost::array<int, TNumDisplacements> mDofDisplacements;

    boost::array<boost::array<double, TNumRotations>,TNumTimeDerivatives+1> mRotations;
    boost::array<int, TNumRotations> mDofRotations;

    boost::array<boost::array<double, TNumTemperatures>,TNumTimeDerivatives+1> mTemperatures;
    boost::array<int, TNumTemperatures> mDofTemperatures;

    boost::array<boost::array<double, TNumNonlocalEqPlasticStrain>,TNumTimeDerivatives+1> mNonlocalEqPlasticStrain;
    boost::array<int, TNumNonlocalEqPlasticStrain> mDofNonlocalEqPlasticStrain;

    boost::array<boost::array<double, TNumNonlocalTotalStrain>,TNumTimeDerivatives+1> mNonlocalTotalStrain;
    boost::array<int, TNumNonlocalTotalStrain> mDofNonlocalTotalStrain;

};

}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
/*
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,1,0,0,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,2,0,0,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,3,0,0,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,1,0,0,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,2,0,0,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,3,0,0,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,2,1,0,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,3,3,0,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,2,1,0,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,3,3,0,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,0,0,1,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<1,0,0,1,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,1,0,1,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,2,0,1,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<0,3,0,1,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,1,0,1,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,2,0,1,0>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::NodeDof<2,3,0,1,0>)))
*/
#endif // ENABLE_SERIALIZATION

#endif //NODEDOF_DEF_H

