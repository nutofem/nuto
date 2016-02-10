// $Id$
#ifndef NODEBASE_H
#define NODEBASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/array.hpp>
#include "nuto/math/EigenBoostSerialization.h"
#else
#include <boost/array.hpp>
#endif  // ENABLE_SERIALIZATION

#include <array>
#include <assert.h>

#include "nuto/math/FullMatrix_Def.h"
#include "nuto/math/FullVector_Def.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeComponent.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

namespace NuTo
{
class NodeDisplacementsMultiscale2D;
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

    //! @brief destructor
    virtual ~NodeBase(){};

    //! @brief assignment operator
    NodeBase& operator=(NodeBase const& rOther) = default;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief sets the global dofs
    //! @param rDOF current maximum DOF, this variable is increased within the routine
    virtual void SetGlobalDofs(int& rDOF);

    //! @brief write dof values to the node (based on global dof number)
    //! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(int rTimeDerivative, const FullVector<double,Eigen::Dynamic>& rActiveDofValues, const FullVector<double,Eigen::Dynamic>& rDependentDofValues);

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(int rTimeDerivative, FullVector<double,Eigen::Dynamic>& rActiveDofValues, FullVector<double,Eigen::Dynamic>& rDependentDofValues) const;

    //! @brief extract all dof numbers from the node (based on global dof number)
    //virtual int* GetGlobalDofs();

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering);

    //! @brief returns the number of time derivatives stored at the node
    //! @return number of derivatives
    virtual int GetNumTimeDerivatives()const;

    //*************************************************
    //************     COORDINATES      ***************
    //*************************************************
    virtual int GetNumCoordinates() const;
    virtual double GetCoordinate(short rComponent) const;

    virtual const Eigen::Matrix<double, 1, 1>& GetCoordinates1D() const;
    virtual const Eigen::Matrix<double, 2, 1>& GetCoordinates2D() const;
    virtual const Eigen::Matrix<double, 3, 1>& GetCoordinates3D() const;
    virtual const Eigen::Matrix<double, Eigen::Dynamic, 1> GetCoordinates() const;

    virtual void SetCoordinates1D(const Eigen::Matrix<double, 1, 1>& rCoordinates);
    virtual void SetCoordinates2D(const Eigen::Matrix<double, 2, 1>& rCoordinates);
    virtual void SetCoordinates3D(const Eigen::Matrix<double, 3, 1>& rCoordinates);
    virtual void SetCoordinates  (const Eigen::Matrix<double, Eigen::Dynamic, 1>& rCoordinates);



    //*************************************************
    //************    DISPLACEMENTS     ***************
    //*************************************************
    virtual int GetNumDisplacements() const;
    virtual int GetDofDisplacement(int rComponent) const;
    virtual double GetDisplacement(short rIndex) const;

    const Eigen::Matrix<double, 1, 1>& GetDisplacements1D() const;
    const Eigen::Matrix<double, 2, 1>& GetDisplacements2D() const;
    const Eigen::Matrix<double, 3, 1>& GetDisplacements3D() const;
    const Eigen::Matrix<double, Eigen::Dynamic, 1> GetDisplacements() const;

    virtual const Eigen::Matrix<double, 1, 1>& GetDisplacements1D(int rTimeDerivative) const;
    virtual const Eigen::Matrix<double, 2, 1>& GetDisplacements2D(int rTimeDerivative) const;
    virtual const Eigen::Matrix<double, 3, 1>& GetDisplacements3D(int rTimeDerivative) const;
    virtual const Eigen::Matrix<double, Eigen::Dynamic, 1> GetDisplacements(int rTimeDerivative) const;

    void SetDisplacements1D(const Eigen::Matrix<double, 1, 1>& rDisplacements);
    void SetDisplacements2D(const Eigen::Matrix<double, 2, 1>& rDisplacements);
    void SetDisplacements3D(const Eigen::Matrix<double, 3, 1>& rDisplacements);
    void SetDisplacements  (const Eigen::Matrix<double, Eigen::Dynamic, 1>& rDisplacements);

    virtual void SetDisplacements1D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rDisplacements);
    virtual void SetDisplacements2D(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rDisplacements);
    virtual void SetDisplacements3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rDisplacements);
    virtual void SetDisplacements  (int rTimeDerivative, const Eigen::Matrix<double, Eigen::Dynamic, 1>& rDisplacements);

    //*************************************************
    //************       ROTATIONS      ***************
    //*************************************************

    virtual int GetNumRotations()const;
    virtual int GetDofRotation(int rComponent)const;
    virtual double GetRotation(short rIndex)const;

    const Eigen::Matrix<double, 1, 1>& GetRotations2D() const;
    const Eigen::Matrix<double, 3, 1>& GetRotations3D() const;
    const Eigen::Matrix<double, Eigen::Dynamic, 1> GetRotations() const;

    virtual const Eigen::Matrix<double, 1, 1>& GetRotations2D(int rTimeDerivative) const;
    virtual const Eigen::Matrix<double, 3, 1>& GetRotations3D(int rTimeDerivative) const;
    virtual const Eigen::Matrix<double, Eigen::Dynamic, 1> GetRotations(int rTimeDerivative) const;

    void SetRotations2D(const Eigen::Matrix<double, 1, 1>& rRotations);
    void SetRotations3D(const Eigen::Matrix<double, 3, 1>& rRotations);
    void SetRotations  (const Eigen::Matrix<double, Eigen::Dynamic, 1>& rRotations);


    virtual void SetRotations2D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rRotations);
    virtual void SetRotations3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rRotations);
    virtual void SetRotations  (int rTimeDerivative, const Eigen::Matrix<double, Eigen::Dynamic, 1>& rRotations);


    //*************************************************
    //************      TEMPERATURE     ***************
    //*************************************************

    virtual int GetNumTemperatures() const;
    virtual int GetDofTemperature() const;

    double GetTemperature() const;
    virtual double GetTemperature(int rTimeDerivative) const;

    void SetTemperature(double rTemperature);
    virtual void SetTemperature(int rTimeDerivative, double rTemperature);

    //*************************************************
    //********  NONLOCAL EQ PLASTIC STRAIN  ***********
    //*************************************************

    virtual int GetNumNonlocalEqPlasticStrain() const;
    virtual int GetDofNonlocalEqPlasticStrain(int rComponent) const;

    const Eigen::Matrix<double, 2, 1>& GetNonlocalEqPlasticStrain() const;
    virtual const Eigen::Matrix<double, 2, 1>& GetNonlocalEqPlasticStrain(int rTimeDerivative) const;

    void SetNonlocalEqPlasticStrain(const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain);
    virtual void SetNonlocalEqPlasticStrain(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain);

    //*************************************************
    //********    NONLOCAL TOTAL STRAIN     ***********
    //*************************************************

    virtual int GetNumNonlocalTotalStrain()const;
    virtual int GetDofNonlocalTotalStrain(int rComponent)const;
    virtual double GetNonlocalTotalStrain(short rIndex)const;

    const Eigen::Matrix<double, 1, 1>& GetNonlocalTotalStrain1D() const;
    const Eigen::Matrix<double, 3, 1>& GetNonlocalTotalStrain2D() const;
    const Eigen::Matrix<double, 6, 1>& GetNonlocalTotalStrain3D() const;
    const Eigen::Matrix<double, Eigen::Dynamic, 1> GetNonlocalTotalStrains() const;

    virtual const Eigen::Matrix<double, 1, 1>& GetNonlocalTotalStrain1D(int rTimeDerivative) const;
    virtual const Eigen::Matrix<double, 3, 1>& GetNonlocalTotalStrain2D(int rTimeDerivative) const;
    virtual const Eigen::Matrix<double, 6, 1>& GetNonlocalTotalStrain3D(int rTimeDerivative) const;
    virtual const Eigen::Matrix<double, Eigen::Dynamic, 1> GetNonlocalTotalStrains(int rTimeDerivative) const;

    void SetNonlocalTotalStrain1D(const Eigen::Matrix<double, 1, 1>& rNonlocalTotalStrain);
    void SetNonlocalTotalStrain2D(const Eigen::Matrix<double, 3, 1>& rNonlocalTotalStrain);
    void SetNonlocalTotalStrain3D(const Eigen::Matrix<double, 6, 1>& rNonlocalTotalStrain);
    void SetNonlocalTotalStrain  (const Eigen::Matrix<double, Eigen::Dynamic, 1>& rNonlocalTotalStrain);

    virtual void SetNonlocalTotalStrain1D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rNonlocalTotalStrain);
    virtual void SetNonlocalTotalStrain2D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rNonlocalTotalStrain);
    virtual void SetNonlocalTotalStrain3D(int rTimeDerivative, const Eigen::Matrix<double, 6, 1>& rNonlocalTotalStrain);
    virtual void SetNonlocalTotalStrain  (int rTimeDerivative, const Eigen::Matrix<double, Eigen::Dynamic, 1>& rNonlocalTotalStrain);

    //*************************************************
    //*******   NONLOCAL EQUIVALENT STRAIN   **********
    //*************************************************

    virtual int GetNumNonlocalEqStrain() const;
    virtual int GetDofNonlocalEqStrain() const;

    double GetNonlocalEqStrain() const;
    virtual double GetNonlocalEqStrain(int rTimeDerivative) const;

    void SetNonlocalEqStrain(double rNonlocalEqStrain);
    virtual void SetNonlocalEqStrain(int rTimeDerivative, double rNonlocalEqStrain);

    //**************************************************
    //*******      WATER VOLUME FRACTION      **********
    //**************************************************

    virtual int GetNumWaterVolumeFraction() const;
    virtual int GetDofWaterVolumeFraction() const;

    double GetWaterVolumeFraction() const;
    virtual double GetWaterVolumeFraction(int rTimeDerivative) const;

    void SetWaterVolumeFraction(double rWaterVolumeFraction);
    virtual void SetWaterVolumeFraction(int rTimeDerivative, double rWaterVolumeFraction);

    //*************************************************
    //*******       RELATIVE HUMIDITY        **********
    //*************************************************

    virtual int GetNumRelativeHumidity() const;
    virtual int GetDofRelativeHumidity() const;

    double GetRelativeHumidity() const;
    virtual double GetRelativeHumidity(int rTimeDerivative) const;

    void SetRelativeHumidity(double rRelativeHumidity);
    virtual void SetRelativeHumidity(int rTimeDerivative, double rRelativeHumidity);


    //! @brief returns the type of node as a string (all the data stored at the node)
    //! @return string
    virtual std::string GetNodeTypeStr()const=0;

#ifdef ENABLE_VISUALIZE
    virtual void Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList) const;

#endif // ENABLE_VISUALIZE

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    virtual NodeBase* Clone()const=0;

protected:
    //the base class of the nodes must not contain any data

};

class less_XCoordinate2D : public std::binary_function<NodeBase*, NodeBase* , bool>
{
public:

    //! @brief sorts the nodes in increasing x-direction
    less_XCoordinate2D()
    {
    }

    bool operator()(NodeBase* nodePtr1, NodeBase* nodePtr2)
    {
        return nodePtr1->GetCoordinate(0) < nodePtr2->GetCoordinate(0);
    }
};

class greater_XCoordinate2D : public std::binary_function<NodeBase*, NodeBase* , bool>
{
public:
    greater_XCoordinate2D()
    {
    }

    bool operator()(NodeBase* nodePtr1, NodeBase* nodePtr2)
    {
        return nodePtr1->GetCoordinate(0) > nodePtr2->GetCoordinate(0);
    }
};

class less_YCoordinate2D : public std::binary_function<NodeBase*, NodeBase* , bool>
{
public:
    less_YCoordinate2D()
    {
    }

    bool operator()(NodeBase* nodePtr1, NodeBase* nodePtr2)
    {
        return nodePtr1->GetCoordinate(1) < nodePtr2->GetCoordinate(1);
    }
};

class greater_YCoordinate2D : public std::binary_function<NodeBase*, NodeBase* , bool>
{
public:
    greater_YCoordinate2D()
    {
    }

    bool operator()(NodeBase* nodePtr1, NodeBase* nodePtr2)
    {
        return nodePtr1->GetCoordinate(1) > nodePtr2->GetCoordinate(1);
    }
};

}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
//namespace boost
//{
//    //! @brief tell boost how to serialize an Eigen::Matrix
//    template<class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
//    inline void serialize(
//        Archive & ar,
//        Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & t,
//        const unsigned int file_version
//    )
//    {
//        ar & boost::serialization::make_array(t.data(), t.size());
//    }
//}
BOOST_CLASS_EXPORT_KEY(NuTo::NodeBase)
#endif // ENABLE_SERIALIZATION
#endif //NODEBASE_H

