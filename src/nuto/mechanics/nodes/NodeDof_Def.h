// $Id$
#ifndef NodeDof_DEF_H
#define NodeDof_DEF_H

#ifdef ENABLE_SERIALIZATION
// serialize
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/utility/identity_type.hpp>
#include <boost/serialization/array.hpp>
#else
#include <boost/array.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeCoordinates.h"
#include "nuto/mechanics/nodes/NodeDisplacements.h"
#include "nuto/mechanics/nodes/NodeRotations.h"
#include "nuto/mechanics/nodes/NodeNonlocalEqPlasticStrain.h"
#include "nuto/mechanics/nodes/NodeNonlocalEqStrain.h"
#include "nuto/mechanics/nodes/NodeRelativeHumidity.h"
#include "nuto/mechanics/nodes/NodeWaterVolumeFraction.h"

#include "nuto/mechanics/nodes/NodeEnum.h"

#include "nuto/math/FullMatrix_Def.h"
#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeComponentBase.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

#define NODE_DOF_TEMPLATE_PARAMETERS int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain, int TNumNonlocalEqStrain, int TNumWaterVolumeFraction, int TNumRelativeHumidity
#define NODE_DOF_TEMPLATE_INITIALIZATION TNumCoordinates, TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures, TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain, TNumWaterVolumeFraction, TNumRelativeHumidity

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for all nodes
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain, int TNumNonlocalEqStrain, int TNumWaterVolumeFraction, int TNumRelativeHumidity>
class NodeDof:
        public NodeCoordinates<TNumCoordinates>,
        public NodeDisplacements<TNumDisplacements, TNumTimeDerivatives>,
        public NodeRotations<TNumRotations, TNumTimeDerivatives>,
        public NodeNonlocalEqPlasticStrain<TNumNonlocalEqPlasticStrain, TNumTimeDerivatives>,
        public NodeNonlocalEqStrain<TNumNonlocalEqStrain, TNumTimeDerivatives>,
        public NodeRelativeHumidity<TNumRelativeHumidity,TNumTimeDerivatives>,
        public NodeWaterVolumeFraction<TNumWaterVolumeFraction,TNumTimeDerivatives>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    NodeDof();

    //! @brief destructor
    ~NodeDof() {}

    //! @brief assignment operator
    NodeDof& operator=(NodeDof const& rOther) = default;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeDof" << std::endl;
#endif

        typedef NodeCoordinates<TNumCoordinates> Coordinates;
        typedef NodeDisplacements<TNumDisplacements, TNumTimeDerivatives> Displacements;
        typedef NodeRotations<TNumRotations, TNumTimeDerivatives> Rotations;
        typedef NodeNonlocalEqPlasticStrain<TNumNonlocalEqPlasticStrain, TNumTimeDerivatives> EqPlasticStrain;
        typedef NodeNonlocalEqStrain<TNumNonlocalEqStrain, TNumTimeDerivatives> EqStrain;
        typedef NodeRelativeHumidity<TNumRelativeHumidity,TNumTimeDerivatives> RelativeHumidity;
        typedef NodeWaterVolumeFraction<TNumWaterVolumeFraction,TNumTimeDerivatives> WaterVolumeFraction;

        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Coordinates);
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Displacements);
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Rotations);
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(EqPlasticStrain);
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(EqStrain);
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RelativeHumidity);
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(WaterVolumeFraction);

/*
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
*/
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeDof" << std::endl;
#endif
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
    int GetNumTimeDerivatives() const override;


    //! @brief returns the type of node as a string (all the data stored at the node)
    //! @return string
    std::string GetNodeTypeStr() const override;


//    //*************************************************
//    //************      TEMPERATURE     ***************
//    //*************************************************
//
//    int GetNumTemperatures() const override;
//    int GetDofTemperature() const override;
//
//    double GetTemperature() const override;
//    double GetTemperature(int rTimeDerivative) const override;
//
//    void SetTemperature(double rTemperature) override;
//    void SetTemperature(int rTimeDerivative, double rTemperature) override;
//
//    //*************************************************
//    //********  NONLOCAL EQ PLASTIC STRAIN  ***********
//    //*************************************************
//
//    int GetNumNonlocalEqPlasticStrain() const override;
//    int GetDofNonlocalEqPlasticStrain(int rComponent) const override;
//
//    const Eigen::Matrix<double, 2, 1>& GetNonlocalEqPlasticStrain() const override;
//    const Eigen::Matrix<double, 2, 1>& GetNonlocalEqPlasticStrain(int rTimeDerivative) const override;
//
//    void SetNonlocalEqPlasticStrain(const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain) override;
//    void SetNonlocalEqPlasticStrain(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain) override;
//
//    //*************************************************
//    //********    NONLOCAL TOTAL STRAIN     ***********
//    //*************************************************
//
//    int GetNumNonlocalTotalStrain()const override;
//    int GetDofNonlocalTotalStrain(int rComponent)const override;
//    double GetNonlocalTotalStrain(short rIndex)const override;
//
//    const Eigen::Matrix<double, 1, 1>& GetNonlocalTotalStrain1D() const override;
//    const Eigen::Matrix<double, 3, 1>& GetNonlocalTotalStrain2D() const override;
//    const Eigen::Matrix<double, 6, 1>& GetNonlocalTotalStrain3D() const override;
//    const Eigen::Matrix<double, Eigen::Dynamic, 1>& GetNonlocalTotalStrains() const override;
//
//    const Eigen::Matrix<double, 1, 1>& GetNonlocalTotalStrain1D(int rTimeDerivative) const override;
//    const Eigen::Matrix<double, 3, 1>& GetNonlocalTotalStrain2D(int rTimeDerivative) const override;
//    const Eigen::Matrix<double, 6, 1>& GetNonlocalTotalStrain3D(int rTimeDerivative) const override;
//    const Eigen::Matrix<double, Eigen::Dynamic, 1>& GetNonlocalTotalStrains(int rTimeDerivative) const override;
//
//    void SetNonlocalTotalStrain1D(const Eigen::Matrix<double, 1, 1>& rNonlocalTotalStrain) override;
//    void SetNonlocalTotalStrain2D(const Eigen::Matrix<double, 3, 1>& rNonlocalTotalStrain) override;
//    void SetNonlocalTotalStrain3D(const Eigen::Matrix<double, 6, 1>& rNonlocalTotalStrain) override;
//
//    void SetNonlocalTotalStrain1D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rNonlocalTotalStrain) override;
//    void SetNonlocalTotalStrain2D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rNonlocalTotalStrain) override;
//    void SetNonlocalTotalStrain3D(int rTimeDerivative, const Eigen::Matrix<double, 6, 1>& rNonlocalTotalStrain) override;
//
//    //*************************************************
//    //*******   NONLOCAL EQUIVALENT STRAIN   **********
//    //*************************************************
//
//    int GetNumNonlocalEqStrain() const override;
//    int GetDofNonlocalEqStrain() const override;
//
//    double GetNonlocalEqStrain() const override;
//    double GetNonlocalEqStrain(int rTimeDerivative) const override;
//
//    void SetNonlocalEqStrain(double rNonlocalEqStrain) override;
//    void SetNonlocalEqStrain(int rTimeDerivative, double rNonlocalEqStrain) override;
//
//    //*************************************************
//    //*******      WATER PHASE FRACTION      **********
//    //*************************************************
//
//    int GetNumWaterPhaseFraction() const override;
//    int GetDofWaterPhaseFraction() const override;
//
//    double GetWaterPhaseFraction() const override;
//    double GetWaterPhaseFraction(int rTimeDerivative) const override;
//
//    void SetWaterPhaseFraction(double rWaterPhaseFraction) override;
//    void SetWaterPhaseFraction(int rTimeDerivative, double rWaterPhaseFraction) override;
//
//    //*************************************************
//    //*******       RELATIVE HUMIDITY        **********
//    //*************************************************
//
//    int GetNumRelativeHumidity() const override;
//    int GetDofRelativeHumidity() const override;
//
//    double GetRelativeHumidity() const override;
//    double GetRelativeHumidity(int rTimeDerivative) const override;
//
//    void SetRelativeHumidity(double rRelativeHumidity) override;
//    void SetRelativeHumidity(int rTimeDerivative, double rRelativeHumidity) override;


    //! @brief returns the type of node as a string (all the data stored at the node)
    //! @return string
    //std::string GetNodeDofypeStr()const override;

    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    //Node::eNodeDofype GetNodeDofype()const;

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    NodeBase* Clone() const override;

protected:

//    std::array<Eigen::Matrix<double, TNumTemperatures, 1>,TNumTimeDerivatives+1> mTemperatures;
//    std::array<int, TNumTemperatures> mDofTemperatures;
//
//    std::array<Eigen::Matrix<double, TNumNonlocalEqPlasticStrain, 1>,TNumTimeDerivatives+1> mNonlocalEqPlasticStrain;
//    std::array<int, TNumNonlocalEqPlasticStrain> mDofNonlocalEqPlasticStrain;
//
//    std::array<Eigen::Matrix<double, TNumNonlocalTotalStrain, 1>,TNumTimeDerivatives+1> mNonlocalTotalStrain;
//    std::array<int, TNumNonlocalTotalStrain> mDofNonlocalTotalStrain;
//
//    std::array<Eigen::Matrix<double, TNumNonlocalEqStrain, 1>,TNumTimeDerivatives+1> mNonlocalEqStrain;
//    std::array<int, TNumNonlocalEqStrain> mDofNonlocalEqStrain;
//
//    std::array<Eigen::Matrix<double, TNumWaterPhaseFraction, 1>,TNumTimeDerivatives+1> mWaterPhaseFraction;
//    std::array<int, TNumWaterPhaseFraction> mDofWaterPhaseFraction;
//
//    std::array<Eigen::Matrix<double, TNumRelativeHumidity, 1>,TNumTimeDerivatives+1> mRelativeHumidity;
//    std::array<int, TNumRelativeHumidity> mDofRelativeHumidity;

};

}//namespace NuTo

#endif //NodeDof_DEF_H

