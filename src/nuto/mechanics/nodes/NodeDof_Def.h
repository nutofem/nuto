#pragma once


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
#include "nuto/mechanics/nodes/NodeTemperature.h"
#include "nuto/mechanics/nodes/NodeWaterVolumeFraction.h"
#include "nuto/mechanics/nodes/NodeDamage.h"

#include "nuto/mechanics/nodes/NodeEnum.h"

#include "nuto/math/FullMatrix_Def.h"
#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeComponentBase.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

#define NODE_DOF_TEMPLATE_PARAMETERS int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperature, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain, int TNumNonlocalEqStrain, int TNumWaterVolumeFraction, int TNumRelativeHumidity, int TNumDamage
#define NODE_DOF_TEMPLATE_INITIALIZATION TNumCoordinates, TNumTimeDerivatives, TNumDisplacements, TNumRotations, TNumTemperature, TNumNonlocalEqPlasticStrain, TNumNonlocalTotalStrain, TNumNonlocalEqStrain, TNumWaterVolumeFraction, TNumRelativeHumidity, TNumDamage

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for all nodes
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperature, int TNumNonlocalEqPlasticStrain, int TNumNonlocalTotalStrain, int TNumNonlocalEqStrain, int TNumWaterVolumeFraction, int TNumRelativeHumidity, int TNumDamage>
class NodeDof:
        public NodeCoordinates<TNumCoordinates>,
        public NodeDisplacements<TNumDisplacements, TNumTimeDerivatives>,
        public NodeRotations<TNumRotations, TNumTimeDerivatives>,
        public NodeNonlocalEqPlasticStrain<TNumNonlocalEqPlasticStrain, TNumTimeDerivatives>,
        public NodeNonlocalEqStrain<TNumNonlocalEqStrain, TNumTimeDerivatives>,
        public NodeRelativeHumidity<TNumRelativeHumidity,TNumTimeDerivatives>,
        public NodeTemperature<TNumTemperature, TNumTimeDerivatives>,
        public NodeWaterVolumeFraction<TNumWaterVolumeFraction,TNumTimeDerivatives>,
        public NodeDamage<TNumDamage,TNumTimeDerivatives>
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

        ar & BOOST_SERIALIZATION_NVP(mTemperature);
        ar & BOOST_SERIALIZATION_NVP(mDofTemperature);

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

    //! @brief sets the global dofs numbers for each dof type
    //! @param rDofNumbers ... map containing the dof type and the current number
    void SetGlobalDofsNumbers(std::map<Node::eDof, int>& rDofNumbers) override;

    //! @brief write dof values to the node (based on global dof number)
    //! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void SetGlobalDofValues(int rTimeDerivative, const FullVector<double,Eigen::Dynamic>& rActiveDofValues, const FullVector<double,Eigen::Dynamic>& rDependentDofValues) override;

    //! @brief extract dof values from the node (based on global dof number)
    //! @param  rTimeDerivative (set eq. displacements=0, velocities=1, accelerations=2
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void GetGlobalDofValues(int rTimeDerivative, FullVector<double,Eigen::Dynamic>& rActiveDofValues, FullVector<double,Eigen::Dynamic>& rDependentDofValues) const override;

    //! @brief write dof values to the node (based on global dof number)
    //! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
    //! @param rDofType ... specific dof type
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void SetGlobalDofValues(int rTimeDerivative, Node::eDof rDofType, const FullVector<double,Eigen::Dynamic>& rActiveDofValues, const FullVector<double,Eigen::Dynamic>& rDependentDofValues) override;

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
    //! @param rDofType ... specific dof type
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void GetGlobalDofValues(int rTimeDerivative, Node::eDof rDofType, FullVector<double,Eigen::Dynamic>& rActiveDofValues, FullVector<double,Eigen::Dynamic>& rDependentDofValues) const override;

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    void RenumberGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering) override;

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rDofType ... specific dof type
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    void RenumberGlobalDofs(Node::eDof rDofType, std::vector<int>& rMappingInitialToNewOrdering) override;

    //! @brief returns the number of time derivatives stored at the node
    //! @return number of derivatives
    int GetNumTimeDerivatives() const override;


    //! @brief returns the type of node as a string (all the data stored at the node)
    //! @return string
    std::string GetNodeTypeStr() const override;

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    NodeBase* Clone() const override;

private:

    //! @brief extracts the appropriate dof value from the active dof vector or the dependent dof vector
    //! @param rDofNumber ... dof number
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    //! @return dof value that corresponds to the rDofNumber
    inline double GetDofValueFromVector(int rDofNumber, const FullVector<double,Eigen::Dynamic>& rActiveDofValues, const FullVector<double,Eigen::Dynamic>& rDependentDofValues) const;

    //! @brief writes the appropriate dof value to the active dof vector or the dependent dof vector
    //! @param rDofNumber ... dof number
    //! @param rDofValue ... dof value
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    inline void WriteNodeValueToVector(int rDofNumber, double rDofValue, FullVector<double,Eigen::Dynamic>& rActiveDofValues, FullVector<double,Eigen::Dynamic>& rDependentDofValues) const;


};

}//namespace NuTo
