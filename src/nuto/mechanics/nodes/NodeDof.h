#pragma once

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/DofHash.h"
#include <unordered_map>


namespace NuTo
{
struct NodeDofInfo
{
    int mDimension;
    int mNumTimeDerivatives;
    bool mIsDof;
};


//! @author Thomas Titscher, BAM
//! @date July 2016
//! @brief ... standard class for all nodes
class NodeDof: public NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    NodeDof() {}
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ctor
    //! @param contains the information for each dof type
    NodeDof(std::map<Node::eDof, NodeDofInfo> rDofInfos);

    //! @brief sets the global dofs numbers for each dof type
    //! @param rDofNumbers ... map containing the dof type and the current number
    void SetGlobalDofsNumbers(std::map<Node::eDof, int>& rDofNumbers) override;

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
    //! @param rDofType ... specific dof type
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    void RenumberGlobalDofs(Node::eDof rDofType, std::vector<int>& rMappingInitialToNewOrdering) override;

    //! @brief returns the number of time derivatives stored at the node
    //! @param rDof ... specific dof type
    //! @return number of derivatives
    int GetNumTimeDerivatives(Node::eDof rDof) const override;

    //! @brief returns if the dof type rDof is an actual degree of freedom (in contrast to COORDINATES)
    bool IsDof(Node::eDof rDof) const override;

    //*************************************************
    //************       ACCESS         ***************
    //*************************************************

    //! @brief returns the total number of dofs
    int GetNumDofs() const override;

    //! @brief returns the number of dofs for a specific dof type
    //! @param rDof ... specific dof type
    int GetNum(Node::eDof rDof) const override;

    //! @brief returns the global dof number for a specific dof type
    //! @param rDof ... specific dof type
    //! @param rComponent ... component index of the dof type (e.g. 0,1,2 for coordinates in 3D)
    int GetDof(Node::eDof rDof, int rComponent) const override;

    //! @brief returns the dof values for a specific dof
    //! @param rDof ... specific dof type
    //! @param rTimeDerivative ... time derivative
    const Eigen::VectorXd& Get(Node::eDof rDof, int rTimeDerivative) const override;

    //! @brief sets the dof values for a specific dof
    //! @param rDof ... specific dof type
    //! @param rTimeDerivative ... time derivative
    //! @param rValue ... dof value
    void Set(Node::eDof rDof, int rTimeDerivative , const Eigen::VectorXd& rValue) override;

    //! @brief returns a set containing all dof types
    std::set<Node::eDof> GetDofTypes() const override;

    //! @brief returns the type of node as a string (all the data stored at the node)
    //! @return string
    std::string GetNodeTypeStr() const override;

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    NodeBase* Clone() const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION


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

    //! @brief stores the dof values (std::vector for time derivatives, VectorXd for values)
    std::unordered_map<Node::eDof, std::vector<Eigen::VectorXd>, Node::eDofHash> mDofValues;

    //! @brief stores the global dof numbers in a dynamic integer vector
    std::unordered_map<Node::eDof, Eigen::VectorXi, Node::eDofHash> mDofNumbers;

};

} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NodeDof)
#endif // ENABLE_SERIALIZATION
