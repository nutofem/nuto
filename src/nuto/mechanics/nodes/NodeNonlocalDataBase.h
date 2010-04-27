#ifndef NODE_NONLOCALDATABASE_H
#define NODE_NONLOCALDATABASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for reference nodes
class NodeNonlocalDataBase : public virtual NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief constructor
    NodeNonlocalDataBase();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the number of nonlocal nodes
    //! @return number of nonlocal nodes
    int GetNumNonlocalNodes()const;

    //! @brief set the nonlocal data
    //! @param rNodes  nonlocal nodes
    //! @param rWeights  nonlocal weights
    //! @param rConstitutive  corresponding constitutive model
    void SetNonlocalData(const std::vector<const NodeBase*>& rNodes,
    		const std::vector<double>& rWeights, const ConstitutiveBase* rConstitutive);

    //! @brief returns the nonlocal nodes
    //! @param rConstitutive  corresponding constitutive model
    //! @return nonlocal nodes
    const std::vector<const NodeBase*>& GetNonlocalNodes(const ConstitutiveBase* rConstitutive)const;

    //! @brief returns the nonlocal weights
    //! @param rConstitutive  corresponding constitutive model
    //! @return nonlocal weights
    const std::vector<double>& GetNonlocalWeights(const ConstitutiveBase* rConstitutive)const;


protected:
    std::vector<const NodeBase*> mNodes; //!< nonlocal nodes
    std::vector<double> mWeights;       //!< nonlocal weights bar(u)=sum_i mWeights_i*u(mNodes_i)
    const ConstitutiveBase* mConstitutive;    //!< this is based on the assumption of a nonlocal formulation for a single material, for more, just make a map of mNodes and mWeights for each constituent
};
}//namespace NuTo
#endif //NODE_NONLOCALDATABASE_H
