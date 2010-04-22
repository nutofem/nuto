// $ld: $ 
// NodeNonlocalDataBase.cpp
// created Apr 22, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include "nuto/mechanics/nodes/NodeNonlocalDataBase.h"
#include "assert.h"

//! @brief constructor
NuTo::NodeNonlocalDataBase::NodeNonlocalDataBase() : NodeBase::NodeBase()
{
   	mConstitutive = 0;
}



#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::NodeNonlocalDataBase::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
           & BOOST_SERIALIZATION_NVP(mNodes)
           & BOOST_SERIALIZATION_NVP(mWeights)
           & BOOST_SERIALIZATION_NVP(mConstitutive);
    }
#endif // ENABLE_SERIALIZATION

    //! @brief returns the number of nonlocal nodes
    //! @return number of nonlocal nodes
    int NuTo::NodeNonlocalDataBase::GetNumNonlocalNodes()const
    {
    	assert(mNodes.size()==mWeights.size());
    	return mNodes.size();
    }

    //! @brief set the nonlocal data
    //! @param rNodes  nonlocal nodes
    //! @param rWeights  nonlocal weights
    //! @param rConstitutive  corresponding constitutive model
    void NuTo::NodeNonlocalDataBase::SetNonlocalData(const std::vector<const NodeBase*>& rNodes,
    		const std::vector<double>& rWeights, const ConstitutiveBase* rConstitutive)
    {
    	mNodes = rNodes;
    	mWeights = rWeights;
    	mConstitutive = rConstitutive;
    }

    //! @brief returns the nonlocal nodes
    //! @param rConstitutive  corresponding constitutive model
    //! @return nonlocal nodes
    const std::vector<const NuTo::NodeBase*>& NuTo::NodeNonlocalDataBase::GetNonlocalNodes(const ConstitutiveBase* rConstitutive)const
	{
    	if (rConstitutive==mConstitutive)
    	    return mNodes;
    	else
    		throw MechanicsException("[NuTo::NodeNonlocalDataBase::GetNonlocalNodes] For this constitutive model no nonlocal data has been calculated.");
	}

    //! @brief returns the nonlocal weights
    //! @param rConstitutive  corresponding constitutive model
    //! @return nonlocal weights
    const std::vector<double>& NuTo::NodeNonlocalDataBase::GetNonlocalWeights(const ConstitutiveBase* rConstitutive)const
    {
    	if (rConstitutive==mConstitutive)
    	    return mWeights;
    	else
    		throw MechanicsException("[NuTo::NodeNonlocalDataBase::GetNonlocalNodes] For this constitutive model no nonlocal data has been calculated.");
    }
