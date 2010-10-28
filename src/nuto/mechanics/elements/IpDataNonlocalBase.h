// $Id$ 
#ifndef IPDATANONLOCALBASE_H_
#define IPDATANONLOCALBASE_H_

#include <vector>
#include "nuto/mechanics/elements/IpDataBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataNonlocalBase : public virtual IpDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	IpDataNonlocalBase() : IpDataBase()
	{}

	virtual ~IpDataNonlocalBase();

	//! @brief adds the weight to an integration point, eventually reallocates the data
	//! @param rNonlocalElement the Element (local number from the nonlocal elements)
	//! @param rNonlocalIp integration point of the nonlocal element
	//! @param rNumIps number of integration points of the nonlocal element (for allocation purpose of not existing)
	//! @param rWeight nonlocal weight
	void SetNonlocalWeight(int rElement,int rNonlocalIp,int rNumIps, double rWeight);


	//! @brief return the nonlocal weights
	//! @param rNonlocalElement nonlocal element (between 0 and nonlocal elements.size stored in nonlocal element data)
	//! @return nonlocal weights
	virtual const std::vector<double>& GetNonlocalWeights(int rNonlocalElement)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
    //! @brief for each nonlocal element (local number stored at the element level) the nonlocal weights (for each integration point of the nonlocal element)
	std::vector<std::vector<double> > mNonlocalWeights;
};
}
#endif /* IPDATANONLOCALBASE_H_ */
