// $Id$ 
#ifndef IPDATABASE_H_
#define IPDATABASE_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include <vector>

namespace NuTo
{
class ConstitutiveBase;
class ConstitutiveStaticDataBase;
class ElementBase;
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public :

	virtual ~IpDataBase();

	virtual void Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive)=0;

	//! @brief adds the weight to an integration point, eventually reallocates the data
	//! @param rNonlocalElement the Element (local number from the nonlocal elements)
	//! @param rNonlocalIp integration point of the nonlocal element
	//! @param rNumIps number of integration points of the nonlocal element (for allocation purpose of not existing)
	//! @param rWeight nonlocal weight
	virtual void SetNonlocalWeight(int rNonlocalElement,int rNonlocalIp, int rNumIps, double rWeight);

	//! @brief return the nonlocal weights
	//! @param rNonlocalElement nonlocal element (between 0 and nonlocal elements.size stored in nonlocal element data)
	//! @return nonlocal weights
	virtual const std::vector<double>& GetNonlocalWeights(int rNonlocalElement)const;

	virtual ConstitutiveStaticDataBase* GetStaticData();

    virtual const ConstitutiveStaticDataBase* GetStaticData()const;

    virtual void SetStaticData(ConstitutiveStaticDataBase* rStaticData);

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleModel(std::string rFileName);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

};
}
#endif /* IPDATABASE_H_ */
