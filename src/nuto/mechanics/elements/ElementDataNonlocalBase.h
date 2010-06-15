// $ld: $ 
#ifndef ELEMENTDATANONLOCALBASE_H_
#define ELEMENTDATANONLOCALBASE_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/ElementDataBase.h"

#ifdef ENABLE_SERIALIZATION
   #include <boost/serialization/vector.hpp>
#else
   #include <vector>
#endif  // ENABLE_SERIALIZATION

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 22, 2010
//! @brief class to store the nonlocal weights for regularized constitutive models
class ElementDataNonlocalBase : public virtual ElementDataBase
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
	ElementDataNonlocalBase();

	virtual ~ElementDataNonlocalBase();

	const std::vector<const ElementBase*>& GetNonlocalElements(const ConstitutiveBase* rConstitutive)const;

    //! @brief adds an element to the nonlocal elements
    //! @param rConstitutive  constitutive model
    //! @return the local element number, the element is either append to the list, or the existing local number is returned
    int AddNonlocalElement(const ElementBase* rElement, const ConstitutiveBase* rConstitutive);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataBase)
           & BOOST_SERIALIZATION_NVP(mConstitutive);
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief adds the nonlocal weight to an integration point
    //! @param rLocalIpNumber local Ip
    //! @param rConstitutive constitutive model for which nonlocal data is to be calculated
    //! @param rNonlocalElement element of the nonlocal ip
    //! @param rNonlocalIp local ip number of the nonlocal ip
    //! @param rWeight weight
//    virtual void AddNonlocalIp(int rLocalIpNumber, const ConstitutiveBase* rConstitutive,
//    		const ElementWithDataBase* rNonlocalElement, int rNonlocalIp, double rWeight);

protected:
    std::vector<const ElementBase*> mNonlocalElements;
    //! @brief constitutive model for which the nonlocal elements are build, for different const. models in one model, just introduce a map<const mConstitutive,std::vector<const ElementBase*> >
    const ConstitutiveBase* mConstitutive;
};
}
#endif /* ELEMENTDATANONLOCALBASE_H_ */
