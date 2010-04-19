/*
 * ElementDataStaticData.h
 *
 *  Created on: Apr 16, 2010
 *      Author: unger3
 */

#ifndef ELEMENTDATASTATICDATA_H_
#define ELEMENTDATASTATICDATA_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <boost/ptr_container/ptr_vector.hpp>
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for elements with a single material per element static data for each integration point
class ElementDataStaticDataBase : public virtual ElementDataBase
{

#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
	ElementDataStaticDataBase(const NuTo::IntegrationTypeBase* rIntegrationType);

	virtual ~ElementDataStaticDataBase();

    ConstitutiveStaticDataBase* GetStaticData(int rIp);

    const ConstitutiveStaticDataBase* GetStaticData(int rIp)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataBase)
           & BOOST_SERIALIZATION_NVP(mStaticData);
    }
#endif  // ENABLE_SERIALIZATION

protected:
     // if the class is extended (e.g. for nonlocal data, plot data, stresses or strains at the integration point
    // please use a vector<IP_data >, where IP_data is a struct comprising all the required data
    // the size of mStaticData corresponds to the number of integration points, so for each integration point a static data object ist stored
    // no ptr_vector is used since it does not (in general) support null pointer
    std::vector<ConstitutiveStaticDataBase*> mStaticData;
};
}
#endif /* ELEMENTDATASTATICDATA_H_ */
