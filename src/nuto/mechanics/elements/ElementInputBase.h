/*
 * ElementInputBase.h
 *
 *  Created on: 11 May 2015
 *      Author: ttitsche
 */

#ifndef ELEMENTINPUTBASE_H_
#define ELEMENTINPUTBASE_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif// ENABLE_SERIALIZATION

namespace NuTo
{

class ElementInputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ElementInputBase();
    virtual ~ElementInputBase();
};
} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementInputBase)
#endif // ENABLE_SERIALIZATION

#endif /* ELEMENTINPUTBASE_H_ */
