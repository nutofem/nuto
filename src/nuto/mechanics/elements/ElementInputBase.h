/*
 * ElementInputBase.h
 *
 *  Created on: 11 May 2015
 *      Author: ttitsche
 */

#ifndef ELEMENTINPUTBASE_H_
#define ELEMENTINPUTBASE_H_

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

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementInputBase)
#endif // ENABLE_SERIALIZATION

} /* namespace NuTo */

#endif /* ELEMENTINPUTBASE_H_ */
