/*
 * CallbackHandlerBase.h
 *
 *  Created on: 24 Mar 2016
 *      Author: ttitsche
 */

#pragma once

#include "nuto/base/Exception.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{

class StructureBase;

//! @author Thomas Titscher, BAM
//! @date March 2016
//! @brief ... abstract class to handle callback routines
class CallbackInterface
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

public:
    CallbackInterface() = default;
    virtual ~CallbackInterface() = default;


    //! @brief Exit function, returns true if a certain state of the structure is reached
    //! @param StructureBase& ... not const since most of the useful methods require an Evaluate, which is non-const
    virtual bool Exit(StructureBase&) const
    {
        throw Exception(__PRETTY_FUNCTION__, "not implemented for this callback.");
    }
};

} /* namespace NuTo */
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::CallbackInterface)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::CallbackInterface)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
