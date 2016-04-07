/*
 * ImplEx.h
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/timeIntegration/ImplicitExplicitBase.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

namespace NuTo
{

class ImplEx: public ImplicitExplicitBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#ifndef SWIG
    template<class Archive> void serialize(Archive & ar, const unsigned int version);
    ImplEx() = default;
#endif// SWIG

#endif  // ENABLE_SERIALIZATION
public:

    ImplEx(StructureBase* rStructure);


    virtual ~ImplEx() = default;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId() const override
    {
        return "ImplEx";
    }
};

} /* namespace NuTo */

