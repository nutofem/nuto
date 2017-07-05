/*
 * CallbackHandlerBase.h
 *
 *  Created on: 24 Mar 2016
 *      Author: ttitsche
 */

#pragma once


// Other:


namespace NuTo
{

class StructureBase;

//! @author Thomas Titscher, BAM
//! @date March 2016
//! @brief ... abstract class to handle callback routines
class CallbackInterface
{

public:
    CallbackInterface() = default;
    virtual ~CallbackInterface() = default;


    //! @brief Exit function, returns true if a certain state of the structure is reached
    //! @param StructureBase& ... not const since most of the useful methods require an Evaluate, which is non-const
    virtual bool Exit(StructureBase&);
};

} /* namespace NuTo */
