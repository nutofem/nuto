// $Id$

#pragma once

#include "mechanics/sections/SectionBase.h"

namespace NuTo
{
class SectionSpring : public SectionBase
{

public:
    //! @brief ... constructor
   SectionSpring();

    //! @brief ... get the section type
    //! @return ... section type
    virtual eSectionType GetType() const override;

    //! @brief ... print information about the section
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const override;

};

}
