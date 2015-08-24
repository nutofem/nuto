// $Id$

#ifndef SECTIONSPRING_H
#define SECTIONSPRING_H

#include "nuto/mechanics/sections/SectionBase.h"

namespace NuTo
{
class SectionSpring : public SectionBase
{

public:
    //! @brief ... constructor
   SectionSpring();

    //! @brief ... get the section type
    //! @return ... section type
    virtual Section::eSectionType GetType() const;

    //! @brief ... print information about the section
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const;

};

}

#endif // SECTIONSPRING_H
