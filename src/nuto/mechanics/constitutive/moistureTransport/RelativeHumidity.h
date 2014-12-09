#ifndef RELATIVEHUMIDITY_H
#define RELATIVEHUMIDITY_H

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
//! @brief ... relative humidity
//! @author Volker Hirthammer, BAM
//! @date December 2014
class RelativeHumidity : public ConstitutiveInputBase, public FullVector<double,1>
{
public:
    RelativeHumidity();

    virtual const RelativeHumidity& GetRelativeHumidity() const override;
};

} // Namespace NuTo



#endif // RELATIVEHUMIDITY_H
