#pragma once

#include "mechanics/sections/SectionTruss.h"

namespace NuTo
{

class SectionVariableTruss : public SectionTruss
{
public:
    //! @brief Create a truss section with a weak spot
    //! @param area Cross-section area outside of weakened range
    //! @param location Global coordinate where the weak spot is located
    //! @param extent Extent of the weakened cross-section
    //! @param reductionRatio Ratio of weakened area to full area
    static std::shared_ptr<SectionVariableTruss> Create(double area, double location, double extent,
                                                        double reductionRatio);

    virtual double GetArea(double coordinate) const override;

private:
    virtual void Info(std::ostream& out) const override;

    SectionVariableTruss(double area, double location, double extent, double reductionRatio);

    //! Global coordinate where the weak spot is located
    double mWeakSpotLocation;

    //! Extent of the weakened cross-section
    double mWeakSpotExtent;

    //! Ratio of weakened area to full area
    double mAlpha;
};
}
