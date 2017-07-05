#pragma once

#include <memory>
#include "mechanics/sections/Section.h"

namespace NuTo
{

//! @brief One-dimensional section
class SectionTruss : public Section
{
public:
    //! @brief Create a new instance of the truss section
    //! @param area Cross-section area
    static std::shared_ptr<SectionTruss> Create(double area);

    virtual double GetArea(double coordinate) const override;

protected:
    virtual void Info(std::ostream& out) const override;

    //! @brief Constructor
    //! @param area Cross-section area
   SectionTruss(double area);

    //! @brief Cross-section area
    double mArea;
};

} // namespace



