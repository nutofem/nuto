#pragma once

#include <memory>
#include "mechanics/sections/Section.h"

namespace NuTo
{

//! @brief Section for the fibre-matrix bond element
class SectionFibreMatrixBond : public NuTo::Section
{
public:
    //! @brief Create a new instance of the fibre-matrix bond section
    static std::shared_ptr<SectionFibreMatrixBond> Create(double circumference);

    //! @brief Get the circumference of the fibre
    virtual double GetCircumference() const override;

private:
    virtual void Info(std::ostream& out) const override;

    //! @brief Constructor
    //! @param circumference Circumference of the fibre
    SectionFibreMatrixBond(double circumference);

    //! @brief Circumference of the fibre
    double mCircumference;
};
}
