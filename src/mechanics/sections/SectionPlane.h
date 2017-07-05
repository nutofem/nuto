#pragma once

#include <memory>
#include "mechanics/sections/Section.h"

namespace NuTo
{

//! @brief Section for two-dimensional elements
class SectionPlane : public NuTo::Section
{
public:
    //! @brief Create a new instance of the plane section
    //! @param thickness Section thickness
    //! @param isPlaneStrain `true` corresponds to plane strain, `false` to plane stress
    static std::shared_ptr<SectionPlane> Create(double thickness, bool isPlaneStrain);

    //! @brief Get the section thickness
    //! @return Section thickness
    double GetThickness() const override;

    virtual bool IsPlaneStrain() const override;

private:
    //! @brief Print information about the section
    void Info(std::ostream& out) const override;

    //! @brief Constructor
    //! @param thickness Section thickness
    //! @param isPlaneStrain `true` corresponds to plane strain, `false` to plane stress
    SectionPlane(double thickness, bool isPlaneStrain);

    //! @brief Section thickness
    double mThickness;

    //! @brief **true** -> plane strain; **false** -> plane stress
    bool mIsPlaneStrain;
};
}
