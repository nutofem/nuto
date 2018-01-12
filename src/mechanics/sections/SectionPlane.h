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

    //! @brief Create an axisymmetric section
    //! @param thickness is set to 1
    static std::shared_ptr<SectionPlane> CreateAxiSymmetric();

    //! @brief Get the section thickness
    //! @return Section thickness
    double GetThickness() const override;

    virtual bool IsPlaneStrain() const override;

    virtual bool IsAxiSymmetric() const override;

private:
    //! @brief Print information about the section
    void Info(std::ostream& out) const override;

    //! @brief Constructor
    //! @param thickness Section thickness
    //! @param isPlaneStrain `true` corresponds to plane strain, `false` to plane stress
    SectionPlane(double thickness, bool isPlaneStrain);

    //! @brief Constructor for axisymmetric
    //! @param thickness = 1
    SectionPlane();

    //! @brief Section thickness
    double mThickness;

    //! @brief **true** -> plane strain; **false** -> plane stress
    bool mIsPlaneStrain;

    //! @brief **true** -> axisymmetric; **false** -> plane train or plane stress
    bool mIsAxiSymmetric;
};
}
