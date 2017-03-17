#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <memory>
#include "mechanics/sections/Section.h"

namespace NuTo
{

//! @brief Section for two-dimensional elements
class SectionPlane : public NuTo::Section
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief Create a new instance of the plane section
    //! @param thickness Section thickness
    //! @param isPlaneStrain `true` corresponds to plane strain, `false` to plane stress
    static std::shared_ptr<SectionPlane> Create(double thickness, bool isPlaneStrain);

    //! @brief Get the section thickness
    //! @return Section thickness
    double GetThickness() const override;

    virtual bool IsPlaneStrain() const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

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
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::SectionPlane)
#endif // ENABLE_SERIALIZATION
