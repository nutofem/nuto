#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <memory>
#include "mechanics/sections/Section.h"

namespace NuTo
{

//! @brief Section for the fibre-matrix bond element
class SectionFibreMatrixBond : public NuTo::Section
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief Create a new instance of the fibre-matrix bond section
    static std::shared_ptr<SectionFibreMatrixBond> Create(double circumference);

    //! @brief Get the circumference of the fibre
    virtual double GetCircumference() const override;

    virtual void Info() const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

private:
    //! @brief Constructor
    //! @param circumference Circumference of the fibre
    SectionFibreMatrixBond(double circumference);

    //! @brief Circumference of the fibre
    double mCircumference;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::SectionFibreMatrixBond)
#endif // ENABLE_SERIALIZATION
