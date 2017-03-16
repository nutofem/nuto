#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <memory>
#include "mechanics/sections/Section.h"

namespace NuTo
{

//! @brief One-dimensional section
class SectionTruss : public Section
{
#ifdef ENABLE_SERIALIZATION
   friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief Create a new instance of the truss section
    //! @param area Cross-section area
    static std::shared_ptr<SectionTruss> Create(double area);

    virtual double GetArea(double coordinate) const override;

    virtual void Info() const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

#endif // ENABLE_SERIALIZATION
protected:
    //! @brief Constructor
    //! @param area Cross-section area
   SectionTruss(double area);

    //! @brief Cross-section area
    double mArea;
};

} // namespace

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::SectionTruss)
#endif // ENABLE_SERIALIZATION



