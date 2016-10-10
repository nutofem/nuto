// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION


namespace NuTo
{
class SectionTruss;
enum class eSectionType;

//! @author Stefan Eckardt, ISM
//! @date October 2009
//! @brief ... standard abstract base class for sections
class SectionBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    SectionBase();

    //! @brief ... destructor
    virtual ~SectionBase(){}

    //! @brief ... get the section type
    //! @return ... section type
    virtual eSectionType GetType() const = 0;

    //! @brief ... get the cross-section area of the section
    //! @return ... section cross-section area
    virtual double GetArea() const;

    //! @brief ... set the cross-section area of the section
    //! @param rArea ... cross-section area
    virtual void SetArea(double rArea);

    //! @brief ... get the thickness of the section
    //! @return ... section thickness
    virtual double GetThickness() const;

    //! @brief ... set the thickness of the section
    //! @return ... section thickness
    virtual void SetThickness(double rThickness);

    //! @brief ... get the circumference of the section
    //! @return ... section circumference
    virtual double GetCircumference() const;

    //! @brief ... set the circumference of the section
    //! @return ... section circumference
    virtual void SetCircumference(double rCircumference);

    //! @brief ... print information about the section
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const;

    //! @brief ... cast the base pointer to an SectionTruss, otherwise throws an exception
    virtual SectionTruss* AsSectionTruss();

    //! @brief ... cast the base pointer to an SectionTruss, otherwise throws an exception
    virtual const SectionTruss* AsSectionTruss() const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

#endif // ENABLE_SERIALIZATION

};

}// namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::SectionBase)
#endif // ENABLE_SERIALIZATION

