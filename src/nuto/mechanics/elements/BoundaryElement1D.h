/*
 * BoundaryElement1D.h
 *
 *  Created on: 4 Jun 2015
 *      Author: ttitsche
 */

#ifndef BOUNDARYELEMENT1D_H_
#define BOUNDARYELEMENT1D_H_

#include "nuto/mechanics/elements/BoundaryElementBase.h"

namespace NuTo
{

class StructureBase;
class ConstitutiveTangentLocal1x1;
class Element1D;
class BoundaryElement1D: public BoundaryElementBase
{
public:
    BoundaryElement1D(const ElementBase* rBaseElement, int rSurfaceEdge);

    virtual ~BoundaryElement1D()
    {

    }

    //! @brief calculates output data for the element
    //! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    //!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
    Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput) override;

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType() const override
    {
        return Element::BOUNDARYELEMENT1D;
    }

    //! @brief returns the local dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! @return local dimension
    int GetLocalDimension() const override
    {
        return 1;
    }

    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const override;

    // getter setter

    double GetBoundaryRelativeHumidity() const;
    void SetBoundaryRelativeHumidity(double rBoundaryRelativeHumidity);

    double GetBoundaryWaterVolumeFraction() const;
    void SetBoundaryWaterVolumeFraction(double rBoundaryWaterVolumeFraction);

    //! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
    const BoundaryElement1D* AsBoundaryElement1D() const override;

    //! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
    BoundaryElement1D* AsBoundaryElement1D() override;

protected:
    //! @brief ... just for serialization
    BoundaryElement1D()
    {
    }

    //! @brief ... relative humidity at the boundary surface
    double mBoundaryRelativeHumidity;

    //! @brief ... water volume fraction at the boundary surface
    double mBoundaryWaterVolumeFraction;

};

} /* namespace NuTo */

#endif /* BOUNDARYELEMENT1D_H_ */
