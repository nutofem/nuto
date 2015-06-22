/*
 * BoundaryElement2D.h
 *
 *  Created on: 9 Jun 2015
 *      Author: ttitsche
 */

#ifndef BOUNDARYELEMENT2D_H_
#define BOUNDARYELEMENT2D_H_

#include "nuto/mechanics/elements/BoundaryElementBase.h"

namespace NuTo
{

class StructureBase;
class ConstitutiveTangentLocal3x1;
class Element2D;
class BoundaryElement2D: public BoundaryElementBase
{
public:
    BoundaryElement2D(const ElementBase* rBaseElement, int rSurfaceEdge);

    virtual ~BoundaryElement2D()
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
        return Element::BOUNDARYELEMENT2D;
    }

    //! @brief returns the global dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! there is also a routine GetLocalDimension, which is e.g. 2 for plane elements and 1 for truss elements
    //! @return global dimension
    int GetGlobalDimension() const override
    {
        return 2;
    }


    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const override;

    // getter setter

    BoundaryType::eType GetBoundaryConditionType() const;
    void SetBoundaryConditionType(BoundaryType::eType rBoundaryConditionType);

    //! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
    const BoundaryElement2D* AsBoundaryElement2D() const override;

    //! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
    BoundaryElement2D* AsBoundaryElement2D() override;

protected:
    //! @brief ... just for serialization
    BoundaryElement2D()
    {
    }

};

} /* namespace NuTo */

#endif /* BOUNDARYELEMENT1D_H_ */
