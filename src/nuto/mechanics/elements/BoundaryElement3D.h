/*
 * BoundaryElement3D.h
 *
 *  Created on: 3 Aug 2015
 *      Author: vhirtham
 */

#ifndef BOUNDARYELEMENT3D_H_
#define BOUNDARYELEMENT3D_H_

#include "nuto/mechanics/elements/BoundaryElementBase.h"

namespace NuTo
{

class StructureBase;
class ConstitutiveTangentLocal3x1;
class Element3D;
class BoundaryElement3D: public BoundaryElementBase
{
public:
    BoundaryElement3D(const ElementBase* rBaseElement, int rSurfaceEdge);

    virtual ~BoundaryElement3D()
    {

    }

    //! @brief calculates output data for the element
    //! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    //!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
    virtual Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput) override;

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType() const override
    {
        return Element::BOUNDARYELEMENT3D;
    }

    //! @brief returns the local dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! @return local dimension
    int GetLocalDimension() const override
    {
        return 3;
    }
    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const override;

    //! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
    const BoundaryElement3D* AsBoundaryElement3D() const override;

    //! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
    BoundaryElement3D* AsBoundaryElement3D() override;

    //! @brief returns true, if the boundary conditions are fulfilled, post-processing
    bool IsBoundaryConditionFulfilled() const;

protected:
    //! @brief ... just for serialization
    BoundaryElement3D()
    {
    }


};

} /* namespace NuTo */

#endif /* BOUNDARYELEMENT3D_H_ */
