/*
 * BoundaryElement2D.h
 *
 *  Created on: 9 Jun 2015
 *      Author: ttitsche
 */

#pragma once

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
    virtual Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput) override;

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType() const override
    {
        return Element::BOUNDARYELEMENT2D;
    }

    //! @brief returns the local dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! @return local dimension
    int GetLocalDimension() const override
    {
        return 2;
    }

    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const override;

    //! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
    const BoundaryElement2D* AsBoundaryElement2D() const override;

    //! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
    BoundaryElement2D* AsBoundaryElement2D() override;

    //! @brief returns true, if the boundary conditions are fulfilled, post-processing
    bool IsBoundaryConditionFulfilled() const;

#ifdef ENABLE_VISUALIZE
    void Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList) override;
#endif // ENABLE_VISUALIZE

protected:
    //! @brief ... just for serialization
    BoundaryElement2D()
    {
    }

};

} /* namespace NuTo */

