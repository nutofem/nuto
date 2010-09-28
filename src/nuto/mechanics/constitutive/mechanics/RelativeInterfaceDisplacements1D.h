// $Id$

#ifndef RELATIVEINTERFACEDISPLACEMENTS1D_H_
#define RELATIVEINTERFACEDISPLACEMENTS1D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
//! @brief ... one-dimensional relative interface displacements
/*!
 *  In the one-dimensional case the vector of relative interface displacements reads
 *  \f[\boldsymbol{\Delta u} = \begin{bmatrix} \Delta u_N \end{bmatrix},\f]
 *  where \f$  \Delta u_N \f$ is the relative interface displacement in normal direction.
 */
//! @author Stefan Eckardt, ISM
//! @date November 2009
class RelativeInterfaceDisplacements1D
{

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    RelativeInterfaceDisplacements1D();

    //! @brief kind of copy constructor
    RelativeInterfaceDisplacements1D(const RelativeInterfaceDisplacements1D& rOther);

    //! @brief ... get number of components
    //! @return ... number of components
    virtual unsigned int GetNumberOfComponents() const;

    //! @brief ... get one-dimensional relative interface displacements
    //! @return ... array of relative interface displacements
    //! @sa mRelativeInterfaceDisplacements
    virtual const double* GetRelativeInterfaceDisplacements1D() const;

    //! @brief ... get one-dimensional relative interface displacements
    //! @param rRelativeInterfaceDisplacements ... array of relative interface displacements
    //! @sa mRelativeInterfaceDisplacements
    virtual void GetRelativeInterfaceDisplacements(NuTo::RelativeInterfaceDisplacements1D& rRelativeInterfaceDisplacements) const;

    //! @brief ... set one-dimensional relative interface displacements
    //! @param rRelativeInterfaceDisplacements ... array of relative interface displacements
    //! @sa mRelativeInterfaceDisplacements
    virtual void SetRelativeInterfaceDisplacements1D(const double* rRelativeInterfaceDisplacements);

    //! @brief ... print information about the strain object
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... array of relative interface displacements
    /*!
     *  The components of the relative interface displacements are stored in column major format \f$ \left[ \Delta u_N \right] \f$
     */
    double mRelativeInterfaceDisplacements;
};

}

#endif // RELATIVEINTERFACEDISPLACEMENTS1D_H_ 
