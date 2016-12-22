// $Id$

#ifndef RELATIVEINTERFACEDISPLACEMENTS3D_H_
#define RELATIVEINTERFACEDISPLACEMENTS3D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
//! @brief ... three-dimensional relative interface displacements
/*!
 *  In the three-dimensional case the vector of relative interface displacements reads
 *  \f[\boldsymbol{\Delta u} = \begin{bmatrix} \Delta u_N \\ \Delta u_{T1} \\  \Delta u_{T2}  \end{bmatrix},\f]
 *  where \f$  \Delta u_N \f$ is the relative interface displacement in normal direction,
 *  and \f$  \Delta u_{T1}, \Delta u_{T2}  \f$ are the relative interface displacements in tangential direction.
 */
//! @author Stefan Eckardt, ISM
//! @date November 2009
class RelativeInterfaceDisplacements3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    RelativeInterfaceDisplacements3D();

    //! @brief kind of copy constructor
    RelativeInterfaceDisplacements3D(const RelativeInterfaceDisplacements3D& rOther);

    //! @brief ... get number of components
    //! @return ... number of components
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get three-dimensional relative interface displacements
    //! @return ... array of relative interface displacements
    //! @sa mRelativeInterfaceDisplacements
    const double* GetRelativeInterfaceDisplacements3D() const;

    //! @brief ... get three-dimensional relative interface displacements
    //! @param rRelativeInterfaceDisplacements ... array of relative interface displacements
    //! @sa mRelativeInterfaceDisplacements
    void GetRelativeInterfaceDisplacements(NuTo::RelativeInterfaceDisplacements3D& rRelativeInterfaceDisplacements) const;

    //! @brief ... set relative interface displacements
    //! @param rRelativeInterfaceDisplacements ... array of relative interface displacements
    //! @sa mRelativeInterfaceDisplacements
    void SetRelativeInterfaceDisplacements3D(const double* rRelativeInterfaceDisplacements);

    //! @brief ... print information about the strain object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const;

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
     *  The components of the relative interface displacements are stored in column major format \f$ \left[ \Delta u_N, \Delta u_{T1}, \Delta u_{T2} \right] \f$
     */
    double mRelativeInterfaceDisplacements[3];
};

}

#endif // RELATIVEINTERFACEDISPLACEMENTS3D_H_ 
