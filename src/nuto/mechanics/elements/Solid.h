// $Id$
#ifndef SOLID_H
#define SOLID_H

#include "nuto/mechanics/elements/ElementBase.h"

namespace NuTo
{
class StructureBase;
class DeformationGradient3D;
class TemperatureGradient3D;
class ConstitutiveTangentLocal3x3;
class ConstitutiveTangentLocal6x6;
class EngineeringStress3D;
class HeatFlux3D;
template <int TNumRows, int TNumColumns> class ConstitutiveTangentLocal;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all isoparametric displacement based finite elements in 3D
class Solid : public ElementBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    Solid(const StructureBase* rStructure,
    		ElementData::eElementDataType rElementDataType,
    		IntegrationType::eIntegrationType rIntegrationType,
    		IpData::eIpDataType rIpDataType
    		);

    //! @brief returns the global dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! there is also a routine GetLocalDimension, which is e.g. 2 for plane elements and 1 for truss elements
    //! @return global dimension
    int GetGlobalDimension()const
    {
        return 3;
    }

    //! @brief calculates output data fo the elmement
    //! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    //!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
    Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rConstitutiveOutput);


    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the stiffness matrix
    //! @param rResult ... coefficient matrix
    //! @param rGlobalDofsRow ... row numbers in global system
    //! @param rGlobalDofsColumn ... column numbers in global system
    //! @param rSymmetry ... matrix is symmetric or not (in the symmetric case the full matrix is also stored
    //Error::eError CalculateCoefficientMatrix_0(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult,
    //        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const;

    //! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the damping matrix
    //Error::eError CalculateCoefficientMatrix_1(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult,
    //        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const;

    //! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the Mass matrix
    //Error::eError CalculateCoefficientMatrix_2(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult,
    //        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const;

    //! @brief calculates the gradient of the internal potential
    //! for a mechanical problem, this corresponds to the internal force vector
    //Error::eError CalculateGradientInternalPotential(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult,
    //                                        std::vector<int>& rGlobalDofs)const;

    //! @brief sets the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @param rSection pointer to section
    //! @return pointer to constitutive law
    void SetSection(const SectionBase* rSection);

    //! @brief returns a pointer to the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @return pointer to section
    const SectionBase* GetSection()const;

    //! @brief Update the static data of an element
    //Error::eError UpdateStaticData(NuTo::Element::eUpdateType rUpdateType);

    //! @brief stores the local coordinates of the nodes
    //! @param localCoordinates vector with already correct size allocated
    //! this can be checked with an assertation
    void CalculateCoordinates(std::vector<double>& rCoordinates)const;

    //! @brief stores the local displacements of the nodes
    //! @param localDisplacements vector with already correct size allocated
    //! this can be checked with an assertation
    void CalculateDisplacements(std::vector<double>& rDisplacements)const;

    //! @brief stores the temperatures of the nodes
    //! @param temperature vector with already correct size allocated
    //! this can be checked with an assertation
    void CalculateTemperatures(std::vector<double>& rTemperatures)const;

    //! @brief calculates the deformation gradient in 3D
    //! @param rRerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rLocalDisp local displacements
    //! @param rDeformationGradient (return value)
    void CalculateDeformationGradient(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                                      const std::vector<double>& rDisp,
                                      DeformationGradient3D& rDeformationGradient)const;

    //! @brief calculates the temperature gradient in 3D
    //! @param rRerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rTemp nodal temperatures
    //! @param rTemperatureGradient (return value)
    void CalculateTemperatureGradient(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
            const std::vector<double>& rTemp,
            TemperatureGradient3D& rTemperatureGradient)const;

    //! @brief Calculates the the inverse of the Jacobian and its determinant
    //! @param rDerivativeShapeFunctions Derivatives of the shape functions (dN1dx, dN1dy, dN1dz, dN2dx, ..
    //! @param rNodeCoordinates Node coordinates (X1,Y1,Z1,X2,Y2,Z2,...
    //! @param rInvJacobian inverse Jacobian matrix (return value)
    //! @param rDetJac determinant of the Jacobian (return value)
    void CalculateJacobian(const std::vector<double>& rDerivativeShapeFunctions,
                           const std::vector<double>& rNodeCoordinates,
                           double rInvJacobian[9],
                           double& rDetJac)const;

    //! @brief returns the number of shape functions
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    virtual int GetNumShapeFunctions()const=0;

    //! @brief returns the coordinates of an integration point
    //! @param rIpNum integration point
    //! @param rCoordinates coordinates to be returned
    virtual void GetLocalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const;

    //! @brief returns the coordinates of an integration point
    //! @param rIpNum integration point
    //! @param rCoordinates coordinates to be returned
    virtual void GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const;


    //! @brief calculates the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
    virtual void CalculateShapeFunctions(const double rLocalCoordinates[3], std::vector<double>& rShapeFunctions)const=0;

    //! @brief calculates the derivative of the shape functions with respect to local coordinates
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsLocal(const double rLocalCoordinates[3], std::vector<double>& rDerivativeShapeFunctions)const=0;

    //! @brief calculates the derivative of the shape functions with respect to global coordinates
    //! @param std::vector<double>& rDerivativeShapeFunctions derivatives of the shape functions
    //! @param rJacInv inverse of the Jacobian
    //! @param rDerivativeShapeFunctionsGlobal derivaties of the shape functions with respect to global coordinates
    //! size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsGlobal(const std::vector<double>& rDerivativeShapeFunctionsLocal, const double rJacInv[9], std::vector<double>& rDerivativeShapeFunctionsGlobal)const;

    //! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
    //! eventually include also area/width of an element
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param ConstitutiveTangentBase constitutive tangent matrix
    //! @param rFactor factor including determinant of Jacobian and IP weight
    //! @param rRow row, where to start to add the submatrix
    //! @param rCol columns, where to start to add the submatrix
    //! @param rCoefficientMatrix to be added to
    void AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                                  const ConstitutiveTangentLocal<6,6>& rConstitutiveTangent, double rFactor,
                                  int rRow, int rCol,
                                  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const;

    //! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
    //! eventually include also area/width of an element (that's the thermal solution)
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param ConstitutiveTangentBase constitutive tangent matrix
    //! @param rFactor factor including determinant of Jacobian and IP weight
    //! @param rRow row, where to start to add the submatrix
    //! @param rCoefficientMatrix to be added to
    void AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                                  const ConstitutiveTangentLocal<3,3>& rConstitutiveTangent, double rFactor,
                                  int rRow, int rCol,
                                  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const;

    //! @brief adds up the internal force vector
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rEngineeringStress stress
    //! @param factor factor including det Jacobian area and integration point weight
    //! @param rRow start row (in case of a multifield problem)
    //! @param rResult resforce vector
    void AddDetJBtSigma(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                                     const EngineeringStress3D& rEngineeringStress,
                                     double rFactor,
                                     int rRow,
                                     FullVector<double,Eigen::Dynamic>& rResult)const;

    //! @brief adds up the internal force vector
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rHeatFlux stress
    //! @param factor factor including det Jacobian area and integration point weight
    //! @param rRow start row (in case of a multifield problem)
    //! @param rResult resforce vector
    void AddDetJBtHeatFlux(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                                     const HeatFlux3D& rHeatFlux,
                                     double rFactor,
                                     int rRow,
                                     FullVector<double,Eigen::Dynamic>& rResult)const;

    //! @brief calculates the integration point data with the current displacements applied
    //! @param rIpDataType data type to be stored for each integration point
    //! @param rIpData return value with dimension (dim of data type) x (numIp)
    //Error::eError GetIpData(NuTo::IpData::eIpStaticDataType rIpDataType, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rIpData)const;

    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const;

    //! @brief ... interpolate three-dimensional global point coordinates from three-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... three-dimensional local point coordinates
    //! @param rGlobalCoordinates ... three-dimension global point coordinates
    void InterpolateCoordinatesFrom3D(double rLocalCoordinates[3], double rGlobalCoordinates[3]) const;

    //! @brief ... interpolate three-dimensional global point displacements from three-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... three-dimensional local point coordinates
    //! @param rGlobalDisplacements ... three-dimension global point displacements
    void InterpolateDisplacementsFrom3D(double rLocalCoordinates[3], double rGlobalDisplacements[3]) const;

    //! @brief ... interpolate three-dimensional global temperature from three-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... three-dimensional local point coordinates
    //! @param rGlobalDisplacements ... three-dimension global point displacements
    void InterpolateTemperatureFrom3D(double rLocalCoordinates[3], double& rTemperature) const;

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @param rVolume  vector for storage of the ip volumes (area in 2D)
    void GetIntegrationPointVolume(std::vector<double>& rVolume)const;

    //! @brief cast the base pointer to an ElementSolid, otherwise throws an exception
    const Solid* AsSolid()const;

    //! @brief cast the base pointer to an ElementSolid, otherwise throws an exception
    Solid* AsSolid();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    Solid(){}

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element volume is negative)
    void CheckElement();

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    //! @param rNumDisp ... number of displacement dofs
    //! @param rNumTemp ... number of temperature dofs
    void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs, int rNumDisp, int rNumTemp) const;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    //! @param rNumDisp ... number of displacement dofs
    //! @param rNumTemp ... number of temperature dofs
    void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs, int rNumDisp, int rNumTemp) const;

    const SectionBase *mSection;
};

} // namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Solid)
#endif // ENABLE_SERIALIZATION

#endif //SOLID_H
