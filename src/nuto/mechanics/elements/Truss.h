// $Id$
#ifndef TRUSS_H
#define TRUSS_H

#include "nuto/mechanics/elements/ElementBase.h"

namespace NuTo
{
class StructureBase;
class DeformationGradient1D;
class ConstitutiveTangentLocal1x1;
class EngineeringStress1D;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all isoparametric displacement based truss finite elements in 1D
class Truss : public ElementBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    Truss(const StructureBase* rStructure,
    		ElementData::eElementDataType rElementDataType,
    		IntegrationType::eIntegrationType rIntegrationType,
    		IpData::eIpDataType rIpDataType
    		);

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the stiffness matrix
    //! @param rResult ... coefficient matrix
    //! @param rGlobalDofsRow ... row numbers in global system
    //! @param rGlobalDofsColumn ... column numbers in global system
    //! @param rSymmetry ... matrix is symmetric or not (in the symmetric case the full matrix is also stored
    virtual void CalculateCoefficientMatrix_0(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const;

    //! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the damping matrix
    virtual void CalculateCoefficientMatrix_1(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const;

    //! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the Mass matrix
    virtual void CalculateCoefficientMatrix_2(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const;

    //! @brief calculates the gradient of the internal potential
    //! for a mechanical problem, this corresponds to the internal force vector
    void CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
                                            std::vector<int>& rGlobalDofs)const;

    //! @brief Update the static data of an element
    void UpdateStaticData(NuTo::Element::eUpdateType rUpdateType);

    //! @brief calculates the local coordinates of the nodes
    //! @param localCoordinates vector with already correct size allocated
    //! this can be checked with an assertation
    virtual void CalculateLocalCoordinates(std::vector<double>& rLocalCoordinates)const=0;

    //! @brief calculates the local displacements of the nodes
    //! @param localDisplacements vector with already correct size allocated
    //! this can be checked with an assertation
    virtual void CalculateLocalDisplacements(std::vector<double>& rLocalDisplacements)const=0;

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

    //! @brief calculates the input of the constitutive relation (probably deformation gradient in 1D, 2D or 3D)
    //! @param rRerivativeShapeFunctions derivatives of the shape functions
    //! @param rLocalDisp local displacements
    //! @param rConstitutiveInput (return value)
    void CalculateDeformationGradient(const std::vector<double>& rDerivativeShapeFunctions,
                                      const std::vector<double>& rLocalCoord, const std::vector<double>& rLocalDisp,
                                      DeformationGradient1D& rDeformationGradient)const;

    //! @brief returns the number of local degrees of freedom
    //! @return number of local degrees of freedom
    virtual int GetNumLocalDofs()const=0;

    //! @brief returns the local dimension of the element
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    int GetLocalDimension()const
    {
        return 1;
    }

    //! @brief returns the number of shape functions
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    virtual int GetNumShapeFunctions()const=0;

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point
    //! @param rCoordinates coordinates to be returned
    virtual void GetLocalIntegrationPointCoordinates(int rIpNum, double& rCoordinates)const;

    //! @brief returns the global coordinates of an integration point
    //! @param rIpNum integration point
    //! @param rCoordinates coordinates to be returned
    virtual void GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const;

    //! @brief calculates the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
    virtual void CalculateShapeFunctions(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const=0;

    //! @brief calculates the derivative of the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctions(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const=0;

    //! @brief returns determinant of the Jacobian
    //! @param derivativeShapeFunctions derivatives of the shape functions
    //! @param localCoord local coordinates
    //! @return determinant of the Jacobian
    virtual double DetJacobian(const std::vector<double>& derivativeShapeFunctions,const std::vector<double>& localCoord)const;

    //! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
    //! eventually include also area/width of an element
    //! @param rDerivativeShapeFunctions derivatives of the shape functions
    //! @param ConstitutiveTangentBase constitutive tangent matrix
    //! @param rFactor factor including area, determinant of Jacobian and IP weight
    //! @param rCoefficientMatrix to be added to
    virtual void AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctions,
                             const ConstitutiveTangentLocal1x1& rConstitutiveTangent, double rFactor,
                             FullMatrix<double>& rCoefficientMatrix)const;

    //! @brief adds up the internal force vector
    //! @param derivativeShapeFunctions derivatives of the shape functions
    //! @param rEngineeringStress stress
    //! @param factor factor including det Jacobian area and integration point weight
    //! @param rResult resforce vector
    virtual void AddDetJBtSigma(const std::vector<double>& rDerivativeShapeFunctions,
                                const EngineeringStress1D& rEngineeringStress, double factor, FullMatrix<double>& rResult)const;


    //! @brief transforms the local matrix to the global system
    //! relevant only for 2D and 3D truss elements
    virtual void BlowLocalMatrixToGlobal(NuTo::FullMatrix<double>& rFullCoefficientMatrix)const=0;

    //! @brief transforms the local vector to the global system
    //! relevant only for 2D and 3D truss elements
    virtual void BlowLocalVectorToGlobal(NuTo::FullMatrix<double>& rFullVector)const=0;

    // calculate list of global dofs related to the entries in the element stiffness matrix
    // rGlobalDofsRow global dofs corresponding to the rows of the matrix
    // rGlobalDofsColumn global dofs corresponding to the columns of the matrix
    virtual void CalculateGlobalDofs(std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const=0;


    //! @brief calculates the integration point data with the current displacements applied
    //! @param rIpDataType data type to be stored for each integration point
    //! @param rIpData return value with dimension (dim of data type) x (numIp)
    void GetIpData(NuTo::IpData::eIpStaticDataType rIpDataType, FullMatrix<double>& rIpData)const;

    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const;

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @param rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    void GetIntegrationPointVolume(std::vector<double>& rVolume)const;

    //! @brief cast the base pointer to an ElementTruss, otherwise throws an exception
    const Truss* AsTruss()const;

    //! @brief cast the base pointer to an ElementTruss, otherwise throws an exception
    Truss* AsTruss();

    //! @brief sets the fine scale model (deserialization from a binary file)
    void SetFineScaleModel(int rIp, std::string rFileName);

    #ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
    //! @brief ... just for serialization
    Truss(){}

    const SectionBase *mSection;

    //! @brief adds to a matrix the product factor * H^tH, where H contains the shape functions
    //! @param rShapeFunctions ... shape functions
    //! @param rFactor factor including area, determinant of Jacobian, IP weight and density
    //! @param rCoefficientMatrix to be added to
    virtual void AddDetJHtH(const std::vector<double>& rShapeFunctions,
                            double rFactor,
                            FullMatrix<double>& rCoefficientMatrix)const;
};

} // namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Truss)
#endif // ENABLE_SERIALIZATION

#endif //TRUSS_H
