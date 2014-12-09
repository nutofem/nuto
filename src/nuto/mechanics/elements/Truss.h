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
class EngineeringStrain1D;
class HeatFlux1D;
class Damage;
class LocalEqPlasticStrain;
class LocalEqStrain;
class NonlocalEqPlasticStrain;
class NonlocalEqStrain;
class RelativeHumidity;
class WaterPhaseFraction;
template <int TNumRows, int TNumColumns> class ConstitutiveTangentLocal;

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

    //! @brief calculates output data fo the elmement
    //! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    //!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
    NuTo::Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput);

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

    //! @brief Update the static data of an element
    //Error::eError UpdateStaticData(NuTo::Element::eUpdateType rUpdateType);

    //! @brief calculates the Kkk matrix
    //! @param shapeFunctions of the ip for all shape functions
    //! @param derivativeShapeFunctions of the ip for all shape functions
    //! @param c nonlocal gradient radius
    //! @param factor multiplication factor (detJ area..)
    //! @param Kkk return matrix with detJ * NtT+cBtB
    void CalculateKkk(const std::vector<double>& shapeFunctions,const std::vector<double>& derivativeShapeFunctions,double nonlocalGradientRadius,double factor,
    		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& Kkk);

    //! @brief add Kkk*kappa+detJ*F (detJ is already included in Kkk)
    //! @param derivativeShapeFunctions of the ip for all shape functions
    //! @param tangentStressNonlocalEqPlasticStrain derivative of the stress with respect to the nonlocal eq plastic strain
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJBtdSigmadNonlocalEqPlasticStrainN(const std::vector<double>& rDerivativeShapeFunctions, ConstitutiveTangentLocal<1,2>& rTangentStressNonlocalEqPlasticStrain,
    		const std::vector<double>& rShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult);

    //! @brief add Kkk*kappa+detJ*F (detJ is already included in Kkk)
    //! @param derivativeShapeFunctions of the ip for all shape functions
    //! @param tangentStressNonlocalEqPlasticStrain derivative of the stress with respect to the nonlocal eq plastic strain
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJBtdSigmadNonlocalTotalStrainN(const std::vector<double>& rDerivativeShapeFunctions, ConstitutiveTangentLocal<1,1>& rTangentStressNonlocalTotalStrain,
    		const std::vector<double>& rShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult);

    //! @brief add detJ B.T dSigma/dnonlocalEqStrain N
    //! @param derivativeShapeFunctions of the ip for all shape functions
    //! @param tangentStressNonlocalEqStrain derivative of the stress with respect to the nonlocal eq strain
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJBtdSigmadNonlocalEqStrainN(const std::vector<double>& rDerivativeShapeFunctions, ConstitutiveTangentLocal<1,1>& rTangentStressNonlocalEqStrain,
            const std::vector<double>& rShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult);


    //! @brief add detJ transpose N dOmega/depsilon B
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rTangentLocalEqPlasticStrainStrain derivative of the local eq plastic strains with respect to the strain
    //! @param rderivativeShapeFunctions of the ip for all shape functions
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJNtdLocalEqPlasticStraindEpsilonB(const std::vector<double>& rShapeFunctions, ConstitutiveTangentLocal<2,1>& rTangentLocalEqPlasticStrainStrain,
    		const std::vector<double>& rDerivativeShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult);

    //! @brief add detJ transpose N dOmega/depsilon B
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rTangentLocalEqStrainStrain derivative of the local eq plastic strains with respect to the strain
    //! @param rderivativeShapeFunctions of the ip for all shape functions
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJNtdLocalEqStraindEpsilonB(const std::vector<double>& rShapeFunctions, ConstitutiveTangentLocal<1,1>& rTangentLocalEqStrainStrain,
            const std::vector<double>& rDerivativeShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult);


    //! @brief add detJ transpose N dOmega/depsilon B
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rTangentLocalEqPlasticStrainStrain derivative of the local eq plastic strains with respect to the strain
    //! @param rderivativeShapeFunctions of the ip for all shape functions
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJNtB(const std::vector<double>& rShapeFunctions,
    		const std::vector<double>& rDerivativeShapeFunctions, double rFactor, int rRow, int rCol,
    		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult);

    //! @brief adds to a matrix the product N^t C N, where N contains the the shape functions and C is the constitutive tangent
    //! eventually include also area/width of an element
    //! @param rShapeFunctions shape functions
    //! @param ConstitutiveTangentBase constitutive tangent matrix
    //! @param rFactor factor including area, determinant of Jacobian and IP weight
    //! @param rRow row, where to start to add the submatrix
    //! @param rCol col, where to start to add the submatrix
    //! @param rCoefficientMatrix to be added to
    void AddDetJNtCN(const std::vector<double>& rShapeFunctions,
                     const ConstitutiveTangentLocal<1,1>& rConstitutiveTangent, double rFactor,
                     int rRow, int rCol,
                     FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const;

    //! @brief calculates the local coordinates of the nodes
    //! @param localCoordinates vector with already correct size allocated
    //! this can be checked with an assertation
    virtual void CalculateLocalCoordinates(std::vector<double>& rLocalCoordinates)const=0;

    //! @brief calculates the local displacements of the nodes
    //! @param time derivative (0 displacements, 1 velocities, 2 accelerations)
    //! @param localDisplacements vector with already correct size allocated
    //! this can be checked with an assertation
    virtual void CalculateLocalDisplacements(int rTimeDerivative, std::vector<double>& rLocalDisplacements)const=0;

    //! @brief stores the temperatures of the nodes
    //! @param time derivative (0 temperature, 1 temperature rate, 2 second time derivative of temperature)
    //! @param temperature vector with already correct size allocated
    //! this can be checked with an assertation
    void CalculateNodalTemperatures(int rTimeDerivative, std::vector<double>& rTemperatures)const;

    //! @brief stores the nonlocal eq plastic strain of the nodes
    //! @param time derivative (0 damage, 1 damage rate, 2 second time derivative of damage)
    //! @param nonlocal eq plastic strain vector with already correct size allocated (2*nodes)
    //! this can be checked with an assertation
    void CalculateNodalNonlocalEqPlasticStrain(int rTimeDerivative, std::vector<double>& rNodalNonlocalEquivalentPlasticStrain)const;

    //! @brief returns the nonlocal eq plastic strain interpolated from the nodal values
    //! @param shapeFunctionsGlobal shape functions
    //! @param rNodeDamage nonlocal eq plastic strain values of the nodes
    //! @param rNonlocalEqentPlasticStrain return value
    void CalculateNonlocalEqPlasticStrain(const std::vector<double>& shapeFunctions,
    		const std::vector<double>& rNodeEquivalentPlasticStrain, NonlocalEqPlasticStrain& rNonlocalEqentPlasticStrain)const;

    //! @brief stores the nonlocal total strain of the nodes
    //! @param time derivative (0 damage, 1 damage rate, 2 second time derivative of damage)
    //! @param nonlocal total strain vector with already correct size allocated (1*nodes)
    //! this can be checked with an assertation
    void CalculateNodalNonlocalTotalStrain(int rTimeDerivative, std::vector<double>& rNodalNonlocalTotalStrain)const;

    //! @brief returns the nonlocal total strain interpolated from the nodal values
    //! @param shapeFunctionsGlobal shape functions
    //! @param rNodeNonlocalTotalStrain nonlocal total strain values of the nodes
    //! @param rNonlocalTotalStrain return value
    void CalculateNonlocalTotalStrain(const std::vector<double>& shapeFunctions,
    		const std::vector<double>& rNodeNonlocalTotalStrain, EngineeringStrain1D& rNonlocalTotalStrain)const;

    //! @brief stores the nonlocal eq strain of the nodes
    //! @param time derivative (0 damage, 1 damage rate, 2 second time derivative of damage)
    //! @param nonlocal eq strain vector with already correct size allocated (1*nodes)
    //! this can be checked with an assertation
    void CalculateNodalNonlocalEqStrain(int rTimeDerivative, std::vector<double>& rNodalNonlocalEqStrain)const;

    //! @brief returns the nonlocal eq strain interpolated from the nodal values
    //! @param shapeFunctions shape functions
    //! @param rNodeNonlocalEqStrain nonlocal eq strain values of the nodes
    //! @param rNonlocalEqStrain return value
    void CalculateNonlocalEqStrain(const std::vector<double>& shapeFunctions,
            const std::vector<double>& rNodeNonlocalEqStrain, NonlocalEqStrain& rNonlocalEqStrain)const;

    //! @brief stores the relative humidity of the nodes
    //! @param time derivative
    //! @param relative humidity vector with already correct size allocated (1*nodes)
    //! this can be checked with an assertation
    void CalculateNodalRelativeHumidity(int rTimeDerivative, std::vector<double>& rNodalRelativeHumidity)const;

    //! @brief returns the relative humidity interpolated from the nodal values
    //! @param shapeFunctions shape functions
    //! @param rNodeRelativeHumidity relative humidity values of the nodes
    //! @param rRelativeHumidity return value
    void CalculateRelativeHumidity(const std::vector<double>& rShapeFunctions,
                                   const std::vector<double>& rNodeRelativeHumidity,
                                   RelativeHumidity& rRelativeHumidity)const;

    //! @brief stores the water phase fraction of the nodes
    //! @param time derivative
    //! @param water phase fraction vector with already correct size allocated (1*nodes)
    //! this can be checked with an assertation
    void CalculateNodalWaterPhaseFraction(int rTimeDerivative, std::vector<double>& rNodalWaterPhaseFraction)const;

    //! @brief returns the water phase fraction interpolated from the nodal values
    //! @param shapeFunctions shape functions
    //! @param rNodeWaterPhaseFraction water phase fraction values of the nodes
    //! @param rWaterPhaseFraction return value
    void CalculateWaterPhaseFraction(const std::vector<double>& rShapeFunctions,
                                     const std::vector<double>& rNodeWaterPhaseFraction,
                                     WaterPhaseFraction& rWaterPhaseFraction)const;

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
    //! @param rRerivativeShapeFunctions derivatives of the shape functions (local coordinates)
    //! @param rLocalDisp local displacements
    //! @param rConstitutiveInput (return value)
    void CalculateDeformationGradient(const std::vector<double>& rDerivativeShapeFunctions,
                                      const std::vector<double>& rLocalDisp,
                                      DeformationGradient1D& rDeformationGradient)const;

    //! @brief returns the local dimension of the element
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    int GetLocalDimension()const
    {
        return 1;
    }

    //! @brief returns the number of shape functions for the interpolation of the nonlocal total strains
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    virtual int GetNumShapeFunctionsNonlocalTotalStrain()const;

    //! @brief returns the number of shape functions for the interpolation of the nonlocal eq strains
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    virtual int GetNumShapeFunctionsNonlocalEqStrain()const;

    //! @brief returns the number of shape functions for the interpolation of the moisture transport
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    virtual int GetNumShapeFunctionsMoistureTransport() const;

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
    virtual void CalculateShapeFunctionsGeometry(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const=0;

    //! @brief calculates the shape functions, uses the geometry shape functions unless implemented differently
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
    virtual void CalculateShapeFunctionsField(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const;

    //! @brief calculates the shape functions, uses the geometry shape functions unless implemented differently
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes
    virtual void CalculateShapeFunctionsNonlocalTotalStrain(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const;

    //! @brief calculates the shape functions, uses the geometry shape functions unless implemented differently
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes
    virtual void CalculateShapeFunctionsNonlocalEqStrain(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const;

    //! @brief calculates the shape functions, uses the geometry shape functions unless implemented differently
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes
    virtual void CalculateShapeFunctionsMoistureTransport(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsGeometry(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const=0;

    //! @brief calculates the derivative of the shape functions, uses the derivatives of the geometry shape function unless implemented differently
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsField(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions, uses the derivatives of the geometry shape function unless implemented differently
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsNonlocalTotalStrain(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions, uses the derivatives of the geometry shape function unless implemented differently
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsNonlocalEqStrain(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions, uses the derivatives of the geometry shape function unless implemented differently
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsMoistureTransport(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const;


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
    //! @param rRow row, where to start to add the submatrix
    //! @param rCol col, where to start to add the submatrix
    //! @param rCoefficientMatrix to be added to
    void AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctions,
                                  const ConstitutiveTangentLocal<1,1>& rConstitutiveTangent, double rFactor,
                                  int rRow, int rCol,
                                  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const ;

    //! @brief adds up the internal force vector
    //! @param derivativeShapeFunctions derivatives of the shape functions
    //! @param rEngineeringStress stress
    //! @param factor factor including det Jacobian area and integration point weight
    //! @param rRow start row (in case of a multifield problem)
    //! @param rResult resforce vector
    void AddDetJBtSigma(const std::vector<double>& rDerivativeShapeFunctions,
                                const EngineeringStress1D& rEngineeringStress, double factor,
                                int rRow,
                                FullVector<double,Eigen::Dynamic>& rResult)const;

    //! @brief adds up the internal force vector
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rHeatFlux stress
    //! @param factor factor including det Jacobian area and integration point weight
    //! @param rRow start row (in case of a multifield problem)
    //! @param rResult resforce vector
    void AddDetJBtHeatFlux(const std::vector<double>& rDerivativeShapeFunctions,
                                     const HeatFlux1D& rHeatFlux,
                                     double rFactor,
                                     int rRow,
                                     FullVector<double,Eigen::Dynamic>& rResult)const;

    //! @brief add Kkk*omega+detJ*F (detJ is already included in Kkk)
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rLocalEqPlasticStrain local eq. plastic strain values
    //! @param rKkk stiffness matrix Kkk
    //! @param rNodeNonlocalEqPlasticStrain nodal nonlocal eq plastic strain values
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJRnonlocalPlasticStrain(const std::vector<double>& rShapeFunctions,const LocalEqPlasticStrain& rLocalEqPlasticStrain, const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkk,
    		const std::vector<double>& rNodeNonlocalEqPlasticStrain, double rFactor, int rRow, FullVector<double,Eigen::Dynamic>& rResult) const;

    //! @brief add Kkk*omega+detJ*F (detJ is already included in Kkk)
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rLocalTotalStrain local total strain values
    //! @param rKkk stiffness matrix Kkk
    //! @param rNodeNonlocalTotalStrain nodal nonlocal total strain values
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJRnonlocalTotalStrain(const std::vector<double>& rShapeFunctions,const EngineeringStrain1D& rLocalTotalStrain,
    		const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkk,
    		const std::vector<double>& rNodeNonlocalTotalStrain, double rFactor, int rRow, FullVector<double,Eigen::Dynamic>& rResult) const;

    //! @brief add Kkk*nonlocalEqStrain-detJ*N.T*localEqStrain (detJ is already included in Kkk)
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rLocalEqStrain local eq. strain values
    //! @param rKkk stiffness matrix Kkk
    //! @param rNodeNonlocalEqStrain nodal nonlocal eq strain values
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJRnonlocalEqStrain(const std::vector<double>& rShapeFunctions,const LocalEqStrain& rLocalEqStrain, const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkk,
            const std::vector<double>& rNodeNonlocalEqStrain, double rFactor, int rRow, FullVector<double,Eigen::Dynamic>& rResult) const;

    
    //! @brief transforms the local matrix to the global system
    //! relevant only for 2D and 3D truss elements
    virtual void BlowLocalMatrixToGlobal(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rFullCoefficientMatrix)const=0;

    //! @brief transforms the local vector to the global system
    //! relevant only for 2D and 3D truss elements
    virtual void BlowLocalVectorToGlobal(NuTo::FullVector<double,Eigen::Dynamic>& rFullVector)const=0;

    //! @brief calculates the integration point data with the current displacements applied
    //! @param rIpDataType data type to be stored for each integration point
    //! @param rIpData return value with dimension (dim of data type) x (numIp)
    //Error::eError GetIpData(NuTo::IpData::eIpStaticDataType rIpDataType, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rIpData)const;

    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const;

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @param rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    void GetIntegrationPointVolume(std::vector<double>& rVolume)const;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual NodeBase* GetNodeNonlocalTotalStrain(int rLocalNodeNumber);

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual const NodeBase* GetNodeNonlocalTotalStrain(int rLocalNodeNumber)const;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual NodeBase* GetNodeNonlocalEqStrain(int rLocalNodeNumber);

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual const NodeBase* GetNodeNonlocalEqStrain(int rLocalNodeNumber)const;

    //! @brief cast the base pointer to an ElementTruss, otherwise throws an exception
    const Truss* AsTruss()const;

    //! @brief cast the base pointer to an ElementTruss, otherwise throws an exception
    Truss* AsTruss();

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

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    //! @param rNumXxxxDofs ... number of Xxxx dofs
    virtual void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs,int rNumDispDofs, int rNumTempDofs, int rNumNonlocalEqPlasticStrainDofs,int rNumNonlocalTotalStrainDofs, int rNumNonlocalEqStrainDofs, int rNumRelativeHumidityDofs, int rNumWaterPhaseFractionDofs) const=0;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    //! @param rNumXxxxDofs ... number of Xxxx dofs
    virtual void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs,int rNumDispDofs, int rNumTempDofs, int rNumNonlocalEqPlasticStrainDofs,int rNumNonlocalTotalStrainDofs, int rNumNonlocalEqStrainDofs, int rNumRelativeHumidityDofs, int rNumWaterPhaseFractionDofs) const;

    //! @brief adds to a matrix the product factor * H^tH, where H contains the shape functions
    //! @param rShapeFunctions ... shape functions
    //! @param rFactor factor including area, determinant of Jacobian, IP weight and density
    //! @param rCoefficientMatrix to be added to
    virtual void AddDetJHtH(const std::vector<double>& rShapeFunctions,
                            double rFactor,
                            FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const;
};

} // namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Truss)
#endif // ENABLE_SERIALIZATION

#endif //TRUSS_H
