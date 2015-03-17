// $Id$
#ifndef PLANE_H
#define PLANE_H

#include "nuto/mechanics/elements/ElementBase.h"

namespace NuTo
{
class StructureBase;
class DeformationGradient2D;
class TemperatureGradient2D;
class LocalEqStrain;
class NonlocalEqStrain;
class ConstitutiveTangentLocal3x3;
class EngineeringStress2D;
class HeatFlux2D;
template <int TNumRows, int TNumColumns> class ConstitutiveTangentLocal;
template <int TNumRows, int TNumColumns> class ConstitutiveTangentNonlocal;
//! @author Jörg F. Unger, ISM
//! @date March 2010
//! @brief ... standard abstract class for all isoparametric displacement based finite elements in 2D
class Plane : public ElementBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    Plane(const StructureBase* rStructure,
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

    //! @brief calculates the local coordinates of the nodes
    //! @param localCoordinates vector with already correct size allocated
    //! this can be checked with an assertation
    virtual void CalculateLocalCoordinates(std::vector<double>& rLocalCoordinates)const=0;

    //! @brief calculates the local displacements of the nodes
    //! @param localDisplacements vector with already correct size allocated
    //! this can be checked with an assertation
    virtual void CalculateLocalDisplacements(std::vector<double>& rLocalDisplacements)const=0;

    //! @brief calculates the nonlocal eq strains of the nodes
    //! @param rNodalNonlocalEqStrains vector with already correct size allocated
    //! this can be checked with an assertation
    virtual void CalculateLocalNonlocalEqStrains(std::vector<double>& rNodeNonlocalEqStrains)const=0;

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

    //! @brief calculates the deformation gradient in 2D
    //! @param rRerivativeShapeFunctions derivatives of the shape functions with respect to local coordinates (plane world coordinates)
    //! @param rLocalDisp local displacements
    //! @param rDeformationGradient (return value)
    void CalculateDeformationGradient(const std::vector<double>& rDerivativeShapeFunctionsLocal,
                                      const std::vector<double>& rDisp,
                                      DeformationGradient2D& rDeformationGradient)const;

    //! @brief stores the temperatures of the nodes
    //! @param temperature vector with already correct size allocated
    //! this can be checked with an assertation
    void CalculateTemperatures(std::vector<double>& rTemperatures)const;

    //! @brief returns the nonlocal eq strain interpolated from the nodal values
    //! @param shapeFunctions shape functions
    //! @param rNodeNonlocalEqStrain nonlocal eq strain values of the nodes
    //! @param rNonlocalEqStrain return value
    void CalculateNonlocalEqStrain(const std::vector<double>& shapeFunctions,
            const std::vector<double>& rNodeNonlocalEqStrain, NonlocalEqStrain& rNonlocalEqStrain)const;

    //! @brief calculates the temperature gradient in 3D
    //! @param rRerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rTemp nodal temperatures
    //! @param rTemperatureGradient (return value)
    void CalculateTemperatureGradient(const std::vector<double>& rDerivativeShapeFunctionsLocal,
            const std::vector<double>& rTemp,
            TemperatureGradient2D& rTemperatureGradient)const;

    //! @brief Calculates the the inverse of the Jacobian and its determinant
    //! @param rDerivativeShapeFunctions Derivatives of the shape functions (dN1dx, dN1dy, dN1dz, dN2dx, ..
    //! @param rNodeCoordinates Node coordinates (X1,Y1,Z1,X2,Y2,Z2,...
    //! @param rInvJacobian inverse Jacobian matrix (return value)
    //! @param rDetJac determinant of the Jacobian (return value)
    void CalculateJacobian(const std::vector<double>& rDerivativeShapeFunctions,
                           const std::vector<double>& rNodeCoordinates,
                           double rInvJacobian[4],
                           double& rDetJac)const;

    //! @brief returns the coordinates of an integration point
    //! @param rIpNum integration point
    //! @param rCoordinates coordinates to be returned
    virtual void GetLocalIntegrationPointCoordinates(int rIpNum, double rCoordinates[2])const;

    //! @brief returns the coordinates of an integration point
    //! @param rIpNum integration point
    //! @param rCoordinates coordinates to be returned
    virtual void GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const;

    //! @brief calculates the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
    virtual void CalculateShapeFunctionsGeometry(const double rLocalCoordinates[2], std::vector<double>& rShapeFunctions)const=0;

    //! @brief calculates the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
    virtual void CalculateShapeFunctionsField(const double rLocalCoordinates[2], std::vector<double>& rShapeFunctions)const;

    //! @brief calculates the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for the nonlocal eq strain nodes, size should already be correct, but can be checked with an assert
    virtual void CalculateShapeFunctionsNonlocalEqStrain(const double rLocalCoordinates[2], std::vector<double>& rShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions with respect to local coordinates
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsGeometryNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const=0;

    //! @brief calculates the derivative of the shape functions with respect to local coordinates
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsFieldNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions with respect to local coordinates
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsNonlocalEqStrainNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const;

    //! @brief calculates the derivative of the shape functions with respect to local coordinates (Jacobi transformation)
    //! @param std::vector<double>& rDerivativeShapeFunctions derivatives of the shape functions
    //! @param rJacInv inverse of the Jacobian
    //! @param rDerivativeShapeFunctionsGlobal derivaties of the shape functions with respect to global coordinates
    //! size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsLocal(const std::vector<double>& rDerivativeShapeFunctionsNatural, const double rJacInv[4], std::vector<double>& rDerivativeShapeFunctionsLocal)const;

    //! @brief returns the number of nodes in this element(field interpolation)
    //! @return number of nodes
    virtual int GetNumNodesField()const;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual NodeBase* GetNodeField(int rLocalNodeNumber);

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual const NodeBase* GetNodeField(int rLocalNodeNumber)const;

    //! @brief returns the number of nodes in this element(nonlocal eq strain interpolation)
    //! @return number of nodes
    virtual int GetNumNodesNonlocalEqStrain()const;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual NodeBase* GetNodeNonlocalEqStrain(int rLocalNodeNumber);

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual const NodeBase* GetNodeNonlocalEqStrain(int rLocalNodeNumber)const;

    //! @brief adds to a matrix the product B^tCBnonlocal, where B contains the derivatives of the shape functions and C is the constitutive tangent and Bnonlocal is the nonlocal B matrix
    //! eventually include also area/width of an element
    //! @param rLocalDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rNonlocalDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param ConstitutiveTangentBase constitutive tangent matrix
    //! @param rFactor factor including determinant of Jacobian and IP weight
    //! @param rCoefficientMatrix to be added to
    //! &param rFirstCol first column of the coefficient matrix to be modified (corresponding to the current nonlocal element)
    void AddDetJBtCB(const std::vector<double>& rLocalDerivativeShapeFunctionsLocal,const std::vector<double>& rNonlocalDerivativeShapeFunctionsLocal,
                                  const ConstitutiveTangentLocal<3,3>& rConstitutiveTangent, double rFactor,
                                  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix, int rFirstRow, int rFirstCol)const;

    //! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
    //! eventually include also area/width of an element (that's the thermal solution)
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param ConstitutiveTangentBase constitutive tangent matrix
    //! @param rFactor factor including determinant of Jacobian and IP weight
    //! @param rRow row, where to start to add the submatrix
    //! @param rCoefficientMatrix to be added to
    void AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                                  const ConstitutiveTangentLocal<2,2>& rConstitutiveTangent, double rFactor,
                                  int rRow, int rCol,
                                  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const;

    //! @brief adds up the internal force vector
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rEngineeringStress stress
    //! @param factor factor including det Jacobian area and integration point weight
    //! @param rRow start row (in case of a multifield problem)
    //! @param rResult resforce vector
    void AddDetJBtSigma(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                        const EngineeringStress2D& rEngineeringStress,
                        double factor,
                        int rRow,
                        FullVector<double,Eigen::Dynamic>& rResult)const;

    //! @brief adds up the internal force vector
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rHeatFlux stress
    //! @param factor factor including det Jacobian area and integration point weight
    //! @param rRow start row (in case of a multifield problem)
    //! @param rResult resforce vector
    void AddDetJBtHeatFlux(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                                     const HeatFlux2D& rHeatFlux,
                                     double rFactor,
                                     int rRow,
                                     FullVector<double,Eigen::Dynamic>& rResult)const;

    //! @brief calculates the Kkk matrix
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rDerivativeShapeFunctions of the ip for all shape functions
    //! @param nonlocal gradient radius xi
    //! @param factor multiplication factor (detJ area..)
    //! @param Kkk return matrix with detJ * (Nt 1/ct N + BtB)
    void CalculateKkkXi(const std::vector<double>& rShapeFunctions,const std::vector<double>& rDerivativeShapeFunctions,double rNonlocalParameterXi,double factor,
            FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& Kkk);

    //! @brief add Kkk*nonlocalEqStrain-detJ*N.T*localEqStrain (detJ is already included in Kkk)
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rLocalEqStrain local eq. strain values
    //! @param rKkk stiffness matrix Kkk
    //! @param rNodeNonlocalEqStrain nodal nonlocal eq strain values
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJRnonlocalEqStrain(const std::vector<double>& rShapeFunctions,const LocalEqStrain& rLocalEqStrain, const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkk,
            const std::vector<double>& rNodeNonlocalEqStrain, double rFactor, int rRow, FullVector<double,Eigen::Dynamic>& rResult) const;

    //! @brief add detJ B.T dSigma/dnonlocalEqStrain N
    //! @param derivativeShapeFunctions of the ip for all shape functions
    //! @param tangentStressNonlocalEqStrain derivative of the stress with respect to the nonlocal eq strain
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJBtdSigmadNonlocalEqStrainN(const std::vector<double>& rDerivativeShapeFunctions, ConstitutiveTangentLocal<3,1>& rTangentStressNonlocalEqStrain,
            const std::vector<double>& rShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult) const;

    //! @brief add detJ transpose N dOmega/depsilon B
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rTangentLocalEqStrain derivative of the local eq strains with respect to the strain
    //! @param rderivativeShapeFunctions of the ip for all shape functions
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJNtdLocalEqStraindEpsilonB(const std::vector<double>& rShapeFunctions, ConstitutiveTangentLocal<3,1>& rTangentLocalEqStrainStrain,
            const std::vector<double>& rDerivativeShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult) const;

    //! @brief calculates the integration point data with the current displacements applied
    //! @param rIpDataType data type to be stored for each integration point
    //! @param rIpData return value with dimension (dim of data type) x (numIp)
    //Error::eError GetIpData(NuTo::IpData::eIpStaticDataType rIpDataType, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rIpData)const;

    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const;

    //! @brief ... interpolate three-dimensional global point coordinates from two-dimensional local point coordinates (element coordinates system)
    //! @param rNaturalCoordinates ... two-dimensional point coordinates in natural coordinate system
    //! @param rGlobalCoordinates ... three-dimension global point coordinates
    void InterpolateCoordinatesFrom2D(double rNaturalCoordinates[2], double rGlobalCoordinates[2]) const override;

    //! @brief ... interpolate three-dimensional global point displacements from three-dimensional local point coordinates (element coordinates system)
    //! @param rTimeDerivative ... time derivative (0 disp, 1 velocities, 2 accelerations)
    //! @param rNaturalCoordinates ... two-dimensional point coordinates in natural coordinate system
    //! @param rGlobalDisplacements ... three-dimension global point displacements
    void InterpolateDisplacementsFrom2D(int rTimeDerivative, double rNaturalCoordinates[2], double rGlobalDisplacements[3]) const override;

    //! @brief ... interpolate three-dimensional global nonlocal eq strain from two-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... two-dimensional local point coordinates
    //! @param rNonlocalEqStrain ... interpolated nonlocal eq strain
    void InterpolateNonlocalEqStrainFrom2D(double rLocalCoordinates[2], double& rNonlocalEqStrain) const override;

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @param rVolume  vector for storage of the ip volumes (area in 2D)
    void GetIntegrationPointVolume(std::vector<double>& rVolume)const;

    //! @brief calculates the area of a plane element via the nodes (probably faster than sum over integration points)
    //! @return Area
    virtual double CalculateArea()const=0;

    //! @brief checks if a node is inside a polygon
    //! @param rPoint (input) ... a pointer to a 2D tuple containing the coordinates
    //! @param rPoint (input) ... a pointer to a vector of 2D tuples representing the polyline
    virtual bool CheckPointInsidePolygon( const std::tuple<double,double> *rPoint, const std::vector<std::tuple<double,double> > * rPolygon)const=0;

    //! @brief cast the base pointer to an Plane, otherwise throws an exception
    const Plane* AsPlane()const;

    //! @brief cast the base pointer to an Plane, otherwise throws an exception
    Plane* AsPlane();

    //! @brief calculate the length of a crack passing through the center of gravity and intersecting two edges
    //! @parameter alpha... angle of the crack
    virtual double CalculateCrackLength2D(double rAlpha)const=0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
    //! @brief ... just for serialization
    Plane(){};

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element volume is negative)
    void CheckElement();

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    //! @param rNumDispDofs number of displacement dofs
    //! @param rNumTempDofs number of temperature dofs
    //! @param rNumNonlocalEqStrainDofs number of nonlocal eq strain dofs
    virtual void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs, int rNumDispDofs, int rNumTempDofs, int rN) const = 0;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    //! @param rNumDispDofs number of displacement dofs
    //! @param rNumTempDofs number of temperature dofs
    //! @param rNumNonlocalEqStrainDofs number of nonlocal eq strain dofs
    virtual void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs, int rNumDispDofs, int rNumTempDofs, int rNumNonlocalEqStrainDofs) const = 0;

    const SectionBase *mSection;
};

} // namespace NuTo

#endif //PLANE_H
