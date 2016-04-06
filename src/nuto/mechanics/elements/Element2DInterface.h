//============================================================================
// Name        : Element2DInterface.cpp
// Author      : Philip Huschke
// Version     : 26 Aug 2015
// Copyright   :
// Description : Element formulation for the interface element proposed by Goodman et al.
//============================================================================

#pragma once

#include "nuto/mechanics/elements/Element2D.h"

namespace NuTo
{


class InterfaceSlip;
class Element2DInterface: public Element2D
{

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    Element2DInterface(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType);

    virtual ~Element2DInterface()
    {
    }
    ;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Element2DInterface" << std::endl;
#endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Element2DInterface" << std::endl;
#endif
    }
#endif // ENABLE_SERIALIZATION

    //! @brief calculates output data for the element
    //! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    //! @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
    Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput) override;

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType() const override;

    //! @brief returns the local dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! @return local dimension
    int GetLocalDimension() const override;

    //! @brief calculates the interface slip of the element, i.e. the relative displacement of the top and bottom part
    const InterfaceSlip CalculateInterfaceSlip(const Eigen::VectorXd& rShapeFunctions, const Eigen::MatrixXd& rNodeDisplacements);

    //! @brief calculates the element stiffness matrix at one integration point
    void AddElementStiffnessMatrix(const ConstitutiveTangentLocal<Eigen::Dynamic, Eigen::Dynamic>& rConstitutiveMatrix, const Eigen::VectorXd& rShapefunctions, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoefficientMatrix, const double rGaussIntegrationFactor);

    //! @brief calculates the jacobian
    double CalculateDetJacobian(const Eigen::MatrixXd& rNodeCoordinates) const;

    //! @brief calculates the internal force vector
    void AddInternalForceVector(const ConstitutiveTangentLocal<Eigen::Dynamic, 1>& rInterfaceStresses, const Eigen::VectorXd& rShapefunctions, NuTo::FullVector<double, Eigen::Dynamic>& rInternalForceVector, const double rGaussIntegrationFactor);

    //! @brief calculates the roation matirx based on the orientation of the element
    Eigen::MatrixXd CalculateRotationMatrix();

    //! @brief calculates the transformation matrix, i.e. blows up the rotation matrix
    Eigen::MatrixXd CalculateTransformationMatrix(unsigned int rGlobalDimension, unsigned int rNumberOfNodes);

//    virtual const Eigen::MatrixXd ExtractNodeValues(int rTimeDerivative, Node::eAttributes) const override;
    using Element2D::ExtractNodeValues;
    void ExtractNodeValues(Eigen::MatrixXd& rNodeValues, int rTimeDerivative, Node::eAttributes rDofType) const override;

    const Eigen::VectorXi CalculateGlobalRowDofs() const;

#ifdef ENABLE_VISUALIZE
    virtual void Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList) override;

    void GetVisualizationCells(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates, unsigned int& NumVisualizationCells, std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType, std::vector<unsigned int>& VisualizationCellsIncidence,
            std::vector<unsigned int>& VisualizationCellsIP) const;
#endif // ENABLE_VISUALIZE

protected:
    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    void CheckElement() override;

    Element2DInterface() = default;

};

} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Element2DInterface)
#endif
