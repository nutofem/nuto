#pragma once

#include "nuto/mechanics/elements/ContinuumElementIGA.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/base/ErrorEnum.h"

namespace NuTo
{
//! @author Peter Otto, BAM
//! @date July, 2016
//! @brief ... class for contimuum elements based on a nurbs discretization
//! @brief ... this elements are common iga elements, but also iga layers

template <int TDim>
class ContinuumElementIGALayer: public ContinuumElementIGA<TDim>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
protected:
    ContinuumElementIGALayer()=default;
#endif // ENABLE_SERIALIZATION

    friend class ContinuumBoundaryElement<TDim>;

public:
    ContinuumElementIGALayer(const NuTo::StructureBase*          rStructure,
                             const std::vector<NuTo::NodeBase*> &rNodes,
                             const Eigen::MatrixXd              &rKnots,
                             const Eigen::VectorXi              &rKnotIDs,
                             ElementData::eElementDataType       rElementDataType,
                             IpData::eIpDataType                 rIpDataType,
                             InterpolationType                  *rInterpolationType);

    ContinuumElementIGALayer(const ContinuumElementIGALayer& ) = default;
    ContinuumElementIGALayer(      ContinuumElementIGALayer&&) = default;

    virtual ~ContinuumElementIGALayer() = default;

    //! @brief calculates output data for the element
    //! @param rInput ... constitutive input map for the constitutive law
    //! @param rOutput ...  coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    eError Evaluate(const ConstitutiveInputMap& rInput, std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rOutput) override;
//    {
//        (void)rInput;
//        (void)rOutput;
//        return eError::NOT_IMPLEMENTED;
//    }

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType() const override;

    //! @brief Calculates the the jacobian of the mapping between the refernce element and parametric space (knots)
    //! @param rKnots ... knots of the element
    //! @return the jacobian
    Eigen::Matrix<double, TDim, TDim> CalculateJacobianParametricSpaceIGA() const override;

    Eigen::MatrixXd CalculateJacobianSurface(const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::VectorXd& rNodeCoordinates) const;

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @return rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    const Eigen::VectorXd GetIntegrationPointVolume() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"Probably not needed.");
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NuTo::NodeBase* GetNode(int rLocalNodeNumber) override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NuTo::NodeBase* GetNode(int rLocalNodeNumber) const override;

    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode) override;

    //! @brief returns the knots of the element
    //! @return reference on the matrix containing the knots
    const Eigen::MatrixXd& GetKnots() const override  {return this->mKnots;}

    //! @brief returns the knots of the element
    //! @return reference on the matrix containing the knots
    const Eigen::VectorXi& GetKnotIDs() const override  {return this->mKnotIDs;}

    Eigen::VectorXd InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const;

    //! @brief returns the position of at rNaturalCoordinates, the interpolation of init and current dofs should be the same
    //! @return reference on the matrix containing the knots
    Eigen::VectorXd InterpolateDofGlobalCurrentConfiguration(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofTypeInit, Node::eDof rDofTypeCurrent) const override;

    Eigen::VectorXd InterpolateDofGlobalSurfaceDerivative(int rTimeDerivative, const Eigen::VectorXd& rParameter, int rDerivative, int rDirection) const override;

    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> InterpolateDofGlobalSurfaceDerivativeTotal(int rTimeDerivative, const Eigen::VectorXd& rParameter, int rDerivative) const override;

//    Eigen::VectorXd InterpolateDofGlobalSurfaceNormal(const Eigen::VectorXd& rParameter) const override;


protected:

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    void CheckElement() override
    {
        // Probably not needed
    }

    virtual void CalculateNMatrixBMatrixDetJacobian(EvaluateDataContinuum<TDim> &data, int rTheIP) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"Probably not needed.");
    }

    virtual double CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"Probably not needed.");
    }


#ifdef ENABLE_SERIALIZATION
private:
    //! @brief serializes the class, this is the load routine
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version);

    //! @brief serializes the class, this is the save routine
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;

    BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif  // ENABLE_SERIALIZATION
};

} /* namespace NuTo */

