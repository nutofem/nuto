#pragma once

#include "nuto/mechanics/elements/ContinuumElement.h"
#include "nuto/mechanics/nodes/NodeEnum.h"


namespace NuTo
{

template <int TDim>
class ContinuumElementIGA: public ContinuumElement<TDim>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
protected:
    ContinuumElementIGA()=default;
#endif // ENABLE_SERIALIZATION

    friend class ContinuumBoundaryElement<TDim>;

public:
    ContinuumElementIGA(const NuTo::StructureBase*          rStructure,
                        const std::vector<NuTo::NodeBase*> &rNodes,
                        const Eigen::MatrixXd              &rKnots,
                        const Eigen::VectorXi              &rKnotIDs,
                        ElementData::eElementDataType       rElementDataType,
                        IpData::eIpDataType                 rIpDataType,
                        InterpolationType                  *rInterpolationType);

    ContinuumElementIGA(const ContinuumElementIGA& ) = default;
    ContinuumElementIGA(      ContinuumElementIGA&&) = default;

    virtual ~ContinuumElementIGA() = default;

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType() const override;

    //! @brief Calculates the the jacobian of the mapping between the refernce element and parametric space (knots)
    //! @param rKnots ... knots of the element
    //! @return the jacobian
    Eigen::Matrix<double, TDim, TDim> CalculateJacobianParametricSpaceIGA() const override;

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @return rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    const Eigen::VectorXd GetIntegrationPointVolume() const override;

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
    const Eigen::MatrixXd& GetKnots() const override  {return mKnots;}

    //! @brief returns the knots of the element
    //! @return reference on the matrix containing the knots
    const Eigen::VectorXi& GetKnotIDs() const override  {return mKnotIDs;}

    Eigen::VectorXd InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const;

    Eigen::VectorXd InterpolateDofGlobalSurfaceDerivative(int rTimeDerivative, int rSurfaceId, const Eigen::VectorXd& rParameter, int rDerivative, int rDirection) const override;

    const ContinuumElementIGA<1>& AsContinuumElementIGA1D() const override
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] Element is not of type ContinuumElementIGA<1>.");}

    const ContinuumElementIGA<2>& AsContinuumElementIGA2D() const override
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] Element is not of type ContinuumElementIGA<2>.");}

    const ContinuumElementIGA<3>& AsContinuumElementIGA3D() const override
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] Element is not of type ContinuumElementIGA<3>.");}

    ContinuumElementIGA<1>& AsContinuumElementIGA1D() override
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] Element is not of type ContinuumElementIGA<1>.");}

    ContinuumElementIGA<2>& AsContinuumElementIGA2D() override
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] Element is not of type ContinuumElementIGA<2>.");}

    ContinuumElementIGA<3>& AsContinuumElementIGA3D() override
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] Element is not of type ContinuumElementIGA<3>.");}


protected:

    //! @brief ... knots (e.g. 2D \f$\begin{pmatrix}\xi_i & \xi_{i+1} \\ \mu_j & \mu_{j+1}\f$)
    Eigen::MatrixXd mKnots;

    //! @brief ... knot ids
    Eigen::VectorXi mKnotIDs;

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    void CheckElement() override;

    virtual void CalculateNMatrixBMatrixDetJacobian(EvaluateDataContinuum<TDim> &data, int rTheIP) const;

    virtual double CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const;


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

