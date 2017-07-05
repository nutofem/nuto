
#pragma once

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/nodes/NodeEnum.h"

#include <eigen3/Eigen/Sparse>
#include <cstring>
#include <fstream>
#include <set>
#include <mechanics/MechanicsEnums.h>
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/nodes/NodeBase.h"

#include "mechanics/feti/ReverseMap.h"
namespace NuTo
{


/**     \author Philip Huschke
 *      \date September 2016
 *      \brief Class for FETI methods
 *
 *
 *
 *      The standard interface problem for FETI methods
 *
 *      \f[
 *      \begin{bmatrix}
 *      \sum_{s=1}^{N_s} \boldsymbol{B}_s \boldsymbol{K}_s^+ \boldsymbol{B}_s^T & \boldsymbol{B}_{1} \boldsymbol{R}_{1}
 * & \dots & \boldsymbol{B}_{{N_s}} \boldsymbol{R}_{{N_s}}
 *      \\
 *      \boldsymbol{R}_{1}^T \boldsymbol{B}_{1}^T  & \boldsymbol{0} & \dots & \boldsymbol{0}
 *      \\
 *      \vdots & \vdots & \ddots & \vdots
 *      \\
 *      \boldsymbol{R}_{{N_s}}^T \boldsymbol{B}_{{N_s}}^T & \boldsymbol{0} & \dots & \boldsymbol{0}
 *      \end{bmatrix}
 *      \begin{bmatrix}
 *      \boldsymbol{\lambda}
 *      \\
 *      \boldsymbol{\alpha}_{1}
 *      \\
 *      \vdots
 *      \\
 *      \boldsymbol{\alpha}_{{N_s}}
 *      \end{bmatrix}
 *      =
 *      \begin{bmatrix}
 *      \sum_{s=1}^{N_s} \boldsymbol{B}_s \boldsymbol{K}_s^+ \boldsymbol{f}_s
 *      \\
 *      \boldsymbol{R}_{{1}}^T \boldsymbol{f}_{{1}}
 *      \\
 *      \vdots
 *      \\
 *      \boldsymbol{R}_{{N_s}}^T \boldsymbol{f}_{{N_s}}
 *      \end{bmatrix}
 *      \label{eq_interface_problem}
 *      \f]
 *
 *
 *      \f$ \boldsymbol{K}_s \f$ stiffness matrix of subdomain s (possibly singular)
 *
 *
 *      \f$ \boldsymbol{B}_s \f$ StructureFETI::mConnectivityMatrix of subdomain s
 *
 *
 *      \f$ \boldsymbol{R}_s \f$ StructureFETI::mRigidBodyModes of subdomain s
 *
 *
 *      \f$ \boldsymbol{\lambda} \f$ Lagrange multiplier
 *
 *
 *      \f$ \boldsymbol{\alpha}_s \f$ Linear combination of rigid body modes of subdomain s
 *
 *
 *      \f$ N_s \f$ Number of subdomains StructureFETI::mNumProcesses
 *
 *
 *
 *
 */
class StructureFeti : public Structure
{
private:
    struct Node
    {
        Eigen::Vector3d mCoordinates;
        int mId;
    };

    struct Element
    {
        int mId;
        std::vector<int> mNodeIds;
    };

    struct Interface
    {
        std::map<int, int> mNodeIdsMap;
        int mValue;
    };


public:
    using NodeList = std::vector<Node>;
    using ElementList = std::vector<Element>;
    using InterfaceList = std::vector<Interface>;
    using Matrix = Eigen::MatrixXd;
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using SparseMatrixMap = std::map<NuTo::Node::eDof, SparseMatrix>;
    using VectorXd = Eigen::VectorXd;

    //! @brief Constructor
    //! @param rDimension   Structural dimension (1,2 or 3)
    StructureFeti(int rDimension);

    ///
    /// \brief AssembleConnectivityMatrix
    ///
    /// Assembles StructureFETI::mConnectivityMatrix
    ///
    /// First it loops over StructureFETI::mInterfaces and inserts
    ///
    /// +1 if StructureFETI::mRank is lower than its neighbour
    /// -1 if StructureFETI::mRank is greater than its neighbour
    ///
    /// this guarantees the continuity of the field variable across interfaces.
    ///
    /// Next it loops over StructureFETI::mBoundaryDofIds and inserts +1 to constraint the corresponding DOF
    ///
    void AssembleConnectivityMatrix();

    ///
    /// \brief ApplyVirtualConstraints
    /// \param nodeIdsBoundaries
    /// \param nodeIdsLoads
    void ApplyVirtualConstraints(const std::vector<int>& nodeIdsBoundaries, const std::vector<int>& nodeIdsLoads);


    ///
    /// \brief ImportMeshJson
    ///
    /// Json mesh file generated with <a href="https://github.com/dbeurle/GmshReader">GmshReader</a>
    ///
    /// \param rFileName
    /// \param interpolationTypeId
    ///
    void ImportMeshJson(std::string rFileName, const int interpolationTypeId);


    /// \brief CalculateRigidBodyModesTotalFETI
    void CalculateRigidBodyModesTotalFETI();

    ///
    /// \brief CheckRigidBodyModes
    /// \param tolerance
    ///
    /// Checks if the calculated rigid body modes are in the null space of \f$ \boldsymbol{K}_s \f$
    ///
    void CheckRigidBodyModes(const StructureOutputBlockMatrix& hessian0, const double tolerance = 1.e-9) const;

    ///
    /// \brief CheckProjectionMatrix
    /// \param tolerance
    ///
    /// Checks if the calculated projection matrix satisfies \f$ \boldsymbol{P} -\boldsymbol{P} \boldsymbol{P} =
    /// \boldsymbol{0} \f$
    ///
    void CheckProjectionMatrix(const double tolerance = 1.e-9) const;

    ///
    /// \brief CheckProjectionMatrix
    /// \param tolerance
    ///
    /// Checks if the calculated projection matrix satisfies \f$ \boldsymbol{P}^T -\boldsymbol{G} = \boldsymbol{0} \f$
    ///
    void CheckProjectionOfCoarseGrid(const double tolerance = 1.e-9) const;

    ///
    /// \brief CheckStiffnessPartitioning
    /// \param hessian0
    /// \param tolerance
    ///
    void CheckStiffnessPartitioning(const StructureOutputBlockMatrix hessian0, const double tolerance = 1.e-9) const;
    void CreateDummy1D();
    void CalculateInterfaceRigidBodyModes();


    ///
    /// \brief Constraints all degrees of freedom of nodes in node group
    /// \param nodeGroupId
    ///
    void ApplyConstraintsTotalFeti(const Group<NodeBase>& nodeGroup);

    ///
    /// \brief Constraints all degrees of freedom
    /// \param dofIds to be constrained
    ///
    void ApplyConstraintsTotalFeti(const std::vector<int>& dofIds);

    ///
    /// \param dofIdAndPrescribedDisplacementMap
    void ApplyPrescribedDisplacements(const std::map<int, double> dofIdAndPrescribedDisplacementMap);

    ///
    void CalculateProjectionMatrix();

    ///
    void CalculateG();
    ///
    std::pair<int, int> CreateRectangularMesh2D(const std::vector<double>& meshDimensions,
                                                const std::vector<int>& numElements)
    {

        assert(GetDimension() == static_cast<int>(meshDimensions.size()) and "Dimensions mismatch");
        assert(GetDimension() == static_cast<int>(numElements.size()) and "Dimensions mismatch");

        mNumInterfaceNodesTotal = (mNumProcesses - 1) * (numElements[1] + 1);

        double deltaX = meshDimensions[0] / mNumProcesses;

        std::vector<double> endPoints;
        endPoints.push_back((mRank + 1) * deltaX);
        endPoints.push_back(meshDimensions[1]);

        std::vector<double> startPoints;
        startPoints.push_back(mRank * deltaX);
        startPoints.push_back(0);

        // create local meshes
        auto importContainer = NuTo::MeshGenerator::Grid(*this, startPoints, endPoints, numElements,
                                                         NuTo::Interpolation::eShapeType::QUAD2D);


        int numInterfaces;
        int globalNodeId = mRank == 0 ? 0 : (mRank - 1) * (numElements[1] + 1);

        if (mRank == 0 or mRank == mNumProcesses - 1)
        {
            numInterfaces = 1;
        }
        else
        {
            numInterfaces = 2;
        }

        mInterfaces.resize(numInterfaces);
        int interfaceId = 0;
        if (mRank > 0)
        {
            AddNodeIdsToInterface(startPoints, globalNodeId, interfaceId, -1);

            ++interfaceId;
        }

        if (mRank < mNumProcesses - 1)
        {
            AddNodeIdsToInterface(endPoints, globalNodeId, interfaceId, 1);
        }

        return importContainer;
    }
    ///
    std::vector<int> CalculateLagrangeMultiplierIds()
    {

        std::vector<int> lagrangeMultiplierDofIds;

        // add lagrange multipliers from interfaces
        for (const auto& interface : mInterfaces)
            for (const auto& nodeIdPair : interface.mNodeIdsMap)
            {
                /// \todo Think about other DOFs. They must be added too!
                auto dofIds = NodeGetDofIds(nodeIdPair.second, NuTo::Node::eDof::DISPLACEMENTS);

                for (const auto& id : dofIds)
                    lagrangeMultiplierDofIds.push_back(id);
            }

        lagrangeMultiplierDofIds.insert(lagrangeMultiplierDofIds.end(), mBoundaryDofIds.begin(), mBoundaryDofIds.end());

        lagrangeMultiplierDofIds.insert(lagrangeMultiplierDofIds.end(), mPrescribedDisplacementDofIds.begin(),
                                        mPrescribedDisplacementDofIds.end());

        return lagrangeMultiplierDofIds;
    }
    /// \brief Assembles vector for multiplicity scaling
    SparseMatrix MultiplicityScaling()
    {
        if (not(GetDimension() == 2))
            throw Exception(__PRETTY_FUNCTION__, "Multiplicity sclaing only implemented for dimension = 2");

        if (GetDofStatus().GetDofTypes().size() > 1)
            throw Exception(__PRETTY_FUNCTION__, "Multiplicity sclaing not implemented for multiple DOFs");

        // \todo special care needs to be taken for multiple dofs
        const double dim = GetDimension();

        // Vector is initialized with ones because it takes care of all DOFs with Dirichlet BCs
        // Interfaces are treated separately
        Eigen::VectorXd multiplicity = Eigen::VectorXd::Ones(mNumLagrangeMultipliers);

        ReverseMap<int> localNodeIdToGlobalNodeIds;
        for (const auto& interface : mInterfaces)
            localNodeIdToGlobalNodeIds.addMap(interface.mNodeIdsMap);

        for (const auto& pair : localNodeIdToGlobalNodeIds)
        {
            for (const auto& ele : pair.second)
            {
                const int numSubdomainsThatShareThisNode = pair.second.size() + 1;
                multiplicity[dim * ele] = 1. / numSubdomainsThatShareThisNode;
                multiplicity[dim * ele + 1] = 1. / numSubdomainsThatShareThisNode;
            }
        }

        SparseMatrix ScalingMatrix(mNumLagrangeMultipliers, mNumLagrangeMultipliers);
        for (int i = 0; i < mNumLagrangeMultipliers; ++i)
            ScalingMatrix.insert(i, i) = multiplicity[i];

        return ScalingMatrix;
    }
    ///
    const Eigen::VectorXd& GetPrescribedDofVector() const
    {
        return mPrescribedDofVector;
    }
    ///
    const Matrix& GetRigidBodyModes() const
    {
        return mRigidBodyModes;
    }
    ///
    const int GetNumRigidBodyModes() const
    {
        return mNumRigidBodyModes;
    }
    ///
    const int GetNumRigidBodyModesTotal() const
    {
        return mNumRigidBodyModesTotal;
    }
    ///
    const Matrix& GetG() const
    {
        return mG;
    }
    ///
    const SparseMatrix& GetConnectivityMatrix() const
    {
        return mConnectivityMatrix;
    }
    ///
    const Matrix& GetProjectionMatrix() const
    {
        return mProjectionMatrix;
    }

private:
    ///
    /// \param interfaceCoordinates
    /// \param globalNodeId
    /// \param interfaceId
    /// \param connectivityValue
    void AddNodeIdsToInterface(const std::vector<double>& interfaceCoordinates, int& globalNodeId, int& interfaceId,
                               const double connectivityValue)
    {
        const int groupNodesInterface = GroupCreate(eGroupId::Nodes);
        GroupAddNodeCoordinateRange(groupNodesInterface, eDirection::X, interfaceCoordinates[0] - 1.e-6,
                                    interfaceCoordinates[0] + 1.e-6);


        std::vector<int> nodeIdsInterface = GroupGetMemberIds(groupNodesInterface);


        mInterfaces[interfaceId].mValue = connectivityValue;


        for (size_t i = 0; i < nodeIdsInterface.size(); ++i)
        {
            mInterfaces[interfaceId].mNodeIdsMap.emplace(globalNodeId, nodeIdsInterface[i]);
            ++globalNodeId;
        }
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  MEMBER VARIABLES
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

public:
    /// \brief MPI_Comm_rank(...)
    const int mRank;

    /// \brief MPI_Comm_size(...)
    const int mNumProcesses;

    int mNumInterfaceNodesTotal = -1337;
    int mNumTotalPrescribedDisplacementDofIds = 0;
    int mNumLagrangeMultipliers = 0;

protected:
    NodeList mNodes;
    ElementList mElements;
    InterfaceList mInterfaces;
    std::set<int> mSubdomainBoundaryNodeIds;

    /// \brief mG Natural coarse space
    Matrix mG;

    /// \brief mProjectionMatrix P = I - G (G^T G)^{-1} G^T
    Matrix mProjectionMatrix;

    int mNumRigidBodyModes = -1337;
    int mNumRigidBodyModesTotal = -1337;

    /// \brief An array of ids of all the degrees of freedom that are constraint in the subdomain
    std::vector<int> mBoundaryDofIds;

    Eigen::VectorXd mPrescribedDofVector;

    /// \brief mPrescribedDisplacementDofIds
    std::vector<int> mPrescribedDisplacementDofIds;

    ///
    /// \brief mGlobalStartIndexBoundaryDofIds
    ///
    /// This variable is important for the assembly of the connectivity matrix. It is equal to the sum of the boundary
    /// DOFs of all subdomains with lower rank. E.g.
    ///
    /// Subdomain 0:  mBoundaryDofIds.size() = 3: =>  mGlobalStartIndexBoundaryDofIds = 0:
    /// Subdomain 1:  mBoundaryDofIds.size() = 6: =>  mGlobalStartIndexBoundaryDofIds = 3:
    /// Subdomain 2:  mBoundaryDofIds.size() = 2: =>  mGlobalStartIndexBoundaryDofIds = 9:
    ///
    int mGlobalStartIndexBoundaryDofIds = -1337;

    int mGlobalStartIndexPrescribedDisplacementDofIds = -1337;

    /// \brief The total number of degrees of freedom that are constraint in all subdomains
    int mNumTotalBoundaryDofIds = -1337;

    Matrix mRigidBodyModes;

    /// \brief mInterfaceRigidBodyModes B_s * R_s
    Matrix mInterfaceRigidBodyModes;

    /// \brief Sparse matrix that restricts the DOFs of one subdomain to its interface boundary. Contains only (0,1,-1)
    SparseMatrix mConnectivityMatrix;
};
} // namespace NuTo
