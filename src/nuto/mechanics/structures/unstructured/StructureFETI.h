
#pragma once

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/nodes/NodeEnum.h"


#include <eigen3/Eigen/Sparse>
#include <cstring>
#include <fstream>
#include <set>


namespace NuTo
{


//! @author Philip Huschke
//! @date September 2016
//! @brief Class for FETI methods
class StructureFETI: public Structure
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

    struct Boundary
    {
        std::map<int, int> mNodeIdsMap;
        std::vector<int> mValues;
    };

    struct Interface
    {
        std::map<int, int> mNodeIdsMap;
        std::vector<int> mNodeIds;
        int mValue;
    };

    using NodeList          = std::vector<Node>;
    using ElementList       = std::vector<Element>;
    using BoundaryList      = std::vector<Boundary>;
    using InterfaceList     = std::vector<Interface>;
    using Matrix            = Eigen::MatrixXd;
    using SparseMatrix      = Eigen::SparseMatrix<double>;
    using SparseMatrixMap   = std::map<NuTo::Node::eDof, SparseMatrix>;

    NodeList        ReadNodeData        (std::ifstream& file);
    ElementList     ReadElementData     (std::ifstream& file);
    BoundaryList    ReadBoundaryData    (std::ifstream& file);
    InterfaceList   ReadInterfaceData   (std::ifstream& file);



public:
    //! @brief Constructor
    //! @param rDimension   Structural dimension (1,2 or 3)
    //! @param rMeshFile    mesh file
    StructureFETI(int rDimension);

    const int mRank;
    const int mNumProcesses;

    ///
    /// \brief FindKeywordInFile
    /// \param file
    /// \param keyword
    ///
    void FindKeywordInFile(std::ifstream &file, std::string keyword);


    Matrix&                      GetRigidBodyModes()           {return mRigidBodyModes;}
    Matrix&                      GetInterfaceRigidBodyModes()  {return mInterfaceRigidBodyModes;}
    SparseMatrix&                GetConnectivityMatrix()       {return mConnectivityMatrix;}

//    Eigen::SparseMatrix<double> &AssembleStiffnessMatrix();
//    void AssembleRigidBodyModes();

    ///
    /// \brief AssembleConnectivityMatrix
    ///
    void AssembleConnectivityMatrix();

    ///
    /// \brief AssembleBoundaryDofIds
    ///
    void AssembleBoundaryDofIds();

    Matrix                      mRigidBodyModes;
    Matrix                      mInterfaceRigidBodyModes;
    SparseMatrix                mConnectivityMatrix;

    ///
    /// \brief mBoundaryRowDofIds
    ///
    Eigen::DiagonalMatrix<double, Eigen::Dynamic>                mBoundaryDofIds;




    Matrix                     mProjectionMatrix;
    NodeList                   mNodes;
    ElementList                mElements;
    BoundaryList               mBoundaries;
    InterfaceList              mInterfaces;
    int                        mNumRigidBodyModes;
    int                        mNumInterfaceNodesTotal = 42;
    std::set<int>              mSubdomainBoundaryNodeIds;


    void ImportMeshJson(std::string rFileName, const int interpolationTypeId);
    void ImportMesh(std::string rFileName, const int interpolationTypeId);

protected:





};
} //namespace NuTo

