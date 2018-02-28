#include <array>
#include <string>


template<int dim>
using zoltanVector_double = Eigen::Matrix<double, dim, 1>;

template<int dim>
using zoltanVector_int = Eigen::Matrix<int, dim, 1>;


template<int valueCount, int dofNumberCount>
class zoltanNode
{
private:
    int mId;
    zoltanVector_double<valueCount> mValues;
    zoltanVector_int<dofNumberCount> mDofNumbers;
    bool mMasterNode;

public:
    zoltanNode() : mId(0), mValues(zoltanVector_double<valueCount>::Constant(valueCount, 0.0)), mMasterNode(false) {}
    zoltanNode(int rId, double rValue) : mId(rId), mValues(zoltanVector_double<valueCount>::Constant(valueCount, rValue)), mMasterNode(false) {}
    zoltanNode(int rId, zoltanVector_double<valueCount> rValues) : mId(rId), mValues(rValues), mMasterNode(false) {}
//    zoltanNode(int rId, Eigen::Vector4d rValues) : mId(rId), mValues(rValues) {}

    zoltanNode& operator=(zoltanNode rNode)
    {
        if (this != &rNode)
        {
            this->mId = rNode.mId;
            this->mValues = rNode.mValues;
            this->mDofNumbers = rNode.mDofNumbers;
            this->mMasterNode = rNode.mMasterNode;
        }
        return *this;
    }

    int getID()
    {
        return mId;
    }

    zoltanVector_double<valueCount> getValues()
    {
        return mValues;
    }

    int getNumValues()
    {
        return valueCount;
    }

    int getNumDofs()
    {
        return dofNumberCount;
    }

    void setDofNumber(int rIndex, int rNumber)
    {
        assert(rIndex < mDofNumbers.rows());
        mDofNumbers[rIndex] = rNumber;
    }

    int getDofNumber(int rIndex)
    {
        assert(rIndex < mDofNumbers.rows());
        return mDofNumbers[rIndex];
    }

    zoltanVector_int<dofNumberCount> getDofNumbers()
    {
        return mDofNumbers;
    }

    void setIsMaster()
    {
        mMasterNode = true;
    }

    bool isMaster()
    {
        return mMasterNode;
    }

    std::string toString()
    {
        std::string output;

        output.append("ID(" + std::to_string(mId) + "), ");
        output.append("dofs(");
        for (int i = 0; i < dofNumberCount; ++i)
        {
            output.append(std::to_string(mDofNumbers(i)));
            if (i < dofNumberCount-1)
                output.append(", ");
        }
        output.append(")");

        return output;
    }

};

class zoltanElement
{
//protected:
public:
    int mId;

    zoltanElement() = default;
    zoltanElement(int rId) : mId(rId) {}

    virtual ~zoltanElement() = default;
};

template<int nodeCount, int nodeValueCount, int nodeDofCount>
class zoltanElementFEM : public zoltanElement
{
private:
//    int mId;
    std::array<zoltanNode<nodeValueCount,nodeDofCount>, nodeCount> mNodes;

public:
    zoltanElementFEM() = default;

////    zoltanElementFEM(std::vector<NuTo::NodeSimple> rNodes, const NuTo::InterpolationQuadLinear& rInterpolation) : mInterpolation(rInterpolation)
//    zoltanElementFEM(std::vector<NuTo::NodeSimple> rNodes, const NuTo::InterpolationSimple& rInterpolation)
//    {
//        mInterpolation = rInterpolation.Clone();
//        for (NuTo::NodeSimple node : rNodes)
//            mNodes.push_back(node);
//    }

////    zoltanElementFEM(int rId, std::vector<NuTo::NodeSimple> rNodes, const NuTo::InterpolationQuadLinear& rInterpolation) : mId(rId),  mInterpolation(rInterpolation)
//    zoltanElementFEM(int rId, std::vector<NuTo::NodeSimple> rNodes, const NuTo::InterpolationSimple& rInterpolation) : mId(rId)
//    {
//        mInterpolation = rInterpolation.Clone();
//        for (NuTo::NodeSimple node : rNodes)
//            mNodes.push_back(node);
//    }

    zoltanElementFEM(std::array<zoltanNode<nodeValueCount, nodeDofCount>, nodeCount> rNodes) : zoltanElement(0)// : mId(0)
    {

        for (int i = 0; i < nodeCount; ++i)
        {
            mNodes[i] = rNodes[i];
        }
    }

    zoltanElementFEM(int rId, std::array<zoltanNode<nodeValueCount, nodeDofCount>, nodeCount> rNodes) : zoltanElement(rId)// : mId(rId)
    {
        for (int i = 0; i < nodeCount; ++i)
        {
            mNodes[i] = rNodes[i];
        }
    }

    zoltanElementFEM(const zoltanElementFEM& rElement) : zoltanElement(rElement.mId)//: mId(rElement.mId)
    {
        for (int i = 0; i < nodeCount; ++i)
        {
            mNodes[i] = rElement.mNodes[i];
        }
    }

    zoltanElementFEM& operator=(zoltanElementFEM rElem)
    {
        if (this != &rElem)
        {
            this->mId = rElem.mId;
            this->mNodes = rElem.mNodes;
        }
        return *this;
    }

//    zoltanElementFEM(const zoltanElementFEM& rElement) : mId(rElement.mId), mInterpolation(rElement.mInterpolation.get()->Clone()), mNodes(rElement.mNodes)
//    {
////        for (NuTo::NodeSimple node : rElement.mNodes)
////        {
////            this->mNodes.push_back(node);
//////            NuTo::NodeSimple newNode(node.GetValues());
//////            for (int i = 0; i < node.mDofNumbers.rows(); ++i)
//////            {
//////                newNode.SetDofNumber(i, node.GetDofNumber(i));
//////            }
//////            this->mNodes.push_back(newNode);
////        }
//    }

    int getID()
    {
        return mId;
    }

    std::array<zoltanNode<nodeValueCount, nodeDofCount>, nodeCount> getNodes()
    {
        return mNodes;
    }

    int getDofDimension()
    {
        return getNode(0).getNumDofs();
    }

    int getNumNodes()
    {
//        return mNodes.size();
        return nodeCount;
    }

//    const std::unique_ptr<NuTo::InterpolationSimple> Interpolation() const
//    {
//        return mInterpolation;
//    }

//    const NuTo::InterpolationSimple& Interpolation()
//    {
//        return *mInterpolation.get();
//    }

    zoltanNode<nodeValueCount, nodeDofCount>& getNode(int i)
    {
        assert(i < static_cast<int>(mNodes.size()));
        return mNodes[i];
    }

    constexpr int getNodeValueCount()
    {
        return nodeValueCount;
    }

    constexpr int getNodeDofCount()
    {
        return nodeDofCount;
    }

};


using zoltanElementFEM_1D = zoltanElementFEM<2, 1, 1>;
using zoltanElementFEM_2D_QUAD = zoltanElementFEM<4, 2, 2>;






