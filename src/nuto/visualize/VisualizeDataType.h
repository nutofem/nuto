// $Id$

#ifndef VISUALIZEDATATYPE_H_
#define VISUALIZEDATATYPE_H_
#include <string>

namespace NuTo
{

//! @brief ... class describing the data type stored at the cells or the points
//! @author Stefan Eckardt, ISM
//! @date November 2009
class VisualizeDataType
{
public:
    enum eDataType
    {
        SCALAR,      //!< scalar data
        VECTOR,      //!< vector data (vector length 3)
        TENSOR,      //!< tensor data (3x3 tensor)
        FIELD        //!< field data (variable length)
    };

    // constructor
    VisualizeDataType(const std::string& rIdent, eDataType rDataType);

    inline std::string GetIdent() const
    {
        return this->mIdent;
    }

    inline bool IsIdent(const std::string& rIdent) const
    {
    	if (rIdent == this->mIdent)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    inline eDataType GetDataType() const
    {
        return this->mDataType;
    }

    inline unsigned int GetNumData() const
    {
        return this->mNumData;
    }

    void SetNumData(unsigned int& rNumData);

private:
    //! @brief ... identifier
    std::string mIdent;

    //! @brief ... type of data
    eDataType mDataType;

    //! @brief ... number of data
    unsigned int mNumData;
};

}

#endif // VISUALIZEDATATYPE_H_
