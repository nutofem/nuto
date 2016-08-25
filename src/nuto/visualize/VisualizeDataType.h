// $Id$

#pragma once
#include <string>

namespace NuTo
{
enum class eVisualizeDataType
{
    SCALAR,      //!< scalar data
    VECTOR,      //!< vector data (vector length 3)
    TENSOR,      //!< tensor data (3x3 tensor)
    FIELD        //!< field data (variable length)
};

//! @brief ... class describing the data type stored at the cells or the points
//! @author Stefan Eckardt, ISM
//! @date November 2009
class VisualizeDataType
{
public:

    // constructor
    VisualizeDataType(const std::string& rIdent, eVisualizeDataType rDataType);

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

    inline eVisualizeDataType GetDataType() const
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
    eVisualizeDataType mDataType;

    //! @brief ... number of data
    unsigned int mNumData;
};

}

