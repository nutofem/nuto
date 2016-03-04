#ifndef EIGEN_BOOST_SERIALIZATION
#define EIGEN_BOOST_SERIALIZATION

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_free.hpp>
#endif  // ENABLE_SERIALIZATION

#include <eigen3/Eigen/Dense>


namespace boost{

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::binary_oarchive & ar, const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version);
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::binary_iarchive & ar, const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version);

template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::xml_oarchive & ar, const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version);
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::xml_iarchive & ar, const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version);

template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::text_oarchive & ar, const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version);
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::text_iarchive & ar, const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version);

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void save(Archive & ar, const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version)
{
    int rows = m.rows();
    int cols = m.cols();
    ar & BOOST_SERIALIZATION_NVP(rows);
    ar & BOOST_SERIALIZATION_NVP(cols);
    ar &  boost::serialization::make_nvp("EigenMatrix", boost::serialization::make_array(m.data(), rows*cols));
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void load(Archive & ar, Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version)
{
    int rows,cols;
    ar & BOOST_SERIALIZATION_NVP(rows);
    ar & BOOST_SERIALIZATION_NVP(cols);
    m.resize(rows,cols);
    ar &  boost::serialization::make_nvp("EigenMatrix", boost::serialization::make_array(m.data(), rows*cols));
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(Archive & ar, Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version)
{
    boost::serialization::split_free(ar,m,version);
}

//template<class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
//inline void serialize(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & m, const unsigned int file_version)
//{
//    ar & boost::serialization::make_array(m.data(), m.size());
//}

}

#endif
