#pragma once

#include <eigen3/Eigen/Dense>


namespace boost
{

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::binary_oarchive& ar,
               const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m, const unsigned int version);
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::binary_iarchive& ar,
               const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m, const unsigned int version);

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::xml_oarchive& ar,
               const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m, const unsigned int version);
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::xml_iarchive& ar,
               const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m, const unsigned int version);

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::text_oarchive& ar,
               const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m, const unsigned int version);
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(boost::archive::text_iarchive& ar,
               const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m, const unsigned int version);

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void save(Archive& ar, const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m,
          const unsigned int version)
{
    int rows = m.rows();
    int cols = m.cols();
    ar& BOOST_SERIALIZATION_NVP(rows);
    ar& BOOST_SERIALIZATION_NVP(cols);
    ar& boost::serialization::make_nvp("EigenMatrix", boost::serialization::make_array(m.data(), rows * cols));
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void load(Archive& ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m,
          const unsigned int version)
{
    int rows, cols;
    ar& BOOST_SERIALIZATION_NVP(rows);
    ar& BOOST_SERIALIZATION_NVP(cols);
    m.resize(rows, cols);
    ar& boost::serialization::make_nvp("EigenMatrix", boost::serialization::make_array(m.data(), rows * cols));
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(Archive& ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m,
               const unsigned int version)
{
    boost::serialization::split_free(ar, m, version);
}


namespace serialization
{

template <class Archive, class TKey, class TUniquePtrClass, class THash>
inline void save(Archive& rArchive, const std::unordered_map<TKey, std::unique_ptr<TUniquePtrClass>, THash>& rMap,
                 const unsigned int /* file_version */)
{

    int size = rMap.size();
    rArchive << boost::serialization::make_nvp("size", size);

    for (const auto& pair : rMap)
    {
        rArchive << boost::serialization::make_nvp("key", pair.first);
        rArchive << boost::serialization::make_nvp("value", pair.second);
    }
}

template <class Archive, class TKey, class TUniquePtrClass, class THash>
inline void load(Archive& rArchive, std::unordered_map<TKey, std::unique_ptr<TUniquePtrClass>, THash>& rMap,
                 const unsigned int /* file_version */)
{

    int size = 0;
    rArchive >> boost::serialization::make_nvp("size", size);

    rMap.clear();
    for (int i = 0; i < size; ++i)
    {
        TKey key;
        std::unique_ptr<TUniquePtrClass> value;

        rArchive >> boost::serialization::make_nvp("key", key);
        rArchive >> boost::serialization::make_nvp("value", value);

        rMap[key] = std::move(value);
    }
}

template <class Archive, class TKey, class TUniquePtrClass, class THash>
inline void serialize(Archive& rArchive, std::unordered_map<TKey, std::unique_ptr<TUniquePtrClass>, THash>& rMap,
                      const unsigned int rFileVersion)
{
    boost::serialization::split_free(rArchive, rMap, rFileVersion);
}

} // namespace serialization
} // namespace boost
