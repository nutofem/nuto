/*
 * MatrixBaseAddons.h
 *
 *  Created on: Apr 11, 2013
 *      Author: junger
 */

#ifndef MATRIXBASEADDONS_H_
#define MATRIXBASEADDONS_H_

inline Scalar at(uint i, uint j) const { return this->operator()(i,j); }
//! @brief elementwise absolute value of the matrix

template<typename OtherDerived>
inline const Eigen::MatrixBase<OtherDerived> abs() const
{
	return this->array().abs();
}


#endif /* MATRIXBASEADDONS_H_ */
