/**
*	@file	ISL_NumericalIntegration.h
*	@brief	[Header file of Numerical Integration Functions], 数值积分；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-06-24
*	@refer
*/


#ifndef PAI_FRAME_ISEISLIB__NUMERICALINTEGRATION_H_
#define PAI_FRAME_ISEISLIB__NUMERICALINTEGRATION_H_

#include "ISL_UserDefine.h"

namespace ISLib {

/**
* @brief	变步长梯形求积法
*
* @param[in]	a			积分下限；
* @param[in]	b			积分上限；
* @param[in]	eps			积分精度(0.000001)；
*
* @return	返回积分的值（double）
*/
ISL_EXPORT double ISL_fffts(double a, double b, double eps);



/**
* @brief	变步长辛卜森求积法
*
* @param[in]	a			积分下限；
* @param[in]	b			积分上限；
* @param[in]	eps			积分精度(0.000001)；
*
* @return	返回积分的值（double）
*/
ISL_EXPORT double ISL_fsimp(double a, double b, double eps);


} /* End Of namespace ISLIB */
#endif /* PAI_FRAME_ISEISLIB__NUMERICALINTEGRATION_H_ */
