/**
*	@file	ISL_SeismicAcquistition.h
*	@brief	[Header file of Seismic Acquistition Functions], 地震采集相关算法；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-07-24
*	@refer	Baidu
*/

#ifndef PAI_FRAME_ISEISLIB__SEISMICACQUISTITION_H_
#define PAI_FRAME_ISEISLIB__SEISMICACQUISTITION_H_


#include "ISL_UserDefine.h"

namespace ISLib {


/**
* @brief	迭加信号响应公式 (ghost reflection stack signal)
*
* @param[in]	r			虚反射接口的反射系数；
* @param[in]	f			激发波的频率；
* @param[in]	ang			激发波传播的方向（与垂直向下的夹角）；
* @param[in]	v			介质层速度；
* @param[in]	d			激发井深；
*
* @return	迭加信号响应（double）
*/
template<typename T>
ISL_EXPORT double ISL_grss(T r, T f, T ang, T v, T d)
{
	double t1 = 1 + r*r + 2*r*cos(4*PI*f*d*cos(ang)/v);
	double b1 = 2*(1 + fabs(r));

	return sqrt(t1 / b1);
}






} /* End Of namespace ISLIB */
#endif /* PAI_FRAME_ISEISLIB__SEISMICACQUISTITION_H_ */
