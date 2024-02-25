/**
*	@file	ISL_Attribute.h
*	@brief	[Header file of Attribute Functions], 地震属性；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-06-26
*	@refer
*/

#ifndef PAI_FRAME_ISEISLIB_ATTRIBUTE_H_
#define PAI_FRAME_ISEISLIB_ATTRIBUTE_H_

#include "ISL_Export.h"

namespace ISLib {

/**
* @brief	计算三瞬属性, 三瞬属性计算，包括瞬时振幅、瞬时相位和瞬时频率
*
* @param[in]	mode		输入的计算类型，0表示计算瞬时振幅，1表示计算瞬时相位，2表示计算瞬时频率
* @param[in]	in			输入的数组
* @param[in]	num			输入的数组的长度
* @param[in]	delta		输入的地震数据的采样间隔，以秒为单位
* @param[out]	out			输出经希尔伯特变换得到的复信号的虚部(外部需分配好内存空间)
*
* @return	无
*/
ISL_EXPORT void ISL_insAttributes(int mode, float *in, int num, float delta, float *out);



/**
* @brief	相位旋转
*
* @param[in]	in			输入的每道的数据；
* @param[in]	num			输入数据的点数；
* @param[in]	angle		旋转相位的度数；
* @param[out]	out			旋转完相位之后输出的数据;
*
* @return	无
*/
ISL_EXPORT void ISL_alterPh(float *in, int num, int angle, float *out);




} /* End Of namespace ISLIB */
#endif /* ISL_ATTRIBUTE_H_ */
