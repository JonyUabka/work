/**
*	@file	ISL_Normalize.h
*	@brief	[Header file of Normalize], 归一化；
*	@see	ISeisLib Manual
*	@author [Song ZhiXiang], 宋志翔；
*	@date	2014-06-11
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_NORMALIZE_H_
#define PAI_FRAME_ISEISLIB_NORMALIZE_H_

#include "ISL_FFT.h"

namespace ISLib {

/**
* @brief	对FFT的计算结果进行归一化
*
* @param[in]	cpfft[nw]	输入的复数数组（个数为nw）
* @param[out]	cpfft[nw]	输出排序过的复数数组（个数为nw）
* @param[in]	nw			输入的复数数组的长度
* @param[in]	type		计算类型
*
* @return		no
*/
ISL_EXPORT void ISL_normalizeFFT(complex *cpfft, int nw, int type);


/**
* @brief	重排子波
*
* @param[in]	nw			输入的复数数组的长度
* @param[in]	*cpfft		输入的复数数组（个数为nw）
* @param[in]	ntfft
* @param[in]	*pfft		输入子波的数组（个数为ntfft）
*
* @param[out]	*pfft		输出子波的数组,输入时请分配好内存空间
* @param[out]	len			输出子波数组的长度
*
* @return		no
*/
ISL_EXPORT void ISL_arrangeWavelet(int nw, complex *cpfft, int ntfft, float *pfft, int &len);



/**
* @brief	将地震子波做频率域归一化处理
*
* @param[in]	inData[waveLen]		输入的数组
* @param[out]	inData[waveLen]		输出排序后的数据
* @param[in]	waveLen				输入数据的采样点数
*
* @return		no
*/
ISL_EXPORT void ISL_normalizeWaveletFre(float *inData, int waveLen);


} /*End of ISLib*/


#endif /* PAI_FRAME_ISEISLIB_NORMALIZE_H_ */
