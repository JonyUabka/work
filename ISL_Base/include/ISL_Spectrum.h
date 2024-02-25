/**
*	@file	ISL_Spectrum.h
*	@brief	[Header file of Spectrum Functions], 频谱分析；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-02-10
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_SPECTRUM_H
#define PAI_FRAME_ISEISLIB_SPECTRUM_H

#include "ISL_FFT.h"
#include "ISL_Interpolation.h"

namespace ISLib {

class ISL_EXPORT Spectrum{

public:
	Spectrum() { ; }
	~Spectrum() { ; }

	/**
	* @brief	计算功率谱  （方法1）
	* @param[in]	data		输入的数组
	* @param[in]	num			输入的数组的个数
	* @param[in]	dtime		输入的地震数据的采样间隔，以微秒为单位
	* @param[in]	flag		输入的标志位，是否进行振幅归一化，1-是，2-否
	*
	* @param[out]	amp			输出的幅值数据
	* @param[out]	freqs		输出的频率数据
	*
	* @return	返回输出数组的个数，个数为零表示计算失败
	*/
	static int ISL_powerSpectrum1(float *data, /* input seismic data    */
							int num, /* number of samples                 */
							float dtime, /* sampling interval in milli-second */
							float *&amps, /* amplitude values                  */
							float *&freqs, /* frequency values                  */
							int flag = 1 /* 是否进行振幅归一化，1-是，2-否    */
							);


	/**
	* @brief	计算功率谱  （方法2）
	* @param[in]	data		输入的数组
	* @param[in]	num			输入的数组的个数
	* @param[in]	dtime		输入的地震数据的采样间隔，以微秒为单位
	* @param[in]	freq_step	频率的步长
	* @param[in]	flag		输入的标志位，是否进行振幅归一化，1-是，2-否
	*
	* @param[out]	amp			输出的幅值数据
	* @param[out]	freqs		输出的频率数据
	*
	* @return	返回输出数组的个数，个数为零表示计算失败
	*/
	static int ISL_powerSpectrum2(float *data, /* input seismic data    */
							int num, /* number of samples                 */
							float dtime, /* sampling interval in milli-second */
							float *&amps, /* amplitude values                  */
							float *&freqs, /* frequency values                  */
							float freq_step, /* frequency step                    */
							int flag = 1 /* 是否进行振幅归一化，1-是，2-否    */
							);

};




}/* End Of namespace ISLIB */
#endif /* ISL_SPECTRUM_H_ */
