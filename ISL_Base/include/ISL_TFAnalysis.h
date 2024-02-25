/**
*	@file	ISL_TFAnalysis.h
*	@brief	[Header file of time-frequency analysis Functions], 时频分析；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang, Chen Ke], 刘百红、杨强、宋志翔、陈科；
*	@date	2014-06-09
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_TFANALYSIS_H_
#define PAI_FRAME_ISEISLIB_TFANALYSIS_H_

#include "ISL_UserDefine.h"

namespace ISLib {

/**
 * @brief	广义S变换
 *
 * @param[in]	inSignal[nPoints]	输入的地震信号；
 * @param[in]	nPoints				输入地震信号的长度；
 * @param[in]	fDT					输入地震信号的采样间隔；(秒)
 * @param[in]	fAlpha				广义S变换的参数(0.3f)；
 * @param[in]	nDelayTime			子波的时间延迟
 * @param[in]	nFreqStrat			起始频率
 * @param[in]	nFreqEnd			终止频率
 * @param[in]	nFreqInterval		频率间隔
 * @param[out]	gstAmpSpec[nPoints][nF]		输出广义S变换求得的时频谱, nF = （起始频率-终止频率）/频率间隔+1
 *
 * @return	no
 */
ISL_EXPORT void ISL_gsTransform(float * inSignal,
						int nPoints,
						float fDT,
						float fAlpha,
						int nDelayTime,
						int nFreqStrat,
						int nFreqEnd,
						int nFreqInterval,
						float ** outAmpSpec);



/**
 * @brief	三参数小波变换
 *
 * @param[in]	inSignal[nPoints]	输入的地震信号；
 * @param[in]	nPoints				输入地震信号的长度；
 * @param[in]	fDT					输入地震信号的采样间隔；(秒)
 * @param[in]	nScales				最大尺度
 * @param[in]	nFreqStrat			起始频率
 * @param[in]	nFreqEnd			终止频率
 * @param[in]	nFreqInterval		频率间隔
 * @param[in]	fWavRe				母小波的实部
 * @param[in]	fWavIm              母小波的虚部
 * @param[in]	nWav                母小波的长度
 * @param[in]	fFre				中心频率
 * @param[out]	outAmpSpec[nPoints][nF]		输出三参数小波变换求得的时频谱,nF = （起始频率-终止频率）/频率间隔+1
 *
 * @return	no
 */
ISL_EXPORT void ISL_triParWTransform(float *inSignal,
							int nPoints,
							float fDT,
							int nFreqStrat,
							int nFreqEnd,
							int nFreqInterval,
							float *fWavRe,
							float *fWavIm,
							int nWav,
							float fFre,
							float **outAmpSpec);




/**
 * @brief	短时傅里叶变换, 输入一个信号，输出信号的短时傅里叶变换的单频谱
 *
 * @param[in]	inSignal[nPoints]		输入的地震信号；
 * @param[in]	nPoints					输入地震信号的长度；
 * @param[in]	fDT						输入地震信号的采样间隔；(秒)
 * @param[in]	inWinData[nWinLength]	输入的STFT的时窗的数据
 * @param[in]	nWinLength				输入的STFT的时窗长度
 * @param[in]	nFreqStrat				起始频率
 * @param[in]	nFreqEnd				终止频率
 * @param[in]	nFreqInterval			频率间隔
 * @param[out]	outAmpSpec[nPoints][nF]		输出求得的单频振幅谱,nF = （起始频率-终止频率）/频率间隔+1
 *
 * @return	no
 */
ISL_EXPORT void ISL_uniFreSTFTransform(float *inSignal,
							int nPoints,
							float fDT,
							float *inWinData,
							int nWinLength,
							int nFreqStrat,
							int nFreqEnd,
							int nFreqInterval,
							float ** outAmpSpec);




/**
 * @brief	小波变换, 输入一个信号，输出信号的小波变换的单频谱
 *
 * @param[in]	inSignal[nPoints]			输入的地震信号；
 * @param[in]	nPoints						输入地震信号的长度
 * @param[in]	fDT							输入地震信号的采样间隔；(秒)
 * @param[in]	freq_start					起始频率；
 * @param[in]	freq_end					终了频率；
 * @param[in]	freq_interval				频率间隔；
 * @param[in]	nWaveletLen					子波的长度
 * @param[in]	nWaveletDT					子波的采样率(秒)
 * @param[in]	nGrade						morlet子波的级别, 一般填6
 * @param[in]	nWFlag						子波类型，1=Rick子波，2=零相位宽带雷克子波，3=Morlet子波
 * @param[out]	outAmpSpec[nPoints][nDot]	输出求得的单频振幅谱
 *
 * @return	no
 */
ISL_EXPORT void ISL_uniFreWTransform(float *inSignal,
							int nPoints,
							float fDT,
							int nFreqStrat,
							int nFreqEnd,
							int nFreqInterval,
							int nWaveletLen,
							float nWaveletDT,
							float ** outAmpSpec,
							int nGrade = 6,
							int nWFlag = 1);



} /* End Of namespace ISLIB */



#endif /* PAI_FRAME_ISEISLIB_TFANALYSIS_H_ */
