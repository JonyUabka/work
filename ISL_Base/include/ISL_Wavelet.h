/**
*	@file	ISL_Wavelet.h
*	@brief	[Header file of Wavelet], 子波；
*	@see	ISeisLib Manual
*	@author [Zhang Yang], 张杨；
*	@date	2014-05-26
*	@refer	SU CWP
*/	

#ifndef PAI_FRAME_ISEISLIB_WAVELET_H
#define PAI_FRAME_ISEISLIB_WAVELET_H

#include "ISL_SU_Matrix.h"
#include "ISL_Convolution.h"
#include "ISL_MathBase.h"
#include "ISL_Interpolation.h"
#include "ISL_Smooth.h"

namespace ISLib {

	/**
	* @brief	标准雷克子波
	*
	* @param[in]	nt			雷克子波的采样点数
	* @param[in]	dt			雷克子波的采样间隔(秒）
	* @param[in]	fpeak		雷克子波的主频
	* @param[in]	mode		1=小波变换专用雷克子波，2=经典雷克子波，3=NEWS雷克子波
	*
	* @param[out]	*wavelet	输出的生成的雷克子波
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_classicRicker (int nt, float dt, float fpeak, float *wavelet, int mode = 1);



	/**
	* @brief	生成宽带雷克子波，零相位的
	*
	* @param[in]	nt			雷克子波的采样点数
	* @param[in]	dt			雷克子波的采样间隔(秒）
	* @param[in]	lowCut		低截频率
	* @param[in]	highCut		高截频率
	*
	* @param[out]	*wavelet	输出的生成的雷克子波
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_broadBandRicker(int nt, float dt, float lowCut, float highCut, float *wavelet);


	/**
	* @brief	SU 雷克子波
	*
	* @param[in]	hlw			half-length of the wavelet including center (samples)
	* @param[in]	dt			sampling interval (s)
	* @param[in]	period		wavelet period (s)
	* @param[in]	ampl		wavelet amplitude
	* @param[in]	distort		wavelet distortion factor
	*
	* @param[out]	*wavelet	Ricker wavelet
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_ricker3 (int hlw, float dt, float period, float ampl, float distort, float *wavelet);



	/**
	* @brief	Morlet小波,本函数用于求取Morlet小波,pis2为Pi的根，pis4为pis的根
	*
	* @param[in]	dt			采样率(秒)
	* @param[in]	mainfreq	输入的Morlet小波的主频
	* @param[in]	nt			输入的Morlet小波的点数
	* @param[in]	level		输入的morlet分量的级数
	*
	* @param[out]	*wavelet	Ricker wavelet
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_morlet(float dt, float mainfreq, int nt, int level, float *wavelet);



	/**
	* @brief	Morlet小波,利用时移、相位、尺度来计算Morlet小波
	*
	* @param[in]	dt			采样率(秒)
	* @param[in]	mainfreq	输入的Morlet小波的主频
	* @param[in]	nt			输入的Morlet小波的点数
	* @param[in]	timeDelay	输入的Morlet小波的时移
	* @param[in]	phase		输入的Morlet小波的相位
	* @param[in]	scale		输入的Morlet小波的尺度
	*
	* @param[out]	*wavelet	Ricker wavelet
	*
	* @return		no
	*/
	ISL_EXPORT bool ISL_morlet2(float dt, float mainfreq, int nt, float timeDelay, float phase, float scale, float *wavelet);


	/**
	* @brief	计算地震子波的幅度
	*
	* @param[in]	pdata		地震信号
	* @param[in]	num			地震信号的长度
	* @param[in]	dt			地震信号采样率(秒)
	* @param[in]	mainfreq	输入的Morlet小波的主频
	* @param[in]	timeDelay	输入的Morlet小波的时移
	* @param[in]	phase		输入的Morlet小波的相位
	* @param[in]	scale		输入的Morlet小波的尺度
	*
	* @return		子波的幅度
	*/
	ISL_EXPORT float ISL_morletAmp(float *pdata, int num, float dt, float meanF, float timeDelay, float phase, float scale);


	/**
	* @brief	修改Morlet子波相关参数,通过匹配原则，不断的搜索子波
	*
	* @param[in]	pdata		地震信号
	* @param[in]	num			地震信号的长度
	* @param[in]	dt			地震信号采样率(秒)
	* @param[in/out]	mainfreq	输入的Morlet小波的主频
	* @param[in/out]	timeDelay	输入的Morlet小波的时移
	* @param[in/out]	phase		输入的Morlet小波的相位
	* @param[in/out]	scale		输入的Morlet小波的尺度
	* @param[in]	offsetT		时移偏移量
	* @param[in]	offsetF		频率偏移量
	* @param[in]	offsetP		相位偏移量
	* @param[in]	offsetS		尺度偏移量
	*
	* @return		无
	*/
	ISL_EXPORT void ISL_refineMorlet(float *pData, float dt, int num,
						float &mainfreq, float &timeDelay, float &phase, float &scale,
						float offsetT, float offsetF, float offsetP, float offsetS);



	/**
	* @brief	计算三参数小波
	*
	* @param[in]	fTal		能量衰减因子大小, default = 5.0
	* @param[in]	fSigma		调谐频率大小, default = 5.0
	* @param[in]	nWav		母小波长度,, default = 64
	* @param[out]	fWavRe		母小波的实部
	* @param[out]	fWavIm		母小波的虚部
	*
	* @return	中心频率(float)
	*/
	ISL_EXPORT float ISL_triParWave(float fTal, float fSigma, int nWav, float *fWavRe, float *fWavIm);



	/**
	* @brief	akb子波
	*
	* @param[in]	nt			雷克子波的采样点数
	* @param[in]	dt			雷克子波的采样间隔(秒）
	* @param[in]	fpeak		雷克子波的主频
	*
	* @param[out]	*wavelet	输出的生成的雷克子波
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_akb_wavelet (int nt, float dt, float fpeak, float *wavelet);



	/**
	* @brief	尖脉冲子波
	*
	* @param[in]	nt			number of time step
	* @param[in]	tindex		time index to locate the spike
	*
	* @param[out]	*wavelet	array[nt] of computed wavelet
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_spike_wavelet (int nt, int tindex, float *wavelet);


	/**
	* @brief	unit_wavelet
	*
	* @param[in]	nt			number of samples in output wavelet
	*
	* @param[out]	*wavelet	array[nt] of computed wavelet
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_unit_wavelet (int nt, float *wavelet);


	/**
	* @brief	零相位子波
	*
	* @param[in]	nt			number of samples in output wavelet
	*
	* @param[out]	*wavelet	array[nt] of computed wavelet
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_zero_wavelet (int nt, float *wavelet);



	/**
	* @brief	贝尔拉格子波
	*
	* @param[in]	nt			number of samples in output wavelet
	* @param[in]	dt			time step
	* @param[in]	fpeak		peak frequency of the Berlage wavelet
	* @param[in]	ampl		wavelet amplitude
	* @param[in]	tn			non-negative time exponent (typically an integer number)
	* @param[in]	decay		non-negative exponential decay factor
	* @param[in]	ipa			initial phase angle in radians
	*
	* @param[out]	*wavelet	array[nt] of computed wavelet,Berlage wavelet
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_berlage_wavelet (int nt, float dt, float fpeak, float ampl, float tn,
							float decay, float ipa, float *wavelet);



	/**
	* @brief	高斯子波
	*
	* @param[in]	nt			number of samples in output wavelet
	* @param[in]	dt			time step
	* @param[in]	fpeak		peak frequency of the Gaussian wavelet
	*
	* @param[out]	*wavelet	array[nt] of computed wavelet, Gaussian wavelet
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_gaussian_wavelet (int nt, float dt, float fpeak, float *wavelet);



	/**
	* @brief	高斯一阶导数子波
	*
	* @param[in]	nt			number of samples in output wavelet
	* @param[in]	dt			time step
	* @param[in]	fpeak		peak frequency of the Gaussian wavelet
	*
	* @param[out]	*wavelet	array[nt] of computed Gaussian wavelet, first derivative
	*
	* @return		no
	*/



	ISL_EXPORT void ISL_gaussderiv_wavelet (int nt, float dt, float fpeak, float *wavelet);

	/**
	* @brief	高斯一阶导数子波
	*
	* @param[in]	dt			sampling interval
	* @param[in]	nt			length of waveform in samples
	* @param[in]	t0			time shift for (pseudo-) causality
	* @param[in]	fpeak		maximum frequency
	* @param[in]	n			order of derivative
	* @param[in]	sign		multiplier for polarity of waveform
	* @param[in]	verbose		flag for diagnostic messages
	*
	* @param[out]	*w			array of size nt containing the waveform
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_deriv_n_gauss(double dt, int nt, double t0, float fpeak, int n, double *w,
						int sign, int verbose);


	/**
	* @brief	Bath 震源爆炸子波（单点）
	*
	* @param[in]	A0
	* @param[in]	r			爆破半径  m
	* @param[in]	Vp			纵波速度   m/s
	* @param[in]	z			传播距离  m
	* @param[in]	t			传播时间  s

	*
	* @return		float  		单个杨点样点值
	*/
	ISL_EXPORT float ISL_bathWavelet_SinglePoint (float A0, float r, float Vp, float z, float t);


	/**
	* @brief	Bath 震源爆炸子波（数组）
	*
	* @param[in]	A0
	* @param[in]	r				爆破半径  m
	* @param[in]	Vp				纵波速度   m/s
	* @param[in]	z				传播距离  m
	* @param[in]	nin				输入数组个数
	* @param[in]	tin[nin]		输入的传播时间（ s）数组，个数为nin
	*
	* @param[out]	tout[nin]		输出的计算结果数组，个数为nin

	*
	* @return		no
	*/
	ISL_EXPORT void ISL_bathWavelet (float A0, float r, float Vp, float z, int nin, float tin[], float tout[]);



	/**
	* @brief	地震子波提取 	根据输入的地震信号提取地震子波
	*
	* @param[in]	vvs[num]		输入的地震数据，数组形式
	* @param[in]	num				输入的地震数据长度，不能小于182
	*
	* @param[out]	len				返回提取的地震子波的长度
	*
	* @return		float * 		返回提取的地震子波，得到后请自己负责释放内存空间,数组的大小为n
	*/
	ISL_EXPORT float *ISL_extractSeisWavelet(float *vvs, int num, int &len);



	/**
	* @brief	地震子波提取 	根据输入的地震信号和反射系数提取地震子波
	*
	* @param[in]	vvs[num]		输入的地震数据，数组形式
	* @param[in]	refs[num]		输入的反射系数序列，数组形式
	* @param[in]	num				输入的地震数据和反射系数的长度，不能小于130
	*
	* @param[out]	len				返回提取的地震子波的长度
	*
	* @return		float * 		返回提取的地震子波，得到后请自己负责释放内存空间,数组的大小为n
	*/
	ISL_EXPORT float *ISL_extractSeisWavelet(float *vvs, float *refs, int num, int &len);


	/**
	* @brief	合成记录
	*
	* @param[in]	wavelet_type		=1 for a spike
	*									=2 for Tong Fei's Ricker wavelet
	*									=3 for Larner's Ricker wavelet
	*									=4 for Akima wavelet
	* @param[in]	nx					number of horizontal samples (traces)
	* @param[in]	nt					number of vertical samples (samples/trace)
	* @param[in]	dt					time sampling interval
	* @param[in]	fpeak				frequency peak for Ricker or Akima wavelets
	* @param[in]	**wfield			array[nx][nt] of reflectivities
	*
	* @param[out]	**wfield			array[nx][nt] of seismic traces
	*
	* @return		no
	*/
	ISL_EXPORT void ISL_convolve_wavelet (int wavelet_type, int nx, int nt, float dt, float fpeak,	float **wfield);




	/**
	* @brief	寻找子波的主频和相位
	*
	* @param[in]	wdata[wnum]			子波数据
	* @param[in]	wnum				子波长度
	* @param[in]	dt					时间采样间隔（秒）
	*
	* @param[out]	freq				得到的主频(HZ)
	* @param[out]	ph					得到的相位
	*
	* @return		int  				返回：0 ，表示一切正常；-1，表示程序异常
	*/
	ISL_EXPORT int ISL_waveletFreqAndPh(float * wdata, int wnum, float dt, float &freq, float &ph);


	/**
	* @brief	高斯窗 Gauss Window
	*
	* @param[out]	out[n]			返回的子波数据
	* @param[in]	n				子波长度
	*
	* @return		无
	*/
	ISL_EXPORT void ISL_gaussWindow(float * out, int n);



	/**
	* @brief	汉宁窗  Hanning Window
	*
	* @param[out]	out[n]			返回的子波数据
	* @param[in]	n				子波长度
	*
	* @return		无
	*/
	ISL_EXPORT void ISL_hanningWindow(float * out, int n);



	/**
	* @brief	海明窗  Hamming Window
	*
	* @param[out]	out[n]			返回的子波数据
	* @param[in]	n				子波长度
	*
	* @return		无
	*/
	ISL_EXPORT void ISL_hammingWindow(float * out, int n);



	/**
	* @brief	Blackman Window
	*
	* @param[out]	out[n]			返回的子波数据
	* @param[in]	n				子波长度
	*
	* @return		无
	*/
	ISL_EXPORT void ISL_blackmanWindow(float * out, int n);


} /*End of ISLib*/

#endif
