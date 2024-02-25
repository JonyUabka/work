/**
*	@file	ISL_Filter.h
*	@brief	[Header file of Filter Functions], 滤波；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-06-03
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_FILTER_H_
#define PAI_FRAME_ISEISLIB_FILTER_H_

#include "ISL_UserDefine.h"

namespace ISLib {

/**
 * @brief	提取滤波因子(for ISL_bandpassFilter)
 *
 * @param[in]	f1			左下；
 * @param[in]	f2			左上；
 * @param[in]	f3			右下；
 * @param[in]	f4			右上；
 * @param[in]	sim			采样间隔(milli-second)
 * @param[in]	il			子波长度
 * @param[out]	filbuf		输出的滤波因子数组
 * @param[in]	scl			time interval scale value，default = 0.001
 *
 * @return	0为正常
 */
ISL_EXPORT int ISL_filterFactor(float f1,			/* left lower corner frequency*/
							float f2,			/* left upper corner frequency */
							float f3,			/* right lower corner frequency*/
							float f4,			/* right upper corner frequency */
							float sim,			/*time interval (milli-second)  */
							int il,				/* filter wavelet length */
							float *&filbuf,		/* filter waelet buffer */
							float scl = 0.001	/* time interval scale value */
							);


/**
 * @brief	滤波或褶积
 *
 * @param[in]	filter		输入数据的采样点；
 * @param[in]	len			输入的子波数组长度；
 * @param[in]	orig		输入的地震数据数组；
 * @param[in]	num			输入的地震数据数组长度
 * @param[out]	out			输出计算得到结果
 * @param[in]	flag		输入的处理类型，1为滤波，2为褶积
 *
 * @return	0为正常
 */
ISL_EXPORT int ISL_filter(float *filter,	/* filter waelet buffer */
				int len,		/* filter wavelet length */
				float *orig,	/* input data */
				int num,		/* data length */
				float *out,		/* output data , if NULL, then out in orig */
				int flag = 1 	/* process flag : 1 - filter ,2 - convolution*/
				);



/**
 * @brief	带通滤波
 * 			Band-Pass or High-Pass: f1 >= 0.0 , f2 > f1 , f3 > f2 , f4 > f3;
 * 			Low-Pass: f1 >= 0.0 , f2 > f1 , f3 = 0.0, f4 = 0.0;
 *
 * @param[in]	vals		输入数据
 * @param[in]	num			输入数据的长度
 * @param[in]	dtime		time interval (milli-second)
 * @param[in]	f1			左下；
 * @param[in]	f2			左上；
 * @param[in]	f3			右下；
 * @param[in]	f4			右上；
 * @param[len]	len			输出的子波的长度
 * @param[in]	scl			time interval scale value
 *
 * @return	子波数组（长度为 len）
 */
ISL_EXPORT float * ISL_bandpassFilter(float *vals,
							int num,
							float dtime,		/* time interval (milli-second)  */
							float f1,           /* left lower corner frequency*/
							float f2,           /* left upper corner frequency */
							float f3,           /* right lower corner frequency*/
							float f4,           /* right upper corner frequency */
							int len = 128,		/* filter wavelet length */
							float scl = 0.001 	/* time interval scale value */
							);


/**
 * @brief	巴特沃斯滤波器 (Butterworth filter)
 * 			compute number of poles and -3 db frequency for a low-pass or high-pass filter,
 * 			given a frequency response constrained at two frequencies.
 *
 * @param[in]	fpass		frequency in pass band at which amplitude is >= apass
 * @param[in]	apass		amplitude in pass band corresponding to frequency fpass
 * @param[in]	fstop		frequency in stop band at which amplitude is <= astop
 * @param[in]	astop		amplitude in stop band corresponding to frequency fstop
 * @param[out]	npoles		number of poles
 * @param[out]	f3db		frequency at which amplitude is sqrt(0.5) (-3 db)
 *
 * @return	NO
 * 	(1) Nyquist frequency equals 0.5
 *
 *  (2) The following conditions must be true:
 *  	(0.0<fpass && fpass<0.5) &&
 *  	(0.0<fstop && fstop<0.5) &&
 *  	(fpass!=fstop) &&
 *  	(0.0<astop && astop<apass && apass<1.0)
 *
 *  (3) if (fpass<fstop)
 *  		a low-pass filter is assumed
 *  	else
 *  		a high-pass filter is assumed
 */
ISL_EXPORT void ISL_bfdesign (float fpass, float apass, float fstop, float astop, int &npoles, float &f3db);


/**
 * @brief	巴特沃斯-高通
 *
 * @param[in]	npoles		number of poles (and zeros); npoles>=0 is required
 * @param[in]	f3db		3 db frequency; nyquist = 0.5; 0.0<=f3db<=0.5 is required
 * @param[in]	n			length of p and q
 * @param[in]	p[]			array[n] to be filtered
 * @param[out]	q[]			filtered array[n] (may be equivalent to p)
 *
 * @return	NO
 */
ISL_EXPORT void ISL_bfhighpass (int npoles, float f3db, int n, float p[], float q[]);



/**
 * @brief	巴特沃斯-低通
 *
 * @param[in]	npoles		number of poles (and zeros); npoles>=0 is required
 * @param[in]	f3db		3 db frequency; nyquist = 0.5; 0.0<=f3db<=0.5 is required
 * @param[in]	n			length of p and q
 * @param[in]	p[]			array[n] to be filtered
 * @param[out]	q[]			filtered array[n] (may be equivalent to p)
 *
 * @return	NO
 */
ISL_EXPORT void ISL_bflowpass (int npoles, float f3db, int n, float p[], float q[]);


/**
 * @brief	离散随机线性系统的卡尔曼滤波
 *
 * @param[in]	n				动态系统的维数
 * @param[in]	m				观测系统的维数
 * @param[in]	k				观测序列的长度
 * @param[in]	f[n][n]			系统状态转移矩阵
 * @param[in]	q[n][n]			模型噪声Wk的协方差阵
 * @param[in]	r[n][n]			观测噪声Vk的协方差阵
 * @param[in]	h[n][n]			观测矩阵
 * @param[in]	y[n][n]			观测向量序列。其中y(i, j)(i=0,1,...,k=1;j=0,1,...,m=1), 表示第i时刻的观测向量的第j个分量。
 * @param[in]	x[n][n]			其中x(0, j)(j=0,1,...,n=1)存放给定的初值，其余各行返回状态向量估值序列。
 *								X(i, j)(i=0,1,...,k=1;j=0,1,...,n=1)表示第i时刻的观测向量的第j个分量。
 * @param[in]	p[n][n]			存放初值P0。
 * @param[out]	p[n][n]			返回最后时刻的估计误差协方差阵。
 * @param[out]	g[n][n]			返回最后时刻的稳定增益矩阵。
 *
 * @return	函数返回整形标志。若返回标志值为0，则说明求增益矩阵过程中求逆失败；
 *  		若返回标志值不为0，则说明正常返回;
 */
ISL_EXPORT int ISL_lman(int n, int m, int k,
			double *f, double *q, double *r, double *h, double *y,
			double *x, double *p, double *g);



/**
 * @brief	alpha-beta-gamma 滤波
 *
 * @param[in]	n			量测数据的点数
 * @param[in]	x[n]		N个等间隔点上的量测值
 * @param[in]	t			采样周期
 * @param[in]	a			滤波器结构参数-alpha
 * @param[in]	b			滤波器结构参数-beta
 * @param[in]	c			滤波器结构参数-gamma
 * @param[out]	y[n]		返回n个等间隔点上的滤波估值
 *
 * @return	NO
 */
ISL_EXPORT void ISL_kabg(int n, double *x, double t, double a, double b, double c, double *y);




/**
* @brief	计算中值滤波,取一定的时窗对取出的数据进行滤波
*
* @param[in]	inData[num]		输入数据；
* @param[in]	num				输入数据的大小；
* @param[in]	mfl				纵向中值滤波去噪宽度 Medium Filter Long
* @param[out]	midValue		输出计算得到的中值;
*
* @return	no
*/
ISL_EXPORT void ISL_windowMidValue(float *inData, int num, int mfl, float *outData);




/**
* @brief	Dip Filter
*
* @param[in]	nt					每道样点数
* @param[in]	ntr					道数
* @param[in]	inData[ntr][nt]		输入数据
* @param[in]	dt					time sampling interval(s)
* @param[in]	dtr					time sampling interval
* @param[in]	nslopes				number of slopes & amplitudes specified
* @param[in]	slopes[4]			slopes at which amplitudes are specified, Array num = 4
* @param[in]	amps[4]				amplitudes corresponding to slopes, Array num = 4
* @param[in]	bias				slope bias
* @param[out]	ntfft				nt after padding for FFT(>=nt)
* @param[out]	ntrfft				ntr after padding for FFT(>=ntr)
* @param[out]	outData				输出数据(数据大小为ntfft*ntrfft，实际输出大小请用nt*ntr)
*
* @return	no
*/
ISL_EXPORT int ISL_dipFilter(	int nt, int ntr, float ** inData,
								float dt, float dtr,
								int nslopes, float *slopes, float *amps,
								float bias,
								int &ntfft, int &ntrfft, float ** & outData );



/**
* @brief	K Filter
*
* @param[in]	nt					每道样点数
* @param[in]	ntr					道数
* @param[in]	inData[ntr][nt]		输入数据
* @param[in]	dt					time sampling interval(s)
* @param[in]	dtr					time sampling interval
* @param[in]	knum				number of points, amps & defining k filter
* @param[in]	k[4]				wavenumber values defining k filter, Array num = 4
* @param[in]	amps[4]				amplitude values defining k filter, Array num = 4
* @param[in]	nphi				number of phi values for polartorect
* @param[out]	ntfft				nt after padding for FFT(>=nt)
* @param[out]	ntrfft				ntr after padding for FFT(>=ntr)
* @param[out]	outData				输出数据(数据大小为ntfft*ntrfft，实际输出大小请用nt*ntr)
*
* @return	no
*/
ISL_EXPORT int ISL_kFilter(		int nt, int ntr, float ** inData,
								float dt, float dtr,
								int knum, float *k, float *amps,
								int nphi,
								int &ntfft, int &ntrfft, float ** & outData );



/**
* @brief	Wiener predictive error filtering
*
* @param[in]	nt					每道样点数
* @param[in]	dt					time sampling interval(s)
* @param[in]	inData[nt]			输入数据
* @param[out]	outData[nt]			输出数据
*
* @return	no
*/
ISL_EXPORT void ISL_pefFilter(	int nt, float dt, float * inData, float * & outData	);



/**
* @brief	applies a zero-phase, sine-squared tapered filter
*
* @param[in]	nt					每道样点数
* @param[in]	dt					time sampling interval(s)
* @param[in]	inData[nt]			输入数据
* @param[in]	fnum				number of points, amps & f
* @param[in]	f					array of filter frequencies
* @param[in]	amp					array of amplitude values
* @param[out]	outData[nt]			输出数据
*
* @return	no
* " Examples of filters:									",
* " Bandpass:   f=10,20,40,50 | ...							",
* " Bandreject: f=10,20,30,40 amps=1.,0.,0.,1. | ..			",
* " Lowpass:    f=10,20,40,50 amps=1.,1.,0.,0. | ...		",
* " Highpass:   f=10,20,40,50 amps=0.,0.,1.,1. | ...		",
* " Notch:      f=10,12.5,35,50,60 amps=1.,.5,0.,.5,1. |..	",
*/
ISL_EXPORT void ISL_taperFilter(	int nt, float dt, float * inData,
									float fnum, float *f, float *amps,
									float * & outData	);


/**
* @brief	applies a zero-phase, sine-squared tapered filter
*
* @param[in]	nt					每道样点数
* @param[in]	dt					time sampling interval(s)
* @param[in]	inData[nt]			输入数据
* @param[in]	nfilter				number of filters specified
* @param[in]	*tf					times at which filters are centered
* @param[in]	**f					4 corner frequency array for each times
* @param[out]	outData[nt]			输出数据
*
* @return	no
*/
ISL_EXPORT void ISL_tvBandFilter(	int nt, float dt, float * inData,
									int nfilter, float *tf, float **f,
									float * & outData	);


} /* End Of namespace ISLIB */
#endif /* PAI_FRAME_ISEISLIB_FILTER_H_ */
