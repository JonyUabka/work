/**
*	@file	ISL_Transform.h
*	@brief	[Header file of Transform Functions], 数值变换；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang, Chen Ke], 刘百红、杨强、宋志翔、陈科；
*	@date	2014-06-03
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_TRANSFORM_H_
#define PAI_FRAME_ISEISLIB_TRANSFORM_H_

#include "ISL_FFT.h"
#include "ISL_Convolution.h"

namespace ISLib {

#define LHHALF 30	/* half-length of Hilbert transform filter*/
#define LH 2*LHHALF+1	/* filter length must be odd */

enum hilbert_type{
	NEWS_hlb = 0,
	SU_hlb
};

/**
 * @brief	希尔伯特变换
 *
 * @param[in]	num			输入数据的采样点数；
 * @param[in]	in			输入的数组；
 * @param[out]	out			输出经希尔伯特变换得到的复信号的虚部
 * @param[in]	type		计算方法， NEWS_hlb = NEWS方法， SU_hlb = SU的方法(自动归一化)
 *
 * @return	no
 */
ISL_EXPORT void ISL_hilbert(int num, float *in, float * out, hilbert_type type = NEWS_hlb);


/**
 * @brief	希尔伯特变换(返回实部与虚部)
 *
 * @param[in]	num			输入数据的采样点数；
 * @param[in]	in			输入的数组；
 * @param[out]	outRe		输出经希尔伯特变换得到的复信号的实部
 * @param[out]	outIm		输出经希尔伯特变换得到的复信号的虚部
 *
 * @return	no
 */
ISL_EXPORT void ISL_hilbert(int num, float * in, float * outRe, float * outIm);


/**
* @brief	沃什（Walsh）变换
*
* @param[in]	p			存放长度为 n=2^k 的给定输入序列；
* @param[in]	n			输入序列的长度；
* @param[in]	k			满足  n=2^k；
* @param[out]	x			返回输入序列, Pi(i=0, 1, ..., n=1)的沃什变换序列;
*
* @return	no
*/
ISL_EXPORT void ISL_kfwt(double *p, int n, int k, double *x);


/**
* @brief	单频小波变换，输入一个信号，输出信号的小波变换的单频谱
*
* @param[in]	in			输入的地震信号；
* @param[in]	num			输入地震信号的长度；
* @param[in]	wavelet		输入的子波信号；
* @param[in]	nw			输入的子波长度；
* @param[out]	ampSpec		输出求得的单频振幅谱；
* @param[out]	phSpec		输出求得的单频相位谱；
*
* @return	无
*/
ISL_EXPORT void ISL_sfCWT(float *in, int num, float *wavelet, int nw, float *ampSpec, float *phSpec);



/**
* @brief	Z变换, 在Z平面单位圆上计算有限长序列x(n)的Z变换的采样值
*
* @param[in]	xr		双精度实型数组，长度大于或等于（n+m-1），并且是2的整数次幂。输入数据的实部；
* @param[in]	xi		双精度实型数组，长度大于或等于（n+m-1），并且是2的整数次幂。输入数据的虚部；
* @param[in]	n		整型变量，输入数据长度（满足n=2^Exp）；
* @param[in]	m		整型变量，频率采样点数；
* @param[in]	f1		双精度实型变量，起始频率；
* @param[in]	f2		双精度实型变量，终止频率；
* @param[out]	xr		输出变换结果的实部；
* @param[out]	xi		输出变换结果的虚部；
*
* @return	无
*/
ISL_EXPORT void ISL_czt(double *xr, double *xi, int n, int m, double f1, double f2);



} /* End Of namespace ISLIB */
#endif /* PAI_FRAME_ISEISLIB_TRANSFORM_H_ */
