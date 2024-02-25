/**
*	@file	ISL_Convolution.h 
*	@brief	[Header file of Convolution Functions], 卷积函数；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-02-26
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_CONVOLUTION_H
#define PAI_FRAME_ISEISLIB_CONVOLUTION_H

#include "ISL_UserDefine.h"

namespace ISLib {

	/**
	* @brief	互相关函数
	*			Compute z = x 与 y 的相关;
	*
	*			ifx + lx - 1
	*		 z[i] = sum x[j]*y[i+j]; i = ifz, ... , ifz + lz - 1
	*			j = ifx
	*			
	* @param[in]	lx		length of x array
	* @param[in]	lfx		sample index of first x
	* @param[in]	*x		array[lx] to be cross-correlated with y	
	* @param[in]	ly		length of y array	
	* @param[in]	ify		sample index of first y
	* @param[in]	*y		array[ly] with which x is to be cross-correlated	
	* @param[in]	lz		length of z array	
	* @param[in]	ifz		sample index of first z	
	* @param[out]	*z		array[lz] containing x cross-correlated with y	
	*
	* @return	no
	*/
	ISL_EXPORT void ISL_xcor(	int lx,	int ifx, float *x, 
					int ly, int ify, float *y,
					int lz, int ifz, float *z);

	
	/**
	* @brief	卷积函数[ISL_conv]的计算核心函数（short X 的优化算法）
	*			外部无需调用！！
	*/
	ISL_EXPORT void ISL_convs(	int lx, int ifx, float *x,
					int ly, int ify, float *y, 
					int lz, int ifz, float *z);


	/**
	* @brief	卷积函数
	*			Compute z = x 与 y 的卷积;
	*
	*			ifx + lx - 1
	*	     z[i] = sum x[j]*y[i+j]; i = ifz, ... , ifz + lz - 1
	*			j = ifx
	*			
	* @param[in]	lx		length of x array
	* @param[in]	lfx		sample index of first x
	* @param[in]	*x		array[lx] to be convolved with y	
	* @param[in]	ly		length of y array	
	* @param[in]	ify		sample index of first y
	* @param[in]	*y		array[ly] with which x is to be convolved	
	* @param[in]	lz		length of z array	
	* @param[in]	ifz		sample index of first z	
	* @param[out]	*z		array[lz] containing x convolved with y	
	*
	* @return	no
	*/
	ISL_EXPORT void ISL_conv(	int lx, int ifx, float *x, 
					int ly, int ify, float *y,
					int lz, int ifz, float *z);



	/**
	* @brief	用FFT计算两个有限长序列的线性卷积
	*
	* @param[in]	x[len]		双精度实型一维数组，长度为len，实际输入数据序列个数为m，输出线性卷积结果。
	* @param[in]	y[n]		双精度实型一维数组，长度为len，实际输入数据序列个数为n。
	* @param[in]	m			序列x(i)的长度。
	* @param[in]	n			序列y(i)的长度
	* @param[out]	len			整型变量，序列x(i)和y(i)线性卷积的长度，len≥m+n-1，且必须是2的整数次幂。
	*
	* @return	no
	*/
	ISL_EXPORT void ISL_conv(double *x, double *y, int m, int n, int len);


	/**
	* @brief	计算合成记录的简单卷积， 输入两个任意长度的离散点序列，返回一个卷积后的结果
	*
	* @param[in]	r		输入的反射系数序列
	* @param[in]	nr		反射系数的维数
	* @param[in]	w		输入的地震子波序列
	* @param[in]	nw		地震子波的维数
	* @param[out]	syn		合成的地震记录
	*
	* @return	no
	*/
	ISL_EXPORT void ISL_synConv(float *r, int nr, float *w, int nw, float *syn);





}/*End of ISLib*/

#endif
