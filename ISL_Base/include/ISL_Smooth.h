/**
*	@file	ISL_Smooth.h 
*	@brief	[Header file of Smooth Functions], 平滑；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-02-10
*	@refer	SU CWP
*/	

#ifndef PAI_FRAME_ISEISLIB_SMOOTH_H
#define PAI_FRAME_ISEISLIB_SMOOTH_H

#include "ISL_SU_Matrix.h"
#include "ISL_Convolution.h"
#include "ISL_MathBase.h"
#include "ISL_Interpolation.h"

namespace ISLib {

class ISL_EXPORT Smooth
{

private:
	/**
	* @brief	此函数内部调用,加权平滑核心算法
	*
	**/
	static int ISL_smooth_weights(int k1, int k2, int pos, int mid, float *ss, float s0 = 0.5);

public:
	/*平滑滤波器*/
	/**
	* @brief	一维高斯平滑滤波器
	*
	* @param[in]	ns		输入数据个数
	* @param[in]	nsr		width (in samples) of the gaussian for which amplitude > 0.5*max amplitude
	* @param[in]	data	输入数据（个数为ns）
	*
	* @param[out]	data	返回平滑后的结果，数组个数为 [ns]
	*
	* @return	no
	*/
	static void ISL_gaussian1d_smoothing (int ns, int nsr, float *data);


	/**
	* @brief	高斯平滑直方图
	*
	* @param[in]	nintlh		输入数据个数
	* @param[in]	pdf			输入数据
	*
	* @param[out]	pdf			返回平滑后的结果，数组个数为 [nintlh]
	*
	* @return	no
	*/
	static void ISL_smooth_histogram (int nintlh, float *pdf);


	/**
	* @brief	一维滚动窗口平均平滑法(生成滤波器，滚动波形)
	*
	* @param[in]	flag	=1 for rectangular window.
	* 						=2 for triangular (weighted) window
	* @param[in]	nl		number of left (past) data points used
	* @param[in]	nr		number of right (future) data points used
	*
	* @param[out]	filter	返回平滑后的结果，数组个数为 [nl+nr+1]
	*
	* @return	no
	*/
	static void ISL_rwa_smoothing_filter (int flag, int nl, int nr, float *filter);

public:
	Smooth(){ ; }
	~Smooth(){ ; }

	/*平滑函数*/

	/**
	* @brief	9点平滑
	*
	* @param[in]	n		输入数组的个数
	* @param[in]	in[]	输入的数组
	*
	* @param[out]	out[]	返回平滑后的结果，数组个数为n
	*
	* @return	no
	*/
	static void ISL_smooth_9_3( int n, float in[], float out[] );


	/**
	* @brief	五点三次平滑
	*
	* @param[in]	n		输入数组的个数, 要求n>=5
	* @param[in]	in[]	输入的数组
	*
	* @param[out]	out[]	返回平滑后的结果，数组个数为n, 进入函数前请分配好内存空间
	*
	* @return	no
	*/
	static void ISL_smooth_5_3( int n, float in[], float out[] );


	/**
	* @brief	多点平滑
	*
	* @param[in]	vlc		输入的数组
	* @param[in]	ns		输入数组的个数
	* @param[in]	spn		平滑点数
	* @param[in]	spn		平滑次数
	*
	* @param[out]	vlc		返回平滑后的结果，数组个数为n
	*
	* @return	no
	*/
	static void ISL_smooth_MuliPoints(float *vlc, int ns, int spn, int cpn);


	/**
	* @brief	 高斯离散点一维平滑
	*
	* @param[in]	nin			输入数组的个数
	* @param[in]	xin[]		输入的X数组，样点对应的坐标
	* @param[in]	yin[]		输入的Y数组，输入数据的样点值，输入个数为nin
	* @param[in]	intr_num	输出数组的个数
	* @param[in]	intr_idx[]	输出的X数组，输出数据样点对应的坐标
	* @param[in]	c			平滑时窗
	*
	* @param[out]	yin[]		返回平滑后的结果，数组个数为intr_num
	*
	* @return	no
	*/
	static void ISL_smooth_disc(int nin, float xin[], float yin[], int intr_num, float intr_idx[], int c);


	/**
	* @brief	 加权平均平滑法
	*
	* @param[in]	vals		输入的数组
	* @param[in]	num			输入数组的个数
	* @param[in]	count		平滑计算周期，最长不能超过65，大于65将自动等于65
	* @param[in]	s0			调整系数
	*
	* @param[out]	vals		返回平滑后的结果，数组个数为num
	*
	* @return	返回计算的周期
	*/
	static int ISL_smooth_WA (float *vals, int num, int count, float s0 = 0.5);


	/**
	* @brief	 帽子平滑法
	*
	* @param[in]	num			输入/输出数组的个数
	* @param[in]	yin			输入的数组
	* @param[in]	step		至少是num的二分之一
	* @param[in]	flag		1 = 9点平滑
	* 							2 = 多点多次平滑
	*
	* @param[out]	yout		返回平滑后的结果，数组个数为num
	*
	* @return	返回计算的周期
	*/
	static void ISL_smooth_CAP(int num, float *yin, float *yout, int step=5, int flag=1);

};


} /*End of ISLib*/

#endif
