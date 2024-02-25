/**
*	@file	ISL_Average.h
*	@brief	[Header file of Average Functions], 平均算法；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-06-24
*/

#ifndef PAI_FRAME_ISEISLIB_AVERAGE_H_
#define PAI_FRAME_ISEISLIB_AVERAGE_H_


#include "ISL_UserDefine.h"

namespace ISLib {

/**
* @brief	算术平均
*
* @param[in]	in			输入的序列；
* @param[in]	num			输入序列的个数；
*
* @return	平均值（double）
*/
template<typename T>
double ISL_arithmeticAve(T *in, long num)
{
	double ave = 0;
	for(long i=0; i<num; i++){
		ave = ave + (double)in[i];
	}
	ave = ave / (double)num;
	return ave;
}



/**
* @brief	均方根平均
*
* @param[in]	in			输入的序列；
* @param[in]	num			输入序列的个数；
*
* @return	平均值（double）
*/
template<typename T>
double ISL_RMSAve(T *in, long num)
{
	double ave = 0;
	for(long i=0; i<num; i++){
		ave = ave + (double)in[i]*(double)in[i];
	}
	ave = sqrt(ave / (double) num);
	return ave;
}



/**
* @brief	几何平均（输入序列必须是正数）
*
* @param[in]	in			输入的序列(序列样点必须>=0)；
* @param[in]	num			输入序列的个数；
*
* @return	平均值（double）
*/
template<typename T>
double ISL_geometricAve(T *in, long num)
{
	double ave = 1;
	for(long i=0; i<num; i++){
		ave = ave * (double)in[i];
	}
	ave = pow(ave, (double)1/(double)num);
	return ave;
}



/**
* @brief	加权平均
*
* @param[in]	in[num]		输入的序列(序列样点必须>=0)；
* @param[in]	inW[num]	输入的加权序列(序列样点必须>=0)；
* @param[in]	num			输入序列的个数；
*
* @return	平均值（double）
*/
template<typename T>
double ISL_weightedAve(T *in, T *inW, long num)
{
	double ave = 0;
	double w = 0;
	for(long i=0; i<num; i++){
		ave = ave + (double)in[i] * (double)inW[i];
		w = w + (double)inW[i];
	}
	ave = ave / w;
	return ave;
}



/**
* @brief	对数平均（输入序列必须是正数）
*
* @param[in]	in			输入的序列(序列样点必须>=0)；
* @param[in]	num			输入序列的个数；
*
* @return	平均值（double）
*/
template<typename T>
double ISL_logarithmicAve(T *in, long num)
{
	double ave = (double)in[0];
	double ds = log((double)in[0]);

	for(long i=1; i<num; i++){
		ave = ave - (double)in[i];
		ds = ds - log((double)in[i]);
	}
	ave = ave / ds;
	return ave;
}


/**
* @brief	调和平均（输入序列必须是正数）
*
* @param[in]	in			输入的序列(序列样点必须>=0)；
* @param[in]	num			输入序列的个数；
*
* @return	平均值（double）
*/
template<typename T>
double ISL_harmonicAve(T *in, long num)
{
	double ave = 0;
	for(long i=0; i<num; i++){
		ave = ave + (double)1.0/(double)in[i];
	}
	ave = ave / (double)num;
	ave = (double)1.0 / ave;
	return ave;
}




} /* End Of namespace ISLIB */
#endif /* PAI_FRAME_ISEISLIB_AVERAGE_H_ */
