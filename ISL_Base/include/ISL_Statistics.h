/**
*	@file	ISL_Statistics.h
*	@brief	[Header file of Statistics Functions], 统计；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-06-04
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_STATISTICS_H_
#define PAI_FRAME_ISEISLIB_STATISTICS_H_

#include "ISL_UserDefine.h"
#include "ISL_Average.h"

namespace ISLib {

/**
 * @brief	半对数-数据相关
 *
 * @param[in]	n			数据点数；
 * @param[in]	x[n]		存放n个数据点。要求所有的y>0；
 * @param[in]	y[n]		存放n个数据点。要求所有的y>0
 * @param[in]	t			指数函数的底。要求t>0；
 * @param[out]	a[7]		返回拟合函数的参数及各种统计量。其意义如下：
 *							a(0)：拟合函数y=bt^(ax)中的b
 *							a(1)：拟合函数y=bt^(ax)中的a
 *							a(2)：偏差平方和q
 *							a(3)：平均标准差s
 *							a(4)：最大偏差umax
 *							a(5)：最小偏差umin
 *							a(6)：偏差平均值u
 *
 * @return	无
 */
ISL_EXPORT void ISL_log1(int n, double *x, double *y, double t, double *a);


/**
 * @brief	对数-数据相关
 *
 * @param[in]	n			数据点数；
 * @param[in]	x[n]		存放n个数据点。要求所有的y>0；
 * @param[in]	y[n]		存放n个数据点。要求所有的y>0
 * @param[out]	a[7]		返回拟合函数的参数及各种统计量。其意义如下：
 *							a(0)：拟合函数y=bt^(ax)中的b
 *							a(1)：拟合函数y=bt^(ax)中的a
 *							a(2)：偏差平方和q
 *							a(3)：平均标准差s,s = sqrt(q/n)
 *							a(4)：最大偏差umax
 *							a(5)：最小偏差umin
 *							a(6)：偏差平均值u
 *
 * @return	无
 */
ISL_EXPORT void ISL_log2(int n, double *x, double *y, double *a);



/**
* @brief	中位数
*
* @param[in]	inData[num]		输入数据；
* @param[in]	num				输入数据的大小；
* @param[out]	midValue		输出计算得到的中值;
*
* @return	中位值在排序后数据中的位置
*/
ISL_EXPORT int ISL_midValue(float *in, int num, float &midValue);



/**
* @brief	方差
*
* @param[in]	in			输入的序列；
* @param[in]	num			输入序列的个数；
*
* @return	返回方差的值（double）
*/
template<typename T>
double ISL_varianceError(T *in, long num)
{
	double ave = ISL_arithmeticAve(in, num);

	double err = 0;
	for(long i=0; i<num; i++){
		err = err + ((double)in[i] - ave) * ((double)in[i] - ave);
	}
	return (err/(double)num);
}



/**
* @brief	均方差（标准差（Standard Deviation））
*
* @param[in]	in			输入的序列；
* @param[in]	num			输入序列的个数；
*
* @return	返回均方差的值（double）
*/
template<typename T>
double ISL_RMSError(T *in, long num)
{
	double err = ISL_varianceError(in, num);
	return sqrt(err);
}


} /* End Of namespace ISLIB */
#endif /* PAI_FRAME_ISEISLIB_STATISTICS_H_ */
