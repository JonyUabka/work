///**
// *	@file	ISL_MathBase.h
// *	@brief	[Header file of Base Mathematics Functions], 基础数学函数；
// *	@see	ISeisLib Manual
// *	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
// *	@date	2014-02-26
// *	@refer	SU CWP
// */

//#ifndef PAI_FRAME_ISEISLIB__MATHBASE_H
//#define PAI_FRAME_ISEISLIB__MATHBASE_H

//#include "ISL_UserDefine.h"

//namespace ISLib
//{

///**
// * @brief	判断浮点数相等
// * @param[in]	f1	浮点数（float or double）
// * @param[in]	f2	浮点数（float or double）
// * @return	正确：1
// *			错误：0
// */
//template<typename T>
//int ISL_floatEqual(T f1, T f2)
//{
//	double dd = fabs((double) (f2 - f1));
//	if (dd < MinDataError)
//		return 1;
//	return 0;
//}

///**
// * @brief	浮点数交换，val1 与 val2的值互换
// * @param[in]	val1	浮点数（float or double）
// * @param[in]	val2	浮点数（float or double）
// * @return	no
// */
//template<typename T>
//void ISL_swapValue(T &val1, T &val2)
//{
//	T vv = val1;
//	val1 = val2;
//	val2 = vv;
//}

///**
// * @brief	二维维数组的行列互换
// *			in与out的大小都是 w * h
// *			flag = 1， 行转列
// *			flag = 2， 列转行
// *
// * @param[in]	in[]	输入数据
// * @param[in]	w	列数 columns
// * @param[in]	h	行数 rows
// * @param[out]	out[]	输出数据
// * @param[in]	flag	标志位，，默认值为1
// *
// * @return	no
// */
//template<typename T>
//void ISL_swapValue(T in[], int w, int h, T out[], int flag = 1)
//{
//	if (flag == 1) {
//		for (int i = 0; i < h; i++) {
//			for (int j = 0; j < w; j++) {
//				out[i + j * h] = in[j + i * w];
//			}
//		}
//	} else {
//		for (int i = 0; i < w; i++) {
//			for (int j = 0; j < h; j++) {
//				out[i + j * w] = in[j + i * h];
//			}
//		}
//	}
//}

//// ========================================================================================
//// === 					排序算法     															===
//// ========================================================================================
///**
// * @brief	冒泡排序（稳定型），排序，适用浮点一维数组
// *
// * @param[in]	val			输入数组
// * @param[in]	num			输入的数组的大小
// * @param[in]	flag		输入的标志位，排序方式；0=从小到大； 1=从大到小；
// *
// * @param[out]	val			输出数组（排序后）
// *
// * @return	no
// */
//template<typename T>
//void ISL_bubbleSort(T * val, long num, int flag = 0)/*flag=0,从小到大； flag=1，从大到小*/
//{
//	if (val == NULL || num <= 0)
//		return;

//	if (flag == 0) {
//		for (long m = num - 1; m > 0; m--) {
//			T temp;
//			for (long n = 0; n < m; n++) {
//				if (val[n] > val[n + 1]) {
//					temp = val[n];
//					val[n] = val[n + 1];
//					val[n + 1] = temp;
//				}
//			}
//		}
//	} else {
//		for (long m = num - 1; m > 0; m--) {
//			T temp;
//			for (long n = 0; n < m; n++) {
//				if (val[n] < val[n + 1]) {
//					temp = val[n];
//					val[n] = val[n + 1];
//					val[n + 1] = temp;
//                }
//            }
//        }
//    }
//}

///**
// * @brief	插入排序（稳定型），排序，适用浮点一维数组
// *
// * @param[in]	val			输入数组
// * @param[in]	num			输入的数组的大小
// * @param[in]	flag		输入的标志位，排序方式；0=从小到大； 1=从大到小；
// *
// * @param[out]	val			输出数组（排序后）
// *
// * @return	no
// */
//template<typename T>
//void ISL_insertSort(T * val, long num, int flag = 0)/*flag=0,从小到大； flag=1，从大到小*/
//{
//	if (val == NULL || num <= 0)
//		return;

//	long i = 0, j = 0;
//	T temp = 0;

//	if (flag == 0) {
//		for (i = 1; i < num; i++) {
//			j = i;
//			temp = val[i];
//			while (j > 0 && temp < val[j - 1]) {
//				val[j] = val[j - 1];
//				j--;
//			}
//			val[j] = temp;
//		}
//	} else {
//		for (i = 1; i < num; i++) {
//			j = i;
//			temp = val[i];
//			while (j > 0 && temp > val[j - 1]) {
//				val[j] = val[j - 1];
//				j--;
//			}
//			val[j] = temp;
//		}
//	}
//}

///**
// * @brief	选择排序（不稳定型），排序，适用浮点一维数组
// *
// * @param[in]	val			输入数组
// * @param[in]	num			输入的数组的大小
// * @param[in]	flag		输入的标志位，排序方式；0=从小到大； 1=从大到小；
// *
// * @param[out]	val			输出数组（排序后）
// *
// * @return	no
// */
//template<typename T>
//void ISL_selectSort(T * val, long num, int flag = 0)/*flag=0,从小到大； flag=1，从大到小*/
//{
//	if (val == NULL || num <= 0)
//		return;

//	if (flag == 0) {
//		for (long i = 0; i < num; i++) {
//			long pos = i;
//			for (long j = i + 1; j < num; j++) {
//				if (val[pos] > val[j]) {
//					pos = j;
//				}
//			}
//			T temp;
//			temp = val[pos];
//			val[pos] = val[i];
//			val[i] = temp;
//		}
//	} else {
//		for (long i = 0; i < num; i++) {
//			long pos = i;
//			for (long j = i + 1; j < num; j++) {
//				if (val[pos] < val[j]) {
//					pos = j;
//				}
//			}
//			T temp;
//			temp = val[pos];
//			val[pos] = val[i];
//			val[i] = temp;
//		}
//	}
//}

//// ========================================================================================
//// === 					数据查找     															===
//// ========================================================================================

///**
// * @brief	查找数组中对应的值，存在则返回数组下标，-1表示出错
// *
// * @param[in]	val			输入数组
// * @param[in]	num			输入的数组的大小
// * @param[in]	a			查找的数值；
// *
// * @return	返回-1为出错，返回数组下标为正确。
// */
//template<typename T>
//int ISL_findValue(T * val, long num, T a)
//{
//	if (val == NULL || num <= 0)
//		return -1;

//	for (long i = 0; i < num; i++) {
//		if (val[i] == a) {
//			return i;
//		}
//		if (ISL_floatEqual(val[i], a) == 1) {
//			return i;
//		}
//	}
//	return -1;
//}


///**
// * @brief	绝对误差
// *
// * @param[in]	c1			输入数组1
// * @param[in]	c2			输入数组2
// * @param[in]	num			数组大小
// * @param[out]	val			返回绝对误差数组（个数为num）
// *
// * @return	无。
// */
//template <typename T>
//void ISL_absoluteError(T * c1, T * c2, int num, T *& val)
//{
//	if(c1 == NULL || c2 == NULL || num<=0) return;

//	if(val){ delete []val; val = NULL; }
//	val = new T[num];
//	if(val == NULL) return ;

//	for(int i=0; i<num; i++){
//		val[i] = fabs(c1[i] - c2[i]);
//	}
//}


///**
// * @brief	 相对误差
// *
// * @param[in]	c1			输入数组1
// * @param[in]	c2			输入数组2
// * @param[in]	num			数组大小
// * @param[out]	val			返回绝对误差数组（个数为num）
// *
// * @return	无。
// */
//template  <typename T>
//void ISL_relativeError(T * c1, T * c2, int num, T *& val)
//{
//	if(c1 == NULL || c2 == NULL || num<=0) return ;

//	if(val){ delete []val; val = NULL; }
//	val = new T[num];
//	if(val == NULL) return;

//	for(int i=0; i<num; i++){
//		val[i] = fabs (fabs(c1[i] - c2[i]) / c1[i]);
//	}
//}


///**
// * @brief	 相关关系
// *
// * @param[in]	c1			输入数组1
// * @param[in]	c2			输入数组2
// * @param[in]	num			数组大小
// * @param[out]	val			返回绝对误差数组（个数为num）
// *
// * @return	无。
// */
//template  <typename T>
//void ISL_correlativity(T * c1, T * c2, int num, T & val)
//{
//	if(c1 == NULL || c2 == NULL || num<=0){
//		val = -9999;
//		return;
//	}

//	T c1_ave, c2_ave;
//	for(int i=0; i<num; i++){
//		c1_ave = c1_ave + c1[i];
//		c2_ave = c2_ave + c2[i];
//	}
//	c1_ave = c1_ave/(T)num;
//	c2_ave = c2_ave/(T)num;

//	T v1, v2, v3;

//	for(int i=0; i<num; i++){
//		v1 = v1 + (c1[i]-c1_ave)*(c2[i]-c2_ave);
//		v2 = v2 + (c1[i]-c1_ave)*(c1[i]-c1_ave);
//		v3 = v3 + (c2[i]-c2_ave)*(c2[i]-c2_ave);
//	}

//	val = v1/(sqrt(v2)*sqrt(v3));
//}


///**
// * @brief	查找数组中的最大值、最小值、极值
// *
// * @param[in]	val			输入数组
// * @param[in]	num			输入的数组的大小
// * @param[out]	max			返回数组中的最大值
// * @param[out]	min			返回数组中的最小值
// * @param[out]	top			返回数值中的极值（绝对值中的最大值）
// *
// * @return	无。
// */
//template<typename T>
//void ISL_findMaxValue(T *val, long num, T &max, T &min, T &top)
//{
//	if (val == NULL || num <= 0)
//	{
//		return;
//	}

//	max = val[0];
//	min = val[0];
//	top = val[0];

//	for (long i = 0; i < num; i++)
//	{
//		if (val[i] > max)
//			max = val[i];
//		if (val[i] < min)
//			min = val[i];
//	}
//	if (fabs(max) > fabs(min))
//		top = fabs(max);
//	else
//		top = fabs(min);
//}


///**
// * @brief	查找数组中的最大值及位置
// *
// * @param[in]	val			输入数组
// * @param[in]	num			输入的数组的大小
// * @param[out]	max			返回数组中的最大值
// * @param[out]	pos			返回数组中的最大值在数组中的位置
// *
// * @return	无。
// */
//template<typename T>
//void ISL_locateMax(T *val, long num, T &max, long &pos)
//{
//	if (val == NULL || num <= 0) {
//		return;
//	}

//	max = val[0];
//	for (long i = 0; i < num; i++) {
//		if (val[i] > max) {
//			max = val[i];
//			pos = i;
//		}
//	}
//}


///**
// * @brief	查找数组中的最小值及位置
// *
// * @param[in]	val			输入数组
// * @param[in]	num			输入的数组的大小
// * @param[out]	min			返回数组中的最小值
// * @param[out]	pos			返回数组中的最小值在数组中的位置
// *
// * @return	无。
// */
//template<typename T>
//void ISL_locateMin(T *val, long num, T &min, long &pos)
//{
//	if (val == NULL || num <= 0) {
//		return;
//	}

//	min = val[0];
//	for (long i = 0; i < num; i++) {
//		if (val[i] < min) {
//			min = val[i];
//			pos = i;
//		}
//	}
//}


///**
// * @brief	查找数组中的极值及位置
// *
// * @param[in]	val			输入数组
// * @param[in]	num			输入的数组的大小
// * @param[out]	top			返回数组中的极值
// * @param[out]	pos			返回数组中的极值在数组中的位置
// *
// * @return	无。
// */
//template<typename T>
//void ISL_locateTop(T *val, long num, T &top, long &pos)
//{
//	if (val == NULL || num <= 0) {
//		return;
//	}

//	T max = val[0];
//	long maxPos = 0;

//	T min = val[0];
//	long minPos = 0;

//	top = val[0];

//	for (long i = 0; i < num; i++) {
//		if (val[i] > max){
//			max = val[i];
//			maxPos = i;
//		}
//		if (val[i] < min){
//			min = val[i];
//			minPos = i;
//		}
//	}
//	if (fabs(max) > fabs(min)){
//		top = fabs(max);
//		pos = maxPos;
//	}else{
//		top = fabs(min);
//		pos = minPos;
//	}
//}

///**
// * @brief	正态分布随机数
// *
// * @return	double 正态分布随机数
// */
//ISL_EXPORT double ISL_normalDistributionRandom();

///**
// * @brief	浮点误差的处理
// * @param[in]	d	输入误差
// * @return	误差标准之内：0
// *			误差标准之外：d>0 返回 1
// *						  d<0 返回 -1
// */
//ISL_EXPORT int ISL_dblcmp(float d);

///**
// * @brief	判断是否为奇数
// * @param[in]	n	输入整型数值
// * @return	true ：为奇数
// *			false ：为偶数
// */
//ISL_EXPORT bool ISL_isOdd(int n);

///**
// * @brief	return float value of sinc(x) for x input as a float
// *			sinc(x) = sin(PI*x)/(PI*x) (float version)
// *
// * @param[in]	x[float]	value at which to evaluate sinc(x)
// * @return	sinc(x)
// */
//ISL_EXPORT float ISL_fsinc(float x);

///**
// * @brief	return double precision sinc(x) for double precision x
// *			sinc(x) = sin(PI*x)/(PI*x) (double version)
// *
// * @param[in]	x[double]	value at which to evaluate sinc(x)
// * @return	sinc(x)
// */
//ISL_EXPORT double ISL_dsinc(double x);

//} /*End of ISLib*/

//#endif
