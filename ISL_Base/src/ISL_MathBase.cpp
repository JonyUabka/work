//*
// *	@file	ISL_MathBase.h
// *	@brief	[Header file of Base Mathematics Functions], 基础数学函数；
// *	@see	ISeisLib Manual
// *	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
// *	@date	2014-02-26
// *	@refer	SU CWP
// */

//#include "ISL_MathBase.h"

//namespace ISLib
//{

///* === 正态分布随机数 === */
//double ISL_normalDistributionRandom()
//{
//	static double V1 = 0, V2 = 0, S = 0;
//	static int phase = 0;
//	double X = 0;

//	if (phase == 0) {
//		do {
//			double U1 = (double) rand() / RAND_MAX;
//			double U2 = (double) rand() / RAND_MAX;

//			V1 = 2 * U1 - 1;
//			V2 = 2 * U2 - 1;
//			S = V1 * V1 + V2 * V2;
//		} while (S >= 1 || S == 0);

//		X = V1 * sqrt(-2 * log(S) / S);
//	} else
//		X = V2 * sqrt(-2 * log(S) / S);

//	phase = 1 - phase;

//	return X;
//}

///* === 浮点误差的处理 === */
//int ISL_dblcmp(float d)
//{
//	if (fabs(d) < EP)
//		return 0;
//	return (d > 0) ? 1 : -1;
//}

///* === 判断是否为奇数，否则是偶数 === */
//bool ISL_isOdd(int n)
//{
//	return (n % 2 != 0);
//}

///**********************************************************************************************************************
// *
// *  功能：return float value of sinc(x) for x input as a float
// *
// *  说明：
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		float 				x					In			value at which to evaluate sinc(x)
// *
// *  返回：sinc(x)
// *
// **********************************************************************************************************************/
//float ISL_fsinc(float x)
//{
//	float pix;

//	if (x == 0.0)
//	{
//		return 1.0;
//	}
//	else
//	{
//		pix = PI * x;
//		return sin(pix) / pix;
//	}
//}

///**********************************************************************************************************************
// *
// *  功能：return double precision sinc(x) for double precision x
// *
// *  说明：
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		double 				x					In			value at which to evaluate sinc(x)
// *
// *  返回：sinc(x)
// *
// **********************************************************************************************************************/
//double ISL_dsinc(double x)
//{
//	double pix;

//	if (x == 0.0)
//	{
//		return 1.0;
//	}
//	else
//	{
//		pix = PI * x;
//		return sin(pix) / pix;
//	}
//}

//}/*End of ISLib
