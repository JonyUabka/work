//*
//*	@file	ISL_Statistics.h
//*	@brief	[Header file of Statistics Functions], 统计；
//*	@see	ISeisLib Manual
//*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
//*	@date	2014-06-04
//*	@refer	SU CWP
//*/

//#include "ISL_Statistics.h"
//#include "ISL_Average.h"

//namespace ISLib {
///* 半对数-数据相关  */
//void ISL_log1(int n, double *x, double *y, double t, double *a)
//{
//	int i;
//	double xx, yy, dx, dxy;
//	xx = 0.0;
//	yy = 0.0;
//	for (i = 0; i <= n - 1; i++) {
//		xx = xx + x[i] / n;
//		yy = yy + log(y[i]) / log(t) / n;
//	}
//	dx = 0.0;
//	dxy = 0.0;

//	for (i = 0; i <= n - 1; i++) {
//		a[2] = x[i] - xx;
//		dx = dx + a[2] * a[2];
//		dxy = dxy + a[2] * (log(y[i]) / log(t) - yy);
//	}

//	a[1] = dxy / dx;
//	a[0] = yy - a[1] * xx;
//	a[0] = a[0] * log(t);
//	a[0] = exp(a[0]);
//	a[2] = 0.0;
//	a[6] = 0.0;
//	a[4] = 0.0;
//	a[5] = 1.0e+30;

//	for (i = 0; i <= n - 1; i++) {
//		a[3] = a[1] * x[i] * log(t);
//		a[3] = a[0] * exp(a[3]);
//		a[2] = a[2] + (y[i] - a[3]) * (y[i] - a[3]);
//		dx = fabs(y[i] - a[3]);
//		if (dx > a[4])
//			a[4] = dx;
//		if (dx < a[5])
//			a[5] = dx;
//		a[6] = a[6] + dx / n;
//	}

//	a[3] = sqrt(a[2] / n);

//	return;
//}

///* 对数-数据相关  */
//void ISL_log2(int n, double *x, double *y, double *a)
//{
//	int i;
//	double xx, yy, dx, dxy;

//	xx = 0.0;
//	yy = 0.0;
//	for (i = 0; i <= n - 1; i++) {
//		xx = xx + log(x[i]) / n;
//		yy = yy + log(y[i]) / n;
//	}

//	dx = 0.0;
//	dxy = 0.0;
//	for (i = 0; i <= n - 1; i++) {
//		a[2] = log(x[i]) - xx;
//		dx = dx + a[2] * a[2];
//		dxy = dxy + a[2] * (log(y[i]) - yy);
//	}

//	a[1] = dxy / dx;
//	a[0] = yy - a[1] * xx;
//	a[0] = exp(a[0]);
//	a[2] = 0.0;
//	a[6] = 0.0;
//	a[4] = 0.0;
//	a[5] = 1.0e+30;

//	for (i = 0; i <= n - 1; i++) {
//		a[3] = a[1] * log(x[i]);
//		a[3] = a[0] * exp(a[3]);
//		a[2] = a[2] + (y[i] - a[3]) * (y[i] - a[3]);
//		dx = fabs(y[i] - a[3]);
//		if (dx > a[4])
//			a[4] = dx;
//		if (dx < a[5])
//			a[5] = dx;
//		a[6] = a[6] + dx / n;
//	}

//	a[3] = sqrt(a[2] / n);

//	return;
//}


///*中位数*/
//int ISL_midValue(float *in, int num, float &midValue)
//{
//	int j, mid = -1;
//	int half = num / 2;
//	float tmp;

//	/////希尔排序，从小到大/////
//	while (half > 0) {
//		for (int i = half; i < num; i++) {
//			tmp = in[i];
//			j = i - half;
//			while (j >= 0 && in[j] > tmp) {
//				in[j + half] = in[j];
//				j -= half;
//			}
//			in[j + half] = tmp;
//		}
//		half /= 2;
//	}

//	mid = num / 2;
//	midValue = in[mid];

//	return mid;
//}

//} /* End Of namespace ISLIB
