/**
*	@file	ISL_NumericalIntegration.cpp
*	@brief	[Source file of Numerical Integration Functions], 数值积分；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-06-24
*	@refer
*/


#include "ISL_NumericalIntegration.h"


namespace ISLib {

double ISL_ffftsf(double x)
{
	double y;
	y = exp(-x * x);
	return (y);
}

/* 变步长梯形求积法 */
double ISL_fffts(double a, double b, double eps)
{
	int n, k;
	double fa, fb, h, t1, p, s, x, t = 0;
	fa = ISL_ffftsf(a);
	fb = ISL_ffftsf(b);
	n = 1;
	h = b - a;
	t1 = h * (fa + fb) / 2.0;
	p = eps + 1.0;
	while (p >= eps) {
		s = 0.0;
		for (k = 0; k <= n - 1; k++) {
			x = a + (k + 0.5) * h;
			s = s + ISL_ffftsf(x);
		}
		t = (t1 + h * s) / 2.0;
		p = fabs(t1 - t);
		t1 = t;
		n = n + n;
		h = h / 2.0;
	}
	return (t);
}

/* 变步长辛卜森求积法 */
double ISL_fsimpf(double x)
{
	double y;
	y = log(1.0 + x) / (1.0 + x * x);
	return (y);
}

double ISL_fsimp(double a, double b, double eps)
{
	int n, k;
	double h, t1, t2, s1, s2 = 0, ep, p, x;
	n = 1;
	h = b - a;
	t1 = h * (ISL_fsimpf(a) + ISL_fsimpf(b)) / 2.0;
	s1 = t1;
	ep = eps + 1.0;
	while (ep >= eps) {
		p = 0.0;
		for (k = 0; k <= n - 1; k++) {
			x = a + (k + 0.5) * h;
			p = p + ISL_fsimpf(x);
		}
		t2 = (t1 + h * p) / 2.0;
		s2 = (4.0 * t2 - t1) / 3.0;
		ep = fabs(s2 - s1);
		t1 = t2;
		s1 = s2;
		n = n + n;
		h = h / 2.0;
	}
	return (s2);
}





} /* End Of namespace ISLIB */
