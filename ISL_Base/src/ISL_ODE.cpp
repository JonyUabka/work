/**
*	@file	ISL_ODE.cpp
*	@brief	[Source file of ODE Functions], 常微分方程数值解；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-08-18
*	@refer
*/

#include "ISL_ODE.h"
#include "ISL_Alloc.h"

namespace ISLib {

// ========== 用改进的Euler法对一阶微分方程组进行定步长全区间积分 =================
void ISL_gelr1::run(double t, double *y, int n, double h, int k, double *z)
{
	double x, *d = NULL;
	d = alloc1double(n);

	for (int i = 0; i <= n - 1; i++)
		z[i * k] = y[i];

	for (int j = 1; j <= k - 1; j++) {
		x = t + (j - 1) * h;
		gelr1f(x, y, n, d);
		for (int i = 0; i <= n - 1; i++)
			y[i] = z[i * k + j - 1] + h * d[i];

		x = t + j * h;
		gelr1f(x, y, n, d);
		for (int i = 0; i <= n - 1; i++)
			d[i] = z[i * k + j - 1] + h * d[i];

		for (int i = 0; i <= n - 1; i++) {
			y[i] = (y[i] + d[i]) / 2.0;
			z[i * k + j] = y[i];
		}
	}
	free(d);
}



// ========== 一阶微分方程的变步长墨森法  =================
void ISL_gmrsn::run(double t, double h, int n, double *y, double eps, int k, double *z)
{
	int i, j, m, nn;
	double aa, bb, x, hh, p, dt, t0, qq, *a, *b, *c, *d, *u, *v;
	a = alloc1double(n);
	b = alloc1double(n);
	c = alloc1double(n);
	d = alloc1double(n);
	u = alloc1double(n);
	v = alloc1double(n);
	aa = t;
	for (i = 0; i <= n - 1; i++)
		z[i * k] = y[i];

	for (i = 1; i <= k - 1; i++) {
		x = aa + (i - 1) * h;
		nn = 1;
		hh = h;
		for (j = 0; j <= n - 1; j++)
			u[j] = y[j];
		p = 1.0 + eps;
		while (p >= eps) {
			for (j = 0; j <= n - 1; j++) {
				v[j] = y[j];
				y[j] = u[j];
			}
			dt = h / nn;
			t = x;
			for (m = 0; m <= nn - 1; m++) {
				gmrsnf(t, y, n, d);
				for (j = 0; j <= n - 1; j++) {
					a[j] = d[j];
					y[j] = y[j] + hh * d[j] / 3.0;
				}
				t0 = t + hh / 3.0;
				gmrsnf(t0, y, n, d);
				for (j = 0; j <= n - 1; j++) {
					b[j] = d[j];
					y[j] = y[j] + hh * (d[j] - a[j]) / 6.0;
				}
				gmrsnf(t0, y, n, d);
				for (j = 0; j <= n - 1; j++) {
					b[j] = d[j];
					bb = (d[j] - 4.0 * (b[j] + a[j] / 4.0) / 9.0) / 8.0;
					y[j] = y[j] + 3.0 * hh * bb;
				}
				t0 = t + hh / 2.0;
				gmrsnf(t0, y, n, d);
				for (j = 0; j <= n - 1; j++) {
					c[j] = d[j];
					qq = d[j] - 15.0 * (b[j] - a[j] / 5.0) / 16.0;
					y[j] = y[j] + 2.0 * hh * qq;
				}
				t0 = t + hh;
				gmrsnf(t0, y, n, d);
				for (j = 0; j <= n - 1; j++) {
					qq = c[j] - 9.0 * (b[j] - 2.0 * a[j] / 9.0) / 8.0;
					qq = d[j] - 8.0 * qq;
					y[j] = y[j] + hh * qq / 6.0;
				}
				t = t + dt;
			}
			p = 0.0;
			for (j = 0; j <= n - 1; j++) {
				qq = fabs(y[j] - v[j]);
				if (qq > p)
					p = qq;
			}
			hh = hh / 2.0;
			nn = nn + nn;
		}
		for (j = 0; j <= n - 1; j++)
			z[j * k + i] = y[j];
	}

	free(a);
	free(b);
	free(c);
	free(d);
	free(u);
	free(v);

	return;
}




// ========== 一阶微分方程的定步长四阶龙格库塔法  =================
void ISL_grkt1::run(double t, double *y, int n, double h, int k, double *z)
{
	int i, j, l;
	double a[4], tt, *b, *d;
	b = alloc1double(n);
	d = alloc1double(n);
	a[0] = h / 2.0;
	a[1] = a[0];
	a[2] = h;
	a[3] = h;

	for (i = 0; i <= n - 1; i++)
		z[i * k] = y[i];

	for (l = 1; l <= k - 1; l++) {
		grkt1f(t, y, n, d);
		for (i = 0; i <= n - 1; i++)
			b[i] = y[i];
		for (j = 0; j <= 2; j++) {
			for (i = 0; i <= n - 1; i++) {
				y[i] = z[i * k + l - 1] + a[j] * d[i];
				b[i] = b[i] + a[j + 1] * d[i] / 3.0;
			}
			tt = t + a[j];
			grkt1f(tt, y, n, d);
		}
		for (i = 0; i <= n - 1; i++)
			y[i] = b[i] + h * d[i] / 6.0;
		for (i = 0; i <= n - 1; i++)
			z[i * k + l] = y[i];
		t = t + h;
	}

	free(b);
	free(d);

	return;
}


}/* End Of namespace ISLIB */


