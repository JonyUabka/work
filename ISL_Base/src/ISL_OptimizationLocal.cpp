/**
*	@file	ISL_OptimizationLocal.cpp
*	@brief	[Source file of ISL_OptimizationLocal], 局部最优化；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-07-1
*	@refer
*/

#include "ISL_OptimizationLocal.h"
#include "ISL_Alloc.h"

namespace ISLib {

// ========== 求非线性方程组的梯度法 =================
int ISL_dsnse::run(double *x, int n, double eps, int js)
{
	int l, j;
	double f, d, s, *y = NULL;
	y = alloc1double(n);
	l = js;
	f = dsnsef(x, y, n);
//	cout<<"dsnsef = "<<f<<endl;

	while (f >= eps) {
		l = l - 1;
		if (l == 0) {
			free(y);
			return (js);
		}
		d = 0.0;
		for (j = 0; j <= n - 1; j++)
			d = d + y[j] * y[j];
		if (d + 1.0 == 1.0) {
			free(y);
			return (-1);
		}
		s = f / d;
		for (j = 0; j <= n - 1; j++)
			x[j] = x[j] - s * y[j];
		f = dsnsef(x, y, n);
	}
	free(y);
	return (js - l);
}



// ========== 求非线性方程组的拟牛顿法 =================
int ISL_dnetn::agaus(double *a, double *b, int n)
{
	int *js = NULL, l, k, i, j, is = 0, p, q;
	double d, t;
	js = alloc1int(n);
	l = 1;
	for (k = 0; k <= n - 2; k++) {
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++) {
				t = fabs(a[i * n + j]);
				if (t > d) {
					d = t;
					js[k] = j;
					is = i;
				}
			}
		if (d + 1.0 == 1.0)
			l = 0;
		else {
			if (js[k] != k)
				for (i = 0; i <= n - 1; i++) {
					p = i * n + k;
					q = i * n + js[k];
					t = a[p];
					a[p] = a[q];
					a[q] = t;
				}
			if (is != k) {
				for (j = k; j <= n - 1; j++) {
					p = k * n + j;
					q = is * n + j;
					t = a[p];
					a[p] = a[q];
					a[q] = t;
				}
				t = b[k];
				b[k] = b[is];
				b[is] = t;
			}
		}
		if (l == 0) {
			free(js);
//			printf("fail\n");
			return (0);
		}
		d = a[k * n + k];
		for (j = k + 1; j <= n - 1; j++) {
			p = k * n + j;
			a[p] = a[p] / d;
		}
		b[k] = b[k] / d;
		for (i = k + 1; i <= n - 1; i++) {
			for (j = k + 1; j <= n - 1; j++) {
				p = i * n + j;
				a[p] = a[p] - a[i * n + k] * a[k * n + j];
			}
			b[i] = b[i] - a[i * n + k] * b[k];
		}
	}
	d = a[(n - 1) * n + n - 1];
	if (fabs(d) + 1.0 == 1.0) {
		free(js);
//		printf("fail\n");
		return (0);
	}
	b[n - 1] = b[n - 1] / d;
	for (i = n - 2; i >= 0; i--) {
		t = 0.0;
		for (j = i + 1; j <= n - 1; j++)
			t = t + a[i * n + j] * b[j];
		b[i] = b[i] - t;
	}
	js[n - 1] = n - 1;
	for (k = n - 1; k >= 0; k--)
		if (js[k] != k) {
			t = b[k];
			b[k] = b[js[k]];
			b[js[k]] = t;
		}
	free(js);
	return (1);
}

int ISL_dnetn::run(int n, double eps, double t, double h, double *x, int k)
{
	int i, j, l;
	double am, z, beta, d, *y, *a, *b;
	y = alloc1double(n);
	a = alloc1double(n*n);
	b = alloc1double(n);
	l = k;
	am = 1.0 + eps;

	while (am >= eps) {
		dnetnf(x, b, n);
		am = 0.0;
		for (i = 0; i <= n - 1; i++) {
			z = fabs(b[i]);
			if (z > am)
				am = z;
		}
		if (am >= eps) {
			l = l - 1;
			if (l == 0) {
				free(y);
				free(b);
				free(a);
//				printf("fail\n");
				return (0);
			}
			for (j = 0; j <= n - 1; j++) {
				z = x[j];
				x[j] = x[j] + h;
				dnetnf(x, y, n);
				for (i = 0; i <= n - 1; i++)
					a[i * n + j] = y[i];
				x[j] = z;
			}
			if (agaus(a, b, n) == 0) {
				free(y);
				free(a);
				free(b);
				return (-1);
			}
			beta = 1.0;
			for (i = 0; i <= n - 1; i++)
				beta = beta - b[i];
			if (fabs(beta) + 1.0 == 1.0) {
				free(y);
				free(a);
				free(b);
//				printf("fail\n");
				return (-2);
			}
			d = h / beta;
			for (i = 0; i <= n - 1; i++)
				x[i] = x[i] - d * b[i];
			h = t * h;
		}
	}
	free(y);
	free(a);
	free(b);
	return (k - l);
}


}/*End of ISLib*/
