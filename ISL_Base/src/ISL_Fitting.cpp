/**
*	@file	ISL_Fitting.cpp
*	@brief	[Header file of Fitting Functions], 拟合；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang, Chen Ke], 刘百红、杨强、宋志翔、陈科；
*	@date	2014-06-04
*	@refer	SU CWP
*/

#include "ISL_Fitting.h"
#include "ISL_Numerical.h"
#include "ISL_Absorption.h"

namespace ISLib {

/*最小二乘曲线拟合,用最小二乘法求给定点的拟合多项式，此函数用于模拟地震记录的振幅谱*/
void ISL_lscf(float *x, float *y, int n, float *cof, int nterm, float *simu)
{
	int i, j, k;
	float z, p, c, g, q, d1, d2;
	float *s, *t, *b;

	s = new float[n];
	t = new float[n];
	b = new float[n];

	if (nterm > n)
		nterm = n;

	if (nterm > 20)
		nterm = 20;

	for (i = 0; i < nterm; i++)
		cof[i] = 0.0;

	z = 0.0;

	for (i = 0; i < n; i++)
		z += x[i] / float(1.0 * n);

	b[0] = 1.0;
	d1 = float(1.0 * n);
	p = 0.0;
	c = 0.0;

	for (i = 0; i < n; i++) {
		p = p + x[i] - z;
		c += y[i];
	}

	c /= d1;
	p /= d1;
	cof[0] = c * b[0];

	if (nterm > 1) {
		t[1] = 1.0;
		t[0] = -p;
		d2 = 0.0;
		c = 0.0;
		g = 0.0;

		for (i = 0; i < n; i++) {
			q = x[i] - z - p;
			d2 = d2 + q * q;
			c += y[i] * q;
			g = g + (x[i] - z) * q * q;
		}

		c /= d2;
		p = g / d2;
		q = d2 / d1;
		d1 = d2;
		cof[1] = c * t[1];
		cof[0] = c * t[0] + cof[0];
	}

	for (j = 2; j < nterm; j++) {
		s[j] = t[j - 1];
		s[j - 1] = -p * t[j - 1] + t[j - 2];

		if (j >= 3)
			for (k = j - 2; k >= 1; k--)
				s[k] = -p * t[k] + t[k - 1] - q * b[k];

		s[0] = -p * t[0] - q * b[0];
		d2 = 0.0;
		c = 0.0;
		g = 0.0;

		for (i = 0; i < n; i++) {
			q = s[j];

			for (k = j - 1; k >= 0; k--)
				q = q * (x[i] - z) + s[k];

			d2 = d2 + q * q;
			c += y[i] * q;
			g = g + (x[i] - z) * q * q;
		}

		c /= d2;
		p = g / d2;
		q = d2 / d1;
		d1 = d2;
		cof[j] = c * s[j];
		t[j] = s[j];

		for (k = j - 1; k >= 0; k--) {
			cof[k] = c * s[k] + cof[k];
			b[k] = t[k];
			t[k] = s[k];
		}
	}

	memset(simu, 0, sizeof(float) * n);

	for (i = 0; i < n; i++) {
		q = cof[nterm - 1];

		for (k = nterm - 2; k >= 0; k--)
			q = cof[k] + q * (x[i] - z);

		simu[i] = q;
		p = q - y[i];

        if (fabs(p)>simu[2])
        	simu[2] = fabs(p);
        simu[0] = p * p;
        simu[1] = fabs(p);
	}

	delete[] s;
	delete[] t;
	delete[] b;
}


/*拟合功率谱函数*/
void ISL_spectrumCurveFitting(float *inFreq, float *inSpec, int num, float deltaFre, float *outSpec)
{
	int nterm = 7;
	float *SpecX = new float[num];
	float *SpecY = new float[num];
	float *simu = new float[num];
	float *Coef = new float[nterm];
	for (int i = 0; i < num; i++) {
		SpecX[i] = inFreq[i];
		if (SpecX[i] == 0)
			SpecX[i] = deltaFre / 2;
		SpecY[i] = inSpec[i];
		if (SpecY[i] == 0)
			SpecY[i] = SpecY[i - 1];
		else
			SpecY[i] = log(SpecY[i] / pow(SpecX[i], 2));
	}

	ISL_lscf(SpecX, SpecY, num, Coef, nterm, simu);
	memset(outSpec, 0, sizeof(float) * num);
	for (int i = 0; i < num; i++)
		outSpec[i] = SpecX[i] * SpecX[i] * exp(simu[i]);

	delete[] SpecX;
	delete[] SpecY;
	delete[] simu;
	delete[] Coef;
}


/* 随机样本分析 */
void ISL_rhis(double *x, int n, double x0, double h, int m, int l, double *dt, int *g, int *q)
{
	int i, j;
	double s;
	dt[0] = 0.0;

	for (i = 0; i <= n - 1; i++)
		dt[0] = dt[0] + x[i] / n;

	dt[1] = 0.0;
	for (i = 0; i <= n - 1; i++)
		dt[1] = dt[1] + (x[i] - dt[0]) * (x[i] - dt[0]);
	dt[1] = dt[1] / n;
	dt[2] = sqrt(dt[1]);

	for (i = 0; i <= m - 1; i++) {
		q[i] = 0;
		s = x0 + (i + 0.5) * h - dt[0];
		s = exp(-s * s / (2.0 * dt[1]));
		g[i] = n * s * h / (dt[2] * 2.5066);
	}
	s = x0 + m * h;

	for (i = 0; i <= n - 1; i++)
		if ((x[i] - x0) >= 0.0)
			if ((s - x[i]) >= 0.0) {
				j = (x[i] - x0) / h;
				q[j] = q[j] + 1;
			}
	if (l == 0)
		return;
	return;
}

/*一元线性回归分析*/
void ISL_sqt1(double *x, double *y, int n, double *a, double *dt)
{
	int i;
	double xx, yy, e, f, q, u, p, umax, umin, s;
	xx = 0.0;
	yy = 0.0;
	for (i = 0; i <= n - 1; i++) {
		xx = xx + x[i] / n;
		yy = yy + y[i] / n;
	}

	e = 0.0;
	f = 0.0;
	for (i = 0; i <= n - 1; i++) {
		q = x[i] - xx;
		e = e + q * q;
		f = f + q * (y[i] - yy);
	}

	a[1] = f / e;
	a[0] = yy - a[1] * xx;
	q = 0.0;
	u = 0.0;
	p = 0.0;
	umax = 0.0;
	umin = 1.0e+30;

	for (i = 0; i <= n - 1; i++) {
		s = a[1] * x[i] + a[0];
		q = q + (y[i] - s) * (y[i] - s);
		p = p + (s - yy) * (s - yy);
		e = fabs(y[i] - s);
		if (e > umax)
			umax = e;
		if (e < umin)
			umin = e;
		u = u + e / n;
	}

	dt[1] = sqrt(q / n);
	dt[0] = q;
	dt[2] = p;
	dt[3] = umax;
	dt[4] = umin;
	dt[5] = u;

	return;
}


/*多元线性回归分析*/
void ISL_sqt2(double *x, double *y, int m, int n, double *a, double *dt, double *v)
{
	int i, j, k, mm;
	double q, e, u, p, yy, s, r, pp, *b;

	b = new double[(m + 1) * (m + 1)];

	mm = m + 1;
	b[mm * mm - 1] = n;
	for (j = 0; j <= m - 1; j++) {
		p = 0.0;
		for (i = 0; i <= n - 1; i++)
			p = p + x[j * n + i];
		b[m * mm + j] = p;
		b[j * mm + m] = p;
	}

	for (i = 0; i <= m - 1; i++)
		for (j = i; j <= m - 1; j++) {
			p = 0.0;
			for (k = 0; k <= n - 1; k++)
				p = p + x[i * n + k] * x[j * n + k];
			b[j * mm + i] = p;
			b[i * mm + j] = p;
		}

	a[m] = 0.0;
	for (i = 0; i <= n - 1; i++)
		a[m] = a[m] + y[i];

	for (i = 0; i <= m - 1; i++) {
		a[i] = 0.0;
		for (j = 0; j <= n - 1; j++)
			a[i] = a[i] + x[i * n + j] * y[j];
	}

	ISL_chlk(b, mm, 1, a);

	yy = 0.0;
	for (i = 0; i <= n - 1; i++)
		yy = yy + y[i] / n;

	q = 0.0;
	e = 0.0;
	u = 0.0;

	for (i = 0; i <= n - 1; i++) {
		p = a[m];
		for (j = 0; j <= m - 1; j++)
			p = p + a[j] * x[j * n + i];
		q = q + (y[i] - p) * (y[i] - p);
		e = e + (y[i] - yy) * (y[i] - yy);
		u = u + (yy - p) * (yy - p);
	}

	s = sqrt(q / n);
	r = sqrt(1.0 - q / e);
	for (j = 0; j <= m - 1; j++) {
		p = 0.0;
		for (i = 0; i <= n - 1; i++) {
			pp = a[m];
			for (k = 0; k <= m - 1; k++)
				if (k != j)
					pp = pp + a[k] * x[k * n + i];
			p = p + (y[i] - pp) * (y[i] - pp);
		}
		v[j] = sqrt(1.0 - q / p);
	}

	dt[0] = q;
	dt[1] = s;
	dt[2] = r;
	dt[3] = u;

	if(b){ delete []b; b = NULL; }

	return;
}


/*逐步回归分析*/
void ISL_sqt3(int n, int k, double *x, double f1, double f2, double eps,
				double *xx, double *b, double *v, double *s, double *dt,
				double *ye, double *yr, double *r)
{
	int i, j, ii, m, imi, imx, l, it;
	double z, phi, sd, vmi, vmx, q, fmi, fmx;
	m = n + 1;
	q = 0.0;

	for (j = 0; j <= n; j++) {
		z = 0.0;
		for (i = 0; i <= k - 1; i++)
			z = z + x[i * m + j] / k;
		xx[j] = z;
	}

	for (i = 0; i <= n; i++)
		for (j = 0; j <= i; j++) {
			z = 0.0;
			for (ii = 0; ii <= k - 1; ii++)
				z = z + (x[ii * m + i] - xx[i]) * (x[ii * m + j] - xx[j]);
			r[i * m + j] = z;
		}

	for (i = 0; i <= n; i++)
		ye[i] = sqrt(r[i * m + i]);
	for (i = 0; i <= n; i++)
		for (j = 0; j <= i; j++) {
			r[i * m + j] = r[i * m + j] / (ye[i] * ye[j]);
			r[j * m + i] = r[i * m + j];
		}

	phi = k - 1.0;
	sd = ye[n] / sqrt(k - 1.0);
	it = 1;
	while (it == 1) {
		it = 0;
		vmi = 1.0e+35;
		vmx = 0.0;
		imi = -1;
		imx = -1;
		for (i = 0; i <= n; i++) {
			v[i] = 0.0;
			b[i] = 0.0;
			s[i] = 0.0;
		}

		for (i = 0; i <= n - 1; i++)
			if (r[i * m + i] >= eps) {
				v[i] = r[i * m + n] * r[n * m + i] / r[i * m + i];
				if (v[i] >= 0.0) {
					if (v[i] > vmx) {
						vmx = v[i];
						imx = i;
					}
				} else {
					b[i] = r[i * m + n] * ye[n] / ye[i];
					s[i] = sqrt(r[i * m + i]) * sd / ye[i];
					if (fabs(v[i]) < vmi) {
						vmi = fabs(v[i]);
						imi = i;
					}
				}
			}

		if (phi != n - 1.0) {
			z = 0.0;
			for (i = 0; i <= n - 1; i++)
				z = z + b[i] * xx[i];
			b[n] = xx[n] - z;
			s[n] = sd;
			v[n] = q;
		} else {
			b[n] = xx[n];
			s[n] = sd;
		}

		fmi = vmi * phi / r[n * m + n];
		fmx = (phi - 1.0) * vmx / (r[n * m + n] - vmx);
		if ((fmi < f2) || (fmx >= f1)) {
			if (fmi < f2) {
				phi = phi + 1.0;
				l = imi;
			} else {
				phi = phi - 1.0;
				l = imx;
			}
			for (i = 0; i <= n; i++)
				if (i != l)
					for (j = 0; j <= n; j++)
						if (j != l)
							r[i * m + j] = r[i * m + j] - (r[l * m + j] / r[l
									* m + l]) * r[i * m + l];
			for (j = 0; j <= n; j++)
				if (j != l)
					r[l * m + j] = r[l * m + j] / r[l * m + l];
			for (i = 0; i <= n; i++)
				if (i != l)
					r[i * m + l] = -r[i * m + l] / r[l * m + l];
			r[l * m + l] = 1.0 / r[l * m + l];
			q = r[n * m + n] * ye[n] * ye[n];
			sd = sqrt(r[n * m + n] / phi) * ye[n];
			dt[0] = sqrt(1.0 - r[n * m + n]);
			dt[1] = (phi * (1.0 - r[n * m + n])) / ((k - phi - 1.0) * r[n * m
					+ n]);
			it = 1;
		}
	}

	for (i = 0; i <= k - 1; i++) {
		z = 0.0;
		for (j = 0; j <= n - 1; j++)
			z = z + b[j] * x[i * m + j];
		ye[i] = b[n] + z;
		yr[i] = x[i * m + n] - ye[i];
	}

	return;
}

} /* End Of namespace ISLIB */
