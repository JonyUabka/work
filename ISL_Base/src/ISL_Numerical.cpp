/**
*	@file	ISL_NUMERICAL.cpp 
*	@brief	[Source file of NUMERICAL], NUMERICAL；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang, He Kai], 刘百红、杨强、宋志翔、何恺；
*	@date	2014-05-22
*	@refer	SU CWP
*/

#include "ISL_Numerical.h"
#include "ISL_Matrix.h"


namespace ISLib
{

/**********************************************************************************************************************
 *
 *  功能：一维多项式求值
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a					In			存放n-1次多项式的n个系数
 *		int					n					In			多项式的项数
 *		double				x					In			指定的自变量值
 *
 *  返回：double， 返回多项式值
 *
**********************************************************************************************************************/
/*  一维多项式求值  */
double ISL_plyv(double * a, int n, double x)
{
	double u;
	u = a[n - 1];
	for (int i = n - 2; i >= 0; i--)
		u = u * x + a[i];

	return u;
}


/**********************************************************************************************************************
 *
 *  功能：一维多项式多组求值
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[n]				In			存放n-1次多项式的n个系数,返回时其值将改变
 *		int					n					In			多项式的项数
 *		double*				x[m]				In			存放给定的m个自变量值
 *		int					m					In			给定自变量的个数
 *		double * 			p[m]				Out			数组大小为m， 返回时存放与给定m个自变量值对应的多项式值, 进入函数前请分配好内存空间
 *
 *  返回：double * ，数组大小为m， 返回时存放与给定m个自变量值对应的多项式值
 *
**********************************************************************************************************************/
/*  一维多项式多组求值  */
void ISL_plys(double *a, int n, double *x, int m, double *p)
{
	int i, j, mm, nn, ll, t, s, kk, k;
	double * b = new double[2 * n], y, z;

	y = a[n - 1];

	for (i = 0; i <= n - 1; i++)
		b[i] = a[i] / y;

	k = log(n - 0.5) / log(2.0) + 1;
	nn = 1;
	for (i = 0; i <= k - 1; i++)
		nn = 2 * nn;
	for (i = n; i < nn - 1; i++)
		b[i] = 0.0;

	b[nn - 1] = 1.0;
	t = nn;
	s = 1;

	for (i = 1; i <= k - 1; i++) {
		t = t / 2;
		mm = -t;

		for (j = 1; j <= s; j++) {
			mm = mm + t + t;
			b[mm - 1] = b[mm - 1] - 1.0;

			for (kk = 2; kk <= t; kk++)
				b[mm - kk] = b[mm - kk] - b[mm - 1] * b[mm + t - kk];
		}
		s = s + s;
	}

	for (kk = 1; kk <= m; kk++) {
		for (i = 0; i <= (nn - 2) / 2; i++)
			a[i] = x[kk - 1] + b[2 * i];
		mm = 1;
		z = x[kk - 1];

		for (i = 1; i <= k - 1; i++) {
			mm = mm + mm;
			ll = mm + mm;
			z = z * z;

			for (j = 0; j <= nn - 1; j = j + ll)
				a[j / 2] = a[j / 2] + a[(j + mm) / 2] * (z + b[j + mm - 1]);
		}
		z = z * z / x[kk - 1];
		if (nn != n)
			a[0] = a[0] - z;
		p[kk - 1] = a[0] * y;
	}

	if(b){ delete []b; b = NULL; }

	return;
}


/**********************************************************************************************************************
 *
 *  功能：二维多项式求值
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[m][n]				In			存放二维次多项式的系数
 *		int					m					In			自变量x的最高次数为m-1
 *		int					n					In			自变量y的最高次数为n-1
 *		double				x					In			给定的x的自变量值
 *		double				y					In			给定的y的自变量值
 *
 *  返回：double， 返回多项式值
 *
**********************************************************************************************************************/
/*  二维多项式求值  */
double ISL_bply(double *a, int m, int n, double x, double y)
{
	int i, j;
	double u, s, xx;

	u = 0.0;
	xx = 1.0;

	for (i = 0; i <= m - 1; i++) {
		s = a[i * n + n - 1] * xx;

		for (j = n - 2; j >= 0; j--)
			s = s * y + a[i * n + j] * xx;

		u = u + s;
		xx = xx * x;
	}
	return u;
}


/**********************************************************************************************************************
 *
 *  功能：复系数多项式求值
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				ar[n]				In			分别存放多项式系数的实部与虚部
 *							ai[n]				In
 *		int					n					In			多项式的项数，其最高次数为n-1
 *		double				x, y				In			给定复数Z的实部与虚部
 *		double &				&u, &v				In/Out		指向返回多项式值p（z）的实部与虚部
 *
 *  返回：无
 *
**********************************************************************************************************************/
/*  复系数多项式求值  */
void ISL_cply(double *ar, double *ai, int n, double x, double y,double &u, double &v)
{
	double p = 0, q = 0, s, t;

	s = ar[n - 1];
	t = ai[n - 1];
	for (int i = n - 2; i >= 0; i--) {
		ISL_cmul(s, t, x, y, p, q);
		s = p + ar[i];
		t = q + ai[i];
	}
	u = s;
	v = t;
	return;
}

void ISL_cmul(double a, double b, double c, double d,double &e, double &f)
{
	double p, q, s;

	p = a * c;
	q = b * d;
	s = (a + b) * (c + d);
	e = p - q;
	f = s - p - q;

	return;
}


/**********************************************************************************************************************
 *
 *  功能：多项式相乘
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				p[m]				In			存放多项式p(x)的系数
 *		int					m					In			多项式p(x)的项数，其最高次数为m-1
 *		double*				q[n]				In			存放多项式q(x)的系数
 *		int					n					In			多项式q(x)的项数，其最高次数为n-1
 *		double*				s[k]				In/Out		返回乘积多项式的系数, 进入函数前请分配好内存空间
 *		int					k					In			乘积多项式s(x)的项数，其最高次数为k-1，其中k=m+n-1
 *
 *  返回：无
 *
**********************************************************************************************************************/
/*  多项式相乘  */
void ISL_pmul(double *p, int m, double *q, int n, double *s, int k)
{
	int i, j;
	for (i = 0; i <= k - 1; i++)
		s[i] = 0.0;

	for (i = 0; i <= m - 1; i++)
		for (j = 0; j <= n - 1; j++)
			s[i + j] = s[i + j] + p[i] * q[j];

	return;
}


/**********************************************************************************************************************
 *
 *  功能：复系数多项式相乘
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				pr[m],pi[m]			In			存放多项式p(x)的系数的实部与虚部
 *		int					m					In			多项式p(x)的项数，其最高次数为m-1
 *		double*				qr[n],qi[n]			In			存放多项式q(x)的系数的实部与虚部
 *		int					n					In			多项式q(x)的项数，其最高次数为n-1
 *		double*				sr[k],si[k]			In/Out		返回乘积多项式的系数的实部与虚部, 进入函数前请分配好内存空间
 *		int					k					In			乘积多项式s(x)的项数，其最高次数为k-1，其中k=m+n-1
 *
 *  返回：无
 *
**********************************************************************************************************************/
/*  复系数多项式相乘  */
void ISL_cpml(double *pr, double * pi, int m, double *qr, double *qi, int n, double *sr, double *si, int k)
{
	int i, j;
	double a, b, c, d, u, v;

	for (i = 0; i <= k - 1; i++) {
		sr[i] = 0.0;
		si[i] = 0.0;
	}

	for (i = 0; i <= m - 1; i++)
		for (j = 0; j <= n - 1; j++) {
			a = pr[i];
			b = pi[i];
			c = qr[j];
			d = qi[j];
			ISL_cmul(a, b, c, d, u, v);
			sr[i + j] = sr[i + j] + u;
			si[i + j] = si[i + j] + v;
		}

	return;
}


/**********************************************************************************************************************
 *
 *  功能：多项式相除
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				p[m]				In			存放多项式p(x)的系数。返回时其中的值将被破坏
 *		int					m					In			多项式p(x)的项数，其最高次数为m-1
 *		double*				q[n]				In			存放多项式q(x)的系数
 *		int					n					In			多项式q(x)的项数，其最高次数为n-1
 *		double*				s[k]				In/Out		返回商多项式s（x）的系数, 进入函数前请分配好内存空间
 *		int					k					In			商多项式s(x)的项数，其最高次数为k-1，其中k=m+n-1
  *		double*				r[l]				In/Out		返回商多项式r（x）的系数, 进入函数前请分配好内存空间
 *		int					l					In			余多项式r(x)的项数，其最高次数为l-1，其中l=n-1
 *
 *  返回：无
 *
**********************************************************************************************************************/
/*  多项式相除  */
void ISL_pdiv(double *p, int m, double *q, int n, double *s, int k, double *r, int l)
{
	int i, j, mm, ll;

	for (i = 0; i <= k - 1; i++)
		s[i] = 0.0;

	if (q[n - 1] + 1.0 == 1.0)
		return;

	ll = m - 1;
	for (i = k; i >= 1; i--) {
		s[i - 1] = p[ll] / q[n - 1];
		mm = ll;

		for (j = 1; j <= n - 1; j++) {
			p[mm - 1] = p[mm - 1] - s[i - 1] * q[n - j - 1];
			mm = mm - 1;
		}
		ll = ll - 1;
	}

	for (i = 0; i <= l - 1; i++)
		r[i] = p[i];

	return;
}


/**********************************************************************************************************************
 *
 *  功能：复系数多项式相除
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				pr[m],pi[m]			In			存放多项式p(x)的系数的实部与虚部。返回时其中的值将被破坏
 *		int					m					In			多项式p(x)的项数，其最高次数为m-1
 *		double*				qr[n],pi[n]			In			存放多项式q(x)的系数的实部与虚部
 *		int					n					In			多项式q(x)的项数，其最高次数为n-1
 *		double*				sr[k],si[k]			In/Out		返回商多项式s（x）的系数的实部与虚部, 进入函数前请分配好内存空间
 *		int					k					In			商多项式s(x)的项数，其最高次数为k-1，其中k=m+n-1
  *		double*				rr[l],ri[l]			In/Out		返回商多项式r（x）的系数的实部与虚部, 进入函数前请分配好内存空间
 *		int					l					In			余多项式r(x)的项数，其最高次数为l-1，其中l=n-1
 *
 *  返回：无
 *
**********************************************************************************************************************/
/*  复系数多项式相除  */
void ISL_cpdv(double *pr, double *pi, int m,
			double *qr, double *qi, int n,
			double *sr, double *si, int k,
			double *rr, double *ri, int l)
{
	int i, j, mm, ll;
	double a, b, c, d, u, v;

	for (i = 0; i <= k - 1; i++) {
		sr[i] = 0.0;
		si[i] = 0.0;
	}

	d = qr[n - 1] * qr[n - 1] + qi[n - 1] * qi[n - 1];
	if (d + 1.0 == 1.0)
		return;

	ll = m - 1;
	for (i = k; i >= 1; i--) {
		a = pr[ll];
		b = pi[ll];
		c = qr[n - 1];
		d = qi[n - 1];
		ISL_cdiv(a, b, c, d, u, v);
		sr[i - 1] = u;
		si[i - 1] = v;
		mm = ll;

		for (j = 1; j <= n - 1; j++) {
			a = sr[i - 1];
			b = si[i - 1];
			c = qr[n - j - 1];
			d = qi[n - j - 1];
			ISL_cmul(a, b, c, d, u, v);
			pr[mm - 1] = pr[mm - 1] - u;
			pi[mm - 1] = pi[mm - 1] - v;
			mm = mm - 1;
		}
		ll = ll - 1;
	}

	for (i = 0; i <= l - 1; i++) {
		rr[i] = pr[i];
		ri[i] = pi[i];
	}
	return;
}

void ISL_cdiv(double a, double b, double c, double d, double &e, double &f)
{
	double p, q, s, w;

	p = a * c;
	q = -b * d;
	s = (a + b) * (c - d);
	w = c * c + d * d;

	if (w + 1.0 == 1.0) {
		e = 1.0e+35 * a / fabs(a);
		f = 1.0e+35 * b / fabs(b);
	} else {
		e = (p - q) / w;
		f = (s - p - q) / w;
	}
	return;
}


/**********************************************************************************************************************
 *
 *  功能：求解实系数方程组的全选主元高斯消去法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[n][n]				In			存放方程组的系数矩阵，返回时将被破坏
 *		double*				b[n]				In/Out		存放方程组右端的常数向量，返回方程组的解向量
 *		int					n					In			方程组的阶数
 *
 *  返回：函数返回整型标志位，若返回标志值为0，则表示程序工作失败（因系数矩阵奇异）；
 *  				若返回标志值不为0；则表示正常返回。
 *
**********************************************************************************************************************/
/*  求解实系数方程组的全选主元高斯消去法  */
int ISL_gaus(double *a, double *b, int n)
{
	int *js=NULL, l, k, i, j, is=0, p, q;
	double d, t;

	js = new int[n];
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
			if(js){ delete []js; js = NULL; }
			printf("fail\n");
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
		if(js){ delete []js; js = NULL; }
		printf("fail\n");
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

	if(js){ delete []js; js = NULL; }
	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：求解实系数方程组的全选主元高斯-约当消去法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[n][n]				In			存放方程组的系数矩阵，返回时将被破坏
 *		double*				b[n][m]				In/Out		存放方程组右端的m组常数向量，返回方程组的m组解向量
 *		int					n					In			方程组的阶数
 *		int					m					In			方程组右端常数向量的组数
 *
 *  返回：函数返回整型标志位，若返回标志值为0，则表示程序工作失败（因系数矩阵奇异）；
 *  				若返回标志值不为0；则表示正常返回。
 *
**********************************************************************************************************************/
/*  求解实系数方程组的全选主元高斯-约当消去法  */
int ISL_gjdn(double *a, double *b, int n, int m)
{
	int *js=NULL, l, k, i, j, is=0, p, q;
	double d, t;

	js = new int[n];
	l = 1;
	for (k = 0; k <= n - 1; k++) {
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
				for (j = 0; j <= m - 1; j++) {
					p = k * m + j;
					q = is * m + j;
					t = b[p];
					b[p] = b[q];
					b[q] = t;
				}
			}
		}

		if (l == 0) {
			if(js){ delete []js; js = NULL; }
			printf("fail\n");
			return (0);
		}
		d = a[k * n + k];
		for (j = k + 1; j <= n - 1; j++) {
			p = k * n + j;
			a[p] = a[p] / d;
		}
		for (j = 0; j <= m - 1; j++) {
			p = k * m + j;
			b[p] = b[p] / d;
		}
		for (j = k + 1; j <= n - 1; j++)
			for (i = 0; i <= n - 1; i++) {
				p = i * n + j;
				if (i != k)
					a[p] = a[p] - a[i * n + k] * a[k * n + j];
			}
		for (j = 0; j <= m - 1; j++)
			for (i = 0; i <= n - 1; i++) {
				p = i * m + j;
				if (i != k)
					b[p] = b[p] - a[i * n + k] * b[k * m + j];
			}
	}

	for (k = n - 1; k >= 0; k--)
		if (js[k] != k)
			for (j = 0; j <= m - 1; j++) {
				p = k * m + j;
				q = js[k] * m + j;
				t = b[p];
				b[p] = b[q];
				b[q] = t;
			}

	if(js){ delete []js; js = NULL; }
	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：求解复系数方程组的全选主元高斯消去法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				ar[n][n]，			In			存放方程组的复系数矩阵的实部和虚部，返回时将被破坏
 *							ai[n][n]
 *		double*				br[n], bi[n			In/Out		存放方程组右端的复常数向量的实部与虚部，返回方程组的解向量的实部与虚部
 *		int					n					In			方程组的阶数
 *
 *  返回：函数返回整型标志位，若返回标志值为0，则表示程序工作失败（因系数矩阵奇异）；
 *  				若返回标志值不为0；则表示正常返回。
 *
**********************************************************************************************************************/
/*  求解复系数方程组的全选主元高斯消去法  */
int ISL_cgas(double *ar, double *ai, int n, double *br, double *bi)
{
	int *js, l, k, i, j, is=0, u, v;
	double p, q, s, d;

	js =  new int[n];
	for (k = 0; k <= n - 2; k++) {
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++) {
				u = i * n + j;
				p = ar[u] * ar[u] + ai[u] * ai[u];
				if (p > d) {
					d = p;
					js[k] = j;
					is = i;
				}
			}
		if (d + 1.0 == 1.0) {
			if(js){ delete []js; js = NULL; }
			printf("err**fail\n");
			return (0);
		}

		if (is != k) {
			for (j = k; j <= n - 1; j++) {
				u = k * n + j;
				v = is * n + j;
				p = ar[u];
				ar[u] = ar[v];
				ar[v] = p;
				p = ai[u];
				ai[u] = ai[v];
				ai[v] = p;
			}
			p = br[k];
			br[k] = br[is];
			br[is] = p;
			p = bi[k];
			bi[k] = bi[is];
			bi[is] = p;
		}

		if (js[k] != k)
			for (i = 0; i <= n - 1; i++) {
				u = i * n + k;
				v = i * n + js[k];
				p = ar[u];
				ar[u] = ar[v];
				ar[v] = p;
				p = ai[u];
				ai[u] = ai[v];
				ai[v] = p;
			}

		v = k * n + k;
		for (j = k + 1; j <= n - 1; j++) {
			u = k * n + j;
			p = ar[u] * ar[v];
			q = -ai[u] * ai[v];
			s = (ar[v] - ai[v]) * (ar[u] + ai[u]);
			ar[u] = (p - q) / d;
			ai[u] = (s - p - q) / d;
		}

		p = br[k] * ar[v];
		q = -bi[k] * ai[v];
		s = (ar[v] - ai[v]) * (br[k] + bi[k]);
		br[k] = (p - q) / d;
		bi[k] = (s - p - q) / d;
		for (i = k + 1; i <= n - 1; i++) {
			u = i * n + k;
			for (j = k + 1; j <= n - 1; j++) {
				v = k * n + j;
				l = i * n + j;
				p = ar[u] * ar[v];
				q = ai[u] * ai[v];
				s = (ar[u] + ai[u]) * (ar[v] + ai[v]);
				ar[l] = ar[l] - p + q;
				ai[l] = ai[l] - s + p + q;
			}
			p = ar[u] * br[k];
			q = ai[u] * bi[k];
			s = (ar[u] + ai[u]) * (br[k] + bi[k]);
			br[i] = br[i] - p + q;
			bi[i] = bi[i] - s + p + q;
		}
	}

	u = (n - 1) * n + n - 1;
	d = ar[u] * ar[u] + ai[u] * ai[u];
	if (d + 1.0 == 1.0) {
		if(js){ delete []js; js = NULL; }
		printf("err**fail\n");
		return (0);
	}

	p = ar[u] * br[n - 1];
	q = -ai[u] * bi[n - 1];
	s = (ar[u] - ai[u]) * (br[n - 1] + bi[n - 1]);
	br[n - 1] = (p - q) / d;
	bi[n - 1] = (s - p - q) / d;
	for (i = n - 2; i >= 0; i--)
		for (j = i + 1; j <= n - 1; j++) {
			u = i * n + j;
			p = ar[u] * br[j];
			q = ai[u] * bi[j];
			s = (ar[u] + ai[u]) * (br[j] + bi[j]);
			br[i] = br[i] - p + q;
			bi[i] = bi[i] - s + p + q;
		}

	js[n - 1] = n - 1;
	for (k = n - 1; k >= 0; k--)
		if (js[k] != k) {
			p = br[k];
			br[k] = br[js[k]];
			br[js[k]] = p;
			p = bi[k];
			bi[k] = bi[js[k]];
			bi[js[k]] = p;
		}

	if(js){ delete []js; js = NULL; }
	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：求解复系数方程组的全选主元高斯-约当消去法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				ar[n][n]			In			存放方程组的复系数矩阵的实部和虚部，返回时将被破坏
 *							ai[n][n]
 *		double*				br[n][m]			In/Out		存放方程组右端的m组常数向量的实部和虚部，返回方程组的m组解向量的实部和虚部
 *							bi[n][m]
 *		int					n					In			方程组的阶数
 *		int					m					In			方程组右端常数向量的组数
 *
 *  返回：函数返回整型标志位，若返回标志值为0，则表示程序工作失败（因系数矩阵奇异）；
 *  				若返回标志值不为0；则表示正常返回。
 *
**********************************************************************************************************************/
/*  求解复系数方程组的全选主元高斯-约当消去法  */
int ISL_cjdn(double *ar, double *ai, double *br, double *bi, int n, int m)
{
	int *js=NULL, l, k, i, j, is=0, u, v;
	double p, q, s, d;
	js = new int[n];

	for (k = 0; k <= n - 1; k++) {
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++) {
				u = i * n + j;
				p = ar[u] * ar[u] + ai[u] * ai[u];
				if (p > d) {
					d = p;
					js[k] = j;
					is = i;
				}
			}

		if (d + 1.0 == 1.0) {
			if(js){ delete []js; js = NULL; }
			printf("err**fail\n");
			return (0);
		}

		if (is != k) {
			for (j = k; j <= n - 1; j++) {
				u = k * n + j;
				v = is * n + j;
				p = ar[u];
				ar[u] = ar[v];
				ar[v] = p;
				p = ai[u];
				ai[u] = ai[v];
				ai[v] = p;
			}
			for (j = 0; j <= m - 1; j++) {
				u = k * m + j;
				v = is * m + j;
				p = br[u];
				br[u] = br[v];
				br[v] = p;
				p = bi[u];
				bi[u] = bi[v];
				bi[v] = p;
			}
		}

		if (js[k] != k)
			for (i = 0; i <= n - 1; i++) {
				u = i * n + k;
				v = i * n + js[k];
				p = ar[u];
				ar[u] = ar[v];
				ar[v] = p;
				p = ai[u];
				ai[u] = ai[v];
				ai[v] = p;
			}

		v = k * n + k;
		for (j = k + 1; j <= n - 1; j++) {
			u = k * n + j;
			p = ar[u] * ar[v];
			q = -ai[u] * ai[v];
			s = (ar[v] - ai[v]) * (ar[u] + ai[u]);
			ar[u] = (p - q) / d;
			ai[u] = (s - p - q) / d;
		}

		for (j = 0; j <= m - 1; j++) {
			u = k * m + j;
			p = br[u] * ar[v];
			q = -bi[u] * ai[v];
			s = (ar[v] - ai[v]) * (br[u] + bi[u]);
			br[u] = (p - q) / d;
			bi[u] = (s - p - q) / d;
		}

		for (i = 0; i <= n - 1; i++)
			if (i != k) {
				u = i * n + k;
				for (j = k + 1; j <= n - 1; j++) {
					v = k * n + j;
					l = i * n + j;
					p = ar[u] * ar[v];
					q = ai[u] * ai[v];
					s = (ar[u] + ai[u]) * (ar[v] + ai[v]);
					ar[l] = ar[l] - p + q;
					ai[l] = ai[l] - s + p + q;
				}
				for (j = 0; j <= m - 1; j++) {
					l = i * m + j;
					v = k * m + j;
					p = ar[u] * br[v];
					q = ai[u] * bi[v];
					s = (ar[u] + ai[u]) * (br[v] + bi[v]);
					br[l] = br[l] - p + q;
					bi[l] = bi[l] - s + p + q;
				}
			}
	}

	for (k = n - 1; k >= 0; k--)
		if (js[k] != k)
			for (j = 0; j <= m - 1; j++) {
				u = k * m + j;
				v = js[k] * m + j;
				p = br[u];
				br[u] = br[v];
				br[v] = p;
				p = bi[u];
				bi[u] = bi[v];
				bi[v] = p;
			}

	if(js){ delete []js; js = NULL; }
	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：求解三对角线方程组的追赶法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				b[m]				In			以行为主，存放三对角矩阵中三条对角线上的元素
 *		int					n					In			方程组的阶数
 *		int					m					In			三对角矩阵三条对角线上的元素个数，其值应为m=3n-2
 *		double*				d[n]				In/Out		存放方程组右端的常数向量，返回方程组的解向量
 *
 *  返回：函数返回整型标志位，若返回标志值小于<0，则表示m的值不正确；
 *  				若返回标志值=0，则表示程序工作失败；
 *  				若返回标志值>0；则表示正常返回。
 *
**********************************************************************************************************************/
/*  求解三对角线方程组的追赶法  */
int ISL_trde(double *b, int n, int m, double *d)
{
	int k, j;
	double s;

	if (m != (3 * n - 2)) {
		printf("err\n");
		return (-2);
	}

	for (k = 0; k <= n - 2; k++) {
		j = 3 * k;
		s = b[j];
		if (fabs(s) + 1.0 == 1.0) {
			printf("fail\n");
			return (0);
		}
		b[j + 1] = b[j + 1] / s;
		d[k] = d[k] / s;
		b[j + 3] = b[j + 3] - b[j + 2] * b[j + 1];
		d[k + 1] = d[k + 1] - b[j + 2] * d[k];
	}

	s = b[3 * n - 3];
	if (fabs(s) + 1.0 == 1.0) {
		printf("fail\n");
		return (0);
	}

	d[n - 1] = d[n - 1] / s;
	for (k = n - 2; k >= 0; k--)
		d[k] = d[k] - b[3 * k + 1] * d[k + 1];

	return (2);
}


/**********************************************************************************************************************
 *
 *  功能：求解一般带型方程组
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				b[n][il]			In/Out		存放带型矩阵A中带区内的元素。返回时将被破坏
 *		double*				d[n][m]				In/Out		存放方程组右端的m组的常数向量。返回方程组的m组解向量
 *		int					n					In			方程组的阶数
 *		int					l					In			系数矩阵的半带宽h
 *		int					il					In			系数矩阵的带宽2h+1。应满足il=2l+1
 *		int					m					In			方程组右端常数向量的组数
 *
 *  返回：函数返回整型标志位，若返回标志值小于<0，则表示m的值不正确；
 *  				若返回标志值=0，则表示程序工作失败；
 *  				若返回标志值>0；则表示正常返回。
 *
**********************************************************************************************************************/
/*  求解一般带型方程组  */
int ISL_band(double *b, double *d, int n, int l, int il, int m)
{
	int ls, k, i, j, is=0, u, v;
	double p, t;
	if (il != (2 * l + 1)) {
		printf("fail\n");
		return (-2);
	}

	ls = l;
	for (k = 0; k <= n - 2; k++) {
		p = 0.0;
		for (i = k; i <= ls; i++) {
			t = fabs(b[i * il]);
			if (t > p) {
				p = t;
				is = i;
			}
		}
		if (p + 1.0 == 1.0) {
			printf("fail\n");
			return (0);
		}
		for (j = 0; j <= m - 1; j++) {
			u = k * m + j;
			v = is * m + j;
			t = d[u];
			d[u] = d[v];
			d[v] = t;
		}
		for (j = 0; j <= il - 1; j++) {
			u = k * il + j;
			v = is * il + j;
			t = b[u];
			b[u] = b[v];
			b[v] = t;
		}
		for (j = 0; j <= m - 1; j++) {
			u = k * m + j;
			d[u] = d[u] / b[k * il];
		}
		for (j = 1; j <= il - 1; j++) {
			u = k * il + j;
			b[u] = b[u] / b[k * il];
		}
		for (i = k + 1; i <= ls; i++) {
			t = b[i * il];
			for (j = 0; j <= m - 1; j++) {
				u = i * m + j;
				v = k * m + j;
				d[u] = d[u] - t * d[v];
			}
			for (j = 1; j <= il - 1; j++) {
				u = i * il + j;
				v = k * il + j;
				b[u - 1] = b[u] - t * b[v];
			}
			u = i * il + il - 1;
			b[u] = 0.0;
		}
		if (ls != (n - 1))
			ls = ls + 1;
	}

	p = b[(n - 1) * il];
	if (fabs(p) + 1.0 == 1.0) {
		printf("fail\n");
		return (0);
	}

	for (j = 0; j <= m - 1; j++) {
		u = (n - 1) * m + j;
		d[u] = d[u] / p;
	}

	ls = 1;
	for (i = n - 2; i >= 0; i--) {
		for (k = 0; k <= m - 1; k++) {
			u = i * m + k;
			for (j = 1; j <= ls; j++) {
				v = i * il + j;
				is = (i + j) * m + k;
				d[u] = d[u] - b[v] * d[is];
			}
		}
		if (ls != (il - 1))
			ls = ls + 1;
	}

	return (2);
}


/**********************************************************************************************************************
 *
 *  功能：求解对称方程组的分解法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[n][n]				In/Out		存放方程组的系数矩阵（应为对称矩阵)。返回时将被破坏
 *		int					n					In			方程组的阶数
 *		int					m					In			方程组右端常数向量的组数
 *		double*				c[n][m]				In/Out		存放方程组右端m组常数向量。返回方程组的m组解向量
 *
 *  返回：函数返回标志值。若返回的标志值小于0，则表示程序工作失败，若返回的标志值大于0，则表示正常返回
 *
**********************************************************************************************************************/
/*  求解对称方程组的分解法  */
int ISL_ldle(double *a, int n, int m, double *c)
{
	int i, j, l, k, u, v, w, k1, k2, k3;
	double p;
	if (fabs(a[0]) + 1.0 == 1.0) {
		printf("fail\n");
		return (-2);
	}

	for (i = 1; i <= n - 1; i++) {
		u = i * n;
		a[u] = a[u] / a[0];
	}

	for (i = 1; i <= n - 2; i++) {
		u = i * n + i;
		for (j = 1; j <= i; j++) {
			v = i * n + j - 1;
			l = (j - 1) * n + j - 1;
			a[u] = a[u] - a[v] * a[v] * a[l];
		}
		p = a[u];
		if (fabs(p) + 1.0 == 1.0) {
			printf("fail\n");
			return (-2);
		}
		for (k = i + 1; k <= n - 1; k++) {
			u = k * n + i;
			for (j = 1; j <= i; j++) {
				v = k * n + j - 1;
				l = i * n + j - 1;
				w = (j - 1) * n + j - 1;
				a[u] = a[u] - a[v] * a[l] * a[w];
			}
			a[u] = a[u] / p;
		}
	}

	u = n * n - 1;
	for (j = 1; j <= n - 1; j++) {
		v = (n - 1) * n + j - 1;
		w = (j - 1) * n + j - 1;
		a[u] = a[u] - a[v] * a[v] * a[w];
	}

	p = a[u];
	if (fabs(p) + 1.0 == 1.0) {
		printf("fail\n");
		return (-2);
	}

	for (j = 0; j <= m - 1; j++)
		for (i = 1; i <= n - 1; i++) {
			u = i * m + j;
			for (k = 1; k <= i; k++) {
				v = i * n + k - 1;
				w = (k - 1) * m + j;
				c[u] = c[u] - a[v] * c[w];
			}
		}

	for (i = 1; i <= n - 1; i++) {
		u = (i - 1) * n + i - 1;
		for (j = i; j <= n - 1; j++) {
			v = (i - 1) * n + j;
			w = j * n + i - 1;
			a[v] = a[u] * a[w];
		}
	}

	for (j = 0; j <= m - 1; j++) {
		u = (n - 1) * m + j;
		c[u] = c[u] / p;
		for (k = 1; k <= n - 1; k++) {
			k1 = n - k;
			k3 = k1 - 1;
			u = k3 * m + j;
			for (k2 = k1; k2 <= n - 1; k2++) {
				v = k3 * n + k2;
				w = k2 * m + j;
				c[u] = c[u] - a[v] * c[w];
			}
			c[u] = c[u] / a[k3 * n + k3];
		}
	}

	return (2);
}


/**********************************************************************************************************************
 *
 *  功能：求解对称正定方程组的平方根法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[n][n]				In/Out		存放对称正定的系数矩阵。返回时其上三角部分存放分解后的U矩
 *		int					n					In			方程组的阶数
 *		int					m					In			方程组右端常数向量的组数
 *		double*				d[n][m]				In/Out		存放方程组右端m组常数向量。返回方程组的m组解向量
 *
 *  返回：函数返回标志值。若返回的标志值小于0，则表示程序工作失败，若返回的标志值大于0，则表示正常返回
 *
**********************************************************************************************************************/
/*  求解对称正定方程组的平方根法  */
int ISL_chlk(double *a, int n, int m, double *d)
{
	int i, j, k, u, v;
	if ((a[0] + 1.0 == 1.0) || (a[0] < 0.0)) {
		printf("fail\n");
		return (-2);
	}
	a[0] = sqrt(a[0]);
	for (j = 1; j <= n - 1; j++)
		a[j] = a[j] / a[0];
	for (i = 1; i <= n - 1; i++) {
		u = i * n + i;
		for (j = 1; j <= i; j++) {
			v = (j - 1) * n + i;
			a[u] = a[u] - a[v] * a[v];
		}
		if ((a[u] + 1.0 == 1.0) || (a[u] < 0.0)) {
			printf("fail\n");
			return (-2);
		}
		a[u] = sqrt(a[u]);
		if (i != (n - 1)) {
			for (j = i + 1; j <= n - 1; j++) {
				v = i * n + j;
				for (k = 1; k <= i; k++)
					a[v] = a[v] - a[(k - 1) * n + i] * a[(k - 1) * n + j];
				a[v] = a[v] / a[u];
			}
		}
	}
	for (j = 0; j <= m - 1; j++) {
		d[j] = d[j] / a[0];
		for (i = 1; i <= n - 1; i++) {
			u = i * n + i;
			v = i * m + j;
			for (k = 1; k <= i; k++)
				d[v] = d[v] - a[(k - 1) * n + i] * d[(k - 1) * m + j];
			d[v] = d[v] / a[u];
		}
	}
	for (j = 0; j <= m - 1; j++) {
		u = (n - 1) * m + j;
		d[u] = d[u] / a[n * n - 1];
		for (k = n - 1; k >= 1; k--) {
			u = (k - 1) * m + j;
			for (i = k; i <= n - 1; i++) {
				v = (k - 1) * n + i;
				d[u] = d[u] - a[v] * d[i * m + j];
			}
			v = (k - 1) * n + k - 1;
			d[u] = d[u] / a[v];
		}
	}
	return (2);
}

/**********************************************************************************************************************
 *
 *  功能：求解大型稀疏方程组
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[n][n]				In/Out		存放对称正定的系数矩阵。返回时其上三角部分存放分解后的U矩
 *		int					n					In			方程组的阶数
 *		double*				b[n][m]				In/Out		存放方程组右端常数向量。返回方程组的解向量
 *
 *  返回：函数返回标志值。若返回的标志值为0，则表示系数矩阵奇异;若返回的标志值不为0，则表示正常返回
 *
**********************************************************************************************************************/
/*  求解大型稀疏方程组  */
int ISL_ggje(double *a, int n, double *b)
{
	int *js=NULL, i, j, k, is=0, u, v;
	double d, t;
	js = new int[n];

	for (k = 0; k <= n - 1; k++) {
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

		if (d + 1.0 == 1.0) {
			if(js){ delete []js; js = NULL; }
			printf("fail\n");
			return (0);
		}

		if (is != k) {
			for (j = k; j <= n - 1; j++) {
				u = k * n + j;
				v = is * n + j;
				t = a[u];
				a[u] = a[v];
				a[v] = t;
			}
			t = b[k];
			b[k] = b[is];
			b[is] = t;
		}

		if (js[k] != k)
			for (i = 0; i <= n - 1; i++) {
				u = i * n + k;
				v = i * n + js[k];
				t = a[u];
				a[u] = a[v];
				a[v] = t;
			}

		t = a[k * n + k];
		for (j = k + 1; j <= n - 1; j++) {
			u = k * n + j;
			if (a[u] != 0.0)
				a[u] = a[u] / t;
		}

		b[k] = b[k] / t;
		for (j = k + 1; j <= n - 1; j++) {
			u = k * n + j;
			if (a[u] != 0.0) {
				for (i = 0; i <= n - 1; i++) {
					v = i * n + k;
					if ((i != k) && (a[v] != 0.0)) {
						is = i * n + j;
						a[is] = a[is] - a[v] * a[u];
					}
				}
			}
		}

		for (i = 0; i <= n - 1; i++) {
			u = i * n + k;
			if ((i != k) && (a[u] != 0.0))
				b[i] = b[i] - a[u] * b[k];
		}
	}

	for (k = n - 1; k >= 0; k--)
		if (k != js[k]) {
			t = b[k];
			b[k] = b[js[k]];
			b[js[k]] = t;
		}

	if(js){ delete []js; js = NULL; }
	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：求解托伯利兹方程组的列文逊方法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				t[n]				In			存放n阶T型矩阵中的元素
 *		int					n					In			方程组的阶数
 *		double*				b[n]				In			存放方程组右端的常数向量
 *		double*				x[n]				Out			返回方程组的解向量
 *
 *  返回：函数返回标志值。若返回的标志值小于0，则表示程序工作失败，若返回的标志值大于0，则表示正常返回
 *
**********************************************************************************************************************/
/*  求解托伯利兹方程组的列文逊方法  */
int ISL_tlvs(double *t, int n, double *b, double *x)
{
	int i, j, k;
	double a, beta, q, c, h, *y, *s;
	s = new double[n];
	y = new double[n];
	a = t[0];

	if (fabs(a) + 1.0 == 1.0) {
		if(s){ delete []s; s = NULL; }
		if(y){ delete []y; y = NULL; }
		printf("fail\n");
		return (-1);
	}

	y[0] = 1.0;
	x[0] = b[0] / a;
	for (k = 1; k <= n - 1; k++) {
		beta = 0.0;
		q = 0.0;
		for (j = 0; j <= k - 1; j++) {
			beta = beta + y[j] * t[j + 1];
			q = q + x[j] * t[k - j];
		}

		if (fabs(a) + 1.0 == 1.0) {
			free(s);
			free(y);
			printf("fail\n");
			return (-1);
		}

		c = -beta / a;
		s[0] = c * y[k - 1];
		y[k] = y[k - 1];
		if (k != 1)
			for (i = 1; i <= k - 1; i++)
				s[i] = y[i - 1] + c * y[k - i - 1];
		a = a + c * beta;

		if (fabs(a) + 1.0 == 1.0) {
			if(s){ delete []s; s = NULL; }
			if(y){ delete []y; y = NULL; }
			printf("fail\n");
			return (-1);
		}

		h = (b[k] - q) / a;
		for (i = 0; i <= k - 1; i++) {
			x[i] = x[i] + h * s[i];
			y[i] = s[i];
		}
		x[k] = h * y[k];
	}

	if(s){ delete []s; s = NULL; }
	if(y){ delete []y; y = NULL; }
	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：高斯-赛德尔迭代法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[n][n]				In			存放方程组的系数矩阵
 *		double*				b[n]				In			存放方程组右端的常数向量
 *		int					n					In			方程组的阶数
 *		double				eps					In			给定的精度要求
 *		double*				x[n]				Out			返回方程组的解向量
 *
 *  返回：函数返回标志值。若返回的标志值小于0，则表示系数矩阵不具有主对角线占绝对优势，若返回的标志值大于0，则表示正常返回
 *
**********************************************************************************************************************/
/*  高斯-赛德尔迭代法  */
int ISL_gsdl(double *a, double *b, int n, double *x, double eps)
{
	int i, j, u, v;
	double p, t, s, q;

	for (i = 0; i <= n - 1; i++) {
		u = i * n + i;
		p = 0.0;
		x[i] = 0.0;
		for (j = 0; j <= n - 1; j++)
			if (i != j) {
				v = i * n + j;
				p = p + fabs(a[v]);
			}
		if (p >= fabs(a[u])) {
			printf("fail\n");
			return (-1);
		}
	}

	p = eps + 1.0;
	while (p >= eps) {
		p = 0.0;
		for (i = 0; i <= n - 1; i++) {
			t = x[i];
			s = 0.0;
			for (j = 0; j <= n - 1; j++)
				if (j != i)
					s = s + a[i * n + j] * x[j];
			x[i] = (b[i] - s) / a[i * n + i];
			q = fabs(x[i] - t) / (1.0 + fabs(x[i]));
			if (q > p)
				p = q;
		}
	}

	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：求解对称正定方程组的共轭梯度法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[n][n]				In			存放对称正定矩阵A
 *		double*				b[n]				In			存放方程组右端的常数向量
 *		int					n					In			方程组的阶数
 *		double				eps					In			给定的精度要求
 *		double*				x[n]				Out			返回方程组的解向量
 *
 *  返回：无
 *
**********************************************************************************************************************/
/*  求解对称正定方程组的共轭梯度法  */
void ISL_grad(double *a, int n, double *b, double eps, double *x)
{
	int i, k;
	double *p = NULL, *r = NULL, *s = NULL, *q = NULL, alpha, beta, d, e;

	p = new double[n];
	r = new double[n];
	s = new double[n];
	q = new double[n];

	for (i = 0; i <= n - 1; i++) {
		x[i] = 0.0;
		p[i] = b[i];
		r[i] = b[i];
	}

	i = 0;
	while (i <= n - 1) {
		ISL_trmul(a, p, n, n, 1, s);

		d = 0.0;
		e = 0.0;
		for (k = 0; k <= n - 1; k++) {
			d = d + p[k] * b[k];
			e = e + p[k] * s[k];
		}

		alpha = d / e;

		for (k = 0; k <= n - 1; k++)
			x[k] = x[k] + alpha * p[k];

		ISL_trmul(a, x, n, n, 1, q);

		d = 0.0;
		for (k = 0; k <= n - 1; k++) {
			r[k] = b[k] - q[k];
			d = d + r[k] * s[k];
		}
		beta = d / e;
		d = 0.0;

		for (k = 0; k <= n - 1; k++)
			d = d + r[k] * r[k];

		d = sqrt(d);

		if (d < eps) {
			if(p){ delete []p; p = NULL; }
			if(r){ delete []r; r = NULL; }
			if(s){ delete []s; s = NULL; }
			if(q){ delete []q; q = NULL; }
			return;
		}

		for (k = 0; k <= n - 1; k++)
			p[k] = r[k] - beta * p[k];

		i = i + 1;
	}
	if(p){ delete []p; p = NULL; }
	if(r){ delete []r; r = NULL; }
	if(s){ delete []s; s = NULL; }
	if(q){ delete []q; q = NULL; }

	return;
}


/**********************************************************************************************************************
 *
 *  功能：求解线性最小二乘问题的豪斯荷尔德变换法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[m][n]				In			存放超定方程组的系数矩阵A。返回时存放QR分解式中的R矩阵
 *		int					m					In			系数矩阵A的行数，m≥n
 *		int					n					In			系数矩阵A的列数，n≤m
 *		double*				b[m]				In			存放方程组右端的常数向量。返回时前n个分量存放方程组的最小二乘解
 *		double*				q[m][m]				Out			返回时存放QR分解式中的正交矩阵Q
 *
 *  返回：函数返回标志值。若返回的标志值为0，则表示程序工作失败（如A列线性相关）；
 *  				若返回的标志值不为0，则表示止常返回 。
 *
**********************************************************************************************************************/
/*  求解线性最小二乘问题的豪斯荷尔德变换法  */
int ISL_gmqr(double *a, int m, int n, double *b, double *q)
{
	int i, j;
	double d, *c = NULL;
	c = new double[n];

	i = ISL_maqr(a, m, n, q);
	if (i == 0) {
		if(c){ delete []c; c = NULL; }
		return (0);
	}

	for (i = 0; i <= n - 1; i++) {
		d = 0.0;
		for (j = 0; j <= m - 1; j++)
			d = d + q[j * m + i] * b[j];
		c[i] = d;
	}

	b[n - 1] = c[n - 1] / a[n * n - 1];
	for (i = n - 2; i >= 0; i--) {
		d = 0.0;
		for (j = i + 1; j <= n - 1; j++)
			d = d + a[i * n + j] * b[j];
		b[i] = (c[i] - d) / a[i * n + i];
	}

	if(c){ delete []c; c = NULL; }

	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：求解线性最小二乘问题的广义逆法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[m][n]				In			存放超定方程组的系数矩阵A。返回时其对角线依次给出奇异值，其余元素为0
 *		int					m					In			系数矩阵A的行数
 *		int					n					In			系数矩阵A的列数
 *		double*				b[m]				In			存放方程组右端的常数向量。
 *		double*				x[n]				Out			返回超定方程组的最小二乘解
 *		double*				aa[n][m]			Out			返回系数矩阵A的广义逆 A'
 *		double				eps					In			奇异值分解中的控制精度要求
 *		double*				u[m][m]				Out			返回系数矩阵A的奇异值分解式中的左奇异向量U
 *		double*				v[n][n]				Out			返回系数矩阵A的奇异值分解式中的右奇异向量Vt
 *		int					ka					In			ka = max(m,n)+1
 *
 *  返回：函数返回标志值。函数返回标志值。若返网的标志值小于0，则表示程序工作失败；
 *  				若返回的标志值大于0，则表示正常返回
 *
**********************************************************************************************************************/
/*  求解线性最小二乘问题的广义逆法  */
int ISL_gmiv(double *a, int m, int n, double *b, double *x, double *aa, double eps,
		double *u, double *v, int ka)
{
	int i, j;
	i = ISL_ginv(a, m, n, aa, eps, u, v, ka);
	if (i < 0)
		return (-1);
	for (i = 0; i <= n - 1; i++) {
		x[i] = 0.0;
		for (j = 0; j <= m - 1; j++)
			x[i] = x[i] + aa[i * m + j] * b[j];
	}
	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：求解病态方程组
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double*				a[m][n]				In			存放方程组的系数矩阵
 *		int					n					In			方程组的解数
 *		double*				b[m]				In			存放方程组右端的常数向量。
 *		double*				x[n]				Out			返回方程组的解向量
 *		double				eps					In			控制精度要求
 *
 *  返回：函数返回标志值。函数返回标志值。若返回的标志值为O，则表示程序工作失败；
 *  			若返回的标志值不为0，则表示正常返回
 *
**********************************************************************************************************************/
/*  求解病态方程组  */
int ISL_bint(double *a, int n, double *b, double eps, double *x)
{
	int i, j, k, kk;
	double *p = NULL, *r = NULL, *e = NULL, q, qq;

	p = new double[n*n];
	r = new double[n];
	e = new double[n];

	i = 60;
	for (k = 0; k <= n - 1; k++)
		for (j = 0; j <= n - 1; j++)
			p[k * n + j] = a[k * n + j];
	for (k = 0; k <= n - 1; k++)
		x[k] = b[k];

	kk = ISL_gaus(p, x, n);
	if (kk == 0) {
		if(p){ delete []p; p = NULL; }
		if(r){ delete []r; r = NULL; }
		if(e){ delete []e; e = NULL; }
		return (kk);
	}

	q = 1.0 + eps;
	while (q >= eps) {
		if (i == 0) {
			if(p){ delete []p; p = NULL; }
			if(r){ delete []r; r = NULL; }
			if(e){ delete []e; e = NULL; }
			return (i);
		}
		i = i - 1;
		ISL_trmul(a, x, n, n, 1, e);

		for (k = 0; k <= n - 1; k++)
			r[k] = b[k] - e[k];

		for (k = 0; k <= n - 1; k++)
			for (j = 0; j <= n - 1; j++)
				p[k * n + j] = a[k * n + j];

		kk = ISL_gaus(p, r, n);
		if (kk == 0) {
			if(p){ delete []p; p = NULL; }
			if(r){ delete []r; r = NULL; }
			if(e){ delete []e; e = NULL; }
			return (kk);
		}

		q = 0.0;
		for (k = 0; k <= n - 1; k++) {
			qq = fabs(r[k]) / (1.0 + fabs(x[k] + r[k]));
			if (qq > q)
				q = qq;
		}

		for (k = 0; k <= n - 1; k++)
			x[k] = x[k] + r[k];
	}

	if(p){ delete []p; p = NULL; }
	if(r){ delete []r; r = NULL; }
	if(e){ delete []e; e = NULL; }

	return (1);
}

// ==============================
// ==============================
// ==============================

} /*End of ISLib*/
