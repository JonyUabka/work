/**
*	@file	ISL_Matrix.cpp
*	@brief	[Source file of Matrix Functions], 矩阵运算；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-06-03
*	@refer	SU CWP
*/

#include "ISL_Matrix.h"

namespace ISLib {

void ISL_ppp(double *a, double *e, double *s, double *v, int m, int n)
{
	int i, j, p, q;
	double d;
	if (m >= n)
		i = n;
	else
		i = m;
	for (j = 1; j <= i - 1; j++) {
		a[(j - 1) * n + j - 1] = s[j - 1];
		a[(j - 1) * n + j] = e[j - 1];
	}
	a[(i - 1) * n + i - 1] = s[i - 1];
	if (m < n)
		a[(i - 1) * n + i] = e[i - 1];
	for (i = 1; i <= n - 1; i++)
		for (j = i + 1; j <= n; j++) {
			p = (i - 1) * n + j - 1;
			q = (j - 1) * n + i - 1;
			d = v[p];
			v[p] = v[q];
			v[q] = d;
		}
	return;
}

void ISL_sss(double *fg, double *cs)
{
	double r, d;
	if ((fabs(fg[0]) + fabs(fg[1])) == 0.0) {
		cs[0] = 1.0;
		cs[1] = 0.0;
		d = 0.0;
	} else {
		d = sqrt(fg[0] * fg[0] + fg[1] * fg[1]);
		if (fabs(fg[0]) > fabs(fg[1])) {
			d = fabs(d);
			if (fg[0] < 0.0)
				d = -d;
		}
		if (fabs(fg[1]) >= fabs(fg[0])) {
			d = fabs(d);
			if (fg[1] < 0.0)
				d = -d;
		}
		cs[0] = fg[0] / d;
		cs[1] = fg[1] / d;
	}
	r = 1.0;
	if (fabs(fg[0]) > fabs(fg[1]))
		r = cs[1];
	else if (cs[0] != 0.0)
		r = 1.0 / cs[0];
	fg[0] = d;
	fg[1] = r;
	return;
}
// ===================================================================================================================

/**********************************************************************************************************************
 *
 *  功能：实矩阵相乘
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[m][n]				In			存放矩阵A的元素
 *		double *			b[n][k]				In			存放矩阵B的元素
 *		int					m					In			矩阵A与乘积矩阵C的行数
 *		int					n					In			矩阵A的列数，矩阵B的行数
 *		int					k					In			矩阵B与乘积矩阵C的行数
 *		double *			c[m][k]				In/Out		返回乘积矩阵C=AB的元素，进入函数前请分配好内存空间
 *
 *  返回：无
 *
**********************************************************************************************************************/
void ISL_trmul(double *a, double *b, int m, int n, int k, double *c)
{
	int i, j, l, u;
	for (i = 0; i <= m - 1; i++) {
		for (j = 0; j <= k - 1; j++) {
			u = i * k + j;
			c[u] = 0.0;
			for (l = 0; l <= n - 1; l++)
				c[u] = c[u] + a[i * n + l] * b[l * k + j];
		}
	}
	return;
}


/**********************************************************************************************************************
 *
 *  功能：复矩阵相乘
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			ar[m][n]			In			存放矩阵A的实部与虚部元素
 *		double *			ai[m][n]			In
 *		double *			br[n][k]			In			存放矩阵B的实部与虚部元素
 *		double *			bi[n][k]			In
 *		int					m					In			矩阵A与乘积矩阵C的行数
 *		int					n					In			矩阵A的列数，矩阵B的行数
 *		int					k					In			矩阵B与乘积矩阵C的行数
 *		double *			cr[m][k]			In/Out		返回乘积矩阵C=AB的实部与虚部元素， 进入函数前请分配好内存空间
 *		double *			ci[m][k]			In/Out
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_tcmul(double *ar, double *ai, double *br, double *bi,int m, int n, int k, double *cr, double *ci)
{
	int i, j, l, u, v, w;
	double p, q, s;

	for (i = 0; i <= m - 1; i++) {
		for (j = 0; j <= k - 1; j++) {
			u = i * k + j;
			cr[u] = 0.0;
			ci[u] = 0.0;

			for (l = 0; l <= n - 1; l++) {
				v = i * n + l;
				w = l * k + j;
				p = ar[v] * br[w];
				q = ai[v] * bi[w];
				s = (ar[v] + ai[v]) * (br[w] + bi[w]);
				cr[u] = cr[u] + p - q;
				ci[u] = ci[u] + s - p - q;
			}
		}
	}
	return;
}


/**********************************************************************************************************************
 *
 *  功能：一般实矩阵 求逆
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[n][n]				In			存放矩阵A。返回时存放其逆矩阵A^-1
 *		int					n					In			矩阵阶数
 *
 *  返回：函数返回整型标志位，0表示A奇异；否则表示正常返回
 *
**********************************************************************************************************************/
int ISL_rinv(double *a, int n)
{
	int i, j, k, l, u, v;
	double d, p;

	int * is = new int[n];
	int * js = new int[n];

	for (k = 0; k <= n - 1; k++) {
		d = 0.0;
		for (i = k; i <= n - 1; i++) {
			for (j = k; j <= n - 1; j++) {
				l = i * n + j;
				p = fabs(a[l]);
				if (p > d) {
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		}

		if (d + 1.0 == 1.0) {
			if (is) {
				delete[] is;
				is = NULL;
			}
			if (js) {
				delete[] js;
				js = NULL;
			}

			printf("error ** not inv \n");
			return 0;
		}

		if (is[k] != k) {
			for (j = 0; j <= n - 1; j++) {
				u = k * n + j;
				v = is[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		}

		if (js[k] != k) {
			for (i = 0; i <= n - 1; i++) {
				u = i * n + k;
				v = i * n + js[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		}

		l = k * n + k;
		a[l] = 1.0 / a[l];

		for (j = 0; j <= n - 1; j++) {
			if (j != k) {
				u = k * n + j;
				a[u] = a[u] * a[l];
			}
		}

		for (i = 0; i <= n - 1; i++) {
			if (i != k) {
				for (j = 0; j <= n - 1; j++) {
					if (j != k) {
						u = i * n + j;
						a[u] = a[u] - a[i * n + k] * a[k * n + j];
					}
				}
			}
		}

		for (i = 0; i <= n - 1; i++) {
			if (i != k) {
				u = i * n + k;
				a[u] = -a[u] * a[l];
			}
		}
	}

	for (k = n - 1; k >= 0; k--) {
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++) {
				u = k * n + j;
				v = js[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; i++) {
				u = i * n + k;
				v = i * n + is[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
	}

	if (is) {
		delete[] is;
		is = NULL;
	}
	if (js) {
		delete[] js;
		js = NULL;
	}

	return 1;
}


/**********************************************************************************************************************
 *
 *  功能：一般复矩阵求逆
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			ar[n][n]			In			存放矩阵A的实部与虚部。返回时存放其逆矩阵A^-1的实部与虚部
 *		double *			ai[n][n]
 *
 *		int					n					In			矩阵阶数
 *
 *  返回：函数返回整型标志位，0表示A奇异；否则表示正常返回
 *
**********************************************************************************************************************/
int ISL_cinv(double *ar, double *ai, int n)
{
	int i, j, k, l, u, v, w;
	double p, q, s, t, d, b;

	int * is = new int[n];
	int * js = new int[n];

	for (k = 0; k <= n - 1; k++) {
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++) {
				u = i * n + j;
				p = ar[u] * ar[u] + ai[u] * ai[u];
				if (p > d) {
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}

		if (d + 1.0 == 1.0) {
			if (is) {
				delete[] is;
				is = NULL;
			}
			if (js) {
				delete[] js;
				js = NULL;
			}
			printf("error**not inv\n");
			return 0;
		}

		if (is[k] != k)
			for (j = 0; j <= n - 1; j++) {
				u = k * n + j;
				v = is[k] * n + j;
				t = ar[u];
				ar[u] = ar[v];
				ar[v] = t;
				t = ai[u];
				ai[u] = ai[v];
				ai[v] = t;
			}

		if (js[k] != k)
			for (i = 0; i <= n - 1; i++) {
				u = i * n + k;
				v = i * n + js[k];
				t = ar[u];
				ar[u] = ar[v];
				ar[v] = t;
				t = ai[u];
				ai[u] = ai[v];
				ai[v] = t;
			}

		l = k * n + k;
		ar[l] = ar[l] / d;
		ai[l] = -ai[l] / d;

		for (j = 0; j <= n - 1; j++)
			if (j != k) {
				u = k * n + j;
				p = ar[u] * ar[l];
				q = ai[u] * ai[l];
				s = (ar[u] + ai[u]) * (ar[l] + ai[l]);
				ar[u] = p - q;
				ai[u] = s - p - q;
			}

		for (i = 0; i <= n - 1; i++)
			if (i != k) {
				v = i * n + k;
				for (j = 0; j <= n - 1; j++)
					if (j != k) {
						u = k * n + j;
						w = i * n + j;
						p = ar[u] * ar[v];
						q = ai[u] * ai[v];
						s = (ar[u] + ai[u]) * (ar[v] + ai[v]);
						t = p - q;
						b = s - p - q;
						ar[w] = ar[w] - t;
						ai[w] = ai[w] - b;
					}
			}

		for (i = 0; i <= n - 1; i++)
			if (i != k) {
				u = i * n + k;
				p = ar[u] * ar[l];
				q = ai[u] * ai[l];
				s = (ar[u] + ai[u]) * (ar[l] + ai[l]);
				ar[u] = q - p;
				ai[u] = p + q - s;
			}
	}

	for (k = n - 1; k >= 0; k--) {
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++) {
				u = k * n + j;
				v = js[k] * n + j;
				t = ar[u];
				ar[u] = ar[v];
				ar[v] = t;
				t = ai[u];
				ai[u] = ai[v];
				ai[v] = t;
			}

		if (is[k] != k)
			for (i = 0; i <= n - 1; i++) {
				u = i * n + k;
				v = i * n + is[k];
				t = ar[u];
				ar[u] = ar[v];
				ar[v] = t;
				t = ai[u];
				ai[u] = ai[v];
				ai[v] = t;
			}
	}

	if (is) {
		delete[] is;
		is = NULL;
	}
	if (js) {
		delete[] js;
		js = NULL;
	}

	return 1;
}


/**********************************************************************************************************************
 *
 *  功能：对称正定矩阵的求逆
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[n][n]				In			存放 对称正定矩阵A。返回时存放其逆矩阵A^-1
 *		int					n					In			矩阵阶数
 *
 *  返回：函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；
 *  				若返回标志值大于0；否则表示正常返回。
 *
**********************************************************************************************************************/
int ISL_ssgj(double *a, int n)
{
	int i, j, k, m;
	double w, g;
	double * b = new double[n];

	for (k = 0; k <= n - 1; k++) {
		w = a[0];
		if (fabs(w) + 1.0 == 1.0) {
			if (b) {
				delete[] b;
				b = NULL;
			}
			printf("fail\n");
			return (-2);
		}
		m = n - k - 1;

		for (i = 1; i <= n - 1; i++) {
			g = a[i * n];
			b[i] = g / w;
			if (i <= m)
				b[i] = -b[i];

			for (j = 1; j <= i; j++)
				a[(i - 1) * n + j - 1] = a[i * n + j] + g * b[j];
		}

		a[n * n - 1] = 1.0 / w;
		for (i = 1; i <= n - 1; i++)
			a[(n - 1) * n + i - 1] = b[i];
	}

	for (i = 0; i <= n - 2; i++)
		for (j = i + 1; j <= n - 1; j++)
			a[i * n + j] = a[j * n + i];

	if (b) {
		delete[] b;
		b = NULL;
	}

	return 2;
}


/**********************************************************************************************************************
 *
 *  功能：托伯利兹矩阵求逆的特兰持方法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			t[n]				In			存放 T型矩阵中的元素。t0, t1, t2,..., t(n-1)
 *		double *			tt[n]				In			后n-1个元素存放 T型矩阵中的元素
 *		int					n					In			T型矩阵阶数
 *		double *			b[n][n]				In/Out		返回T型矩阵的逆矩阵
 *
 *  返回：函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；
 *  				若返回标志值大于0；否则表示正常返回。
 *
**********************************************************************************************************************/
int ISL_trch(double *t, double *tt, int n, double *b)
{
	int i, j, k;
	double a, s, *c, *r, *p;
	c = new double[n];
	r = new double[n];
	p = new double[n];

	if (fabs(t[0]) + 1.0 == 1.0) {
		if (c) {
			delete[] c;
			c = NULL;
		}
		if (r) {
			delete[] r;
			r = NULL;
		}
		if (p) {
			delete[] p;
			p = NULL;
		}
		printf("fail\n");
		return (-1);
	}

	a = t[0];
	c[0] = tt[1] / t[0];
	r[0] = t[1] / t[0];

	for (k = 0; k <= n - 3; k++) {
		s = 0.0;
		for (j = 1; j <= k + 1; j++)
			s = s + c[k + 1 - j] * tt[j];

		s = (s - tt[k + 2]) / a;
		for (i = 0; i <= k; i++)
			p[i] = c[i] + s * r[k - i];

		c[k + 1] = -s;
		s = 0.0;
		for (j = 1; j <= k + 1; j++)
			s = s + r[k + 1 - j] * t[j];

		s = (s - t[k + 2]) / a;
		for (i = 0; i <= k; i++) {
			r[i] = r[i] + s * c[k - i];
			c[k - i] = p[k - i];
		}
		r[k + 1] = -s;
		a = 0.0;
		for (j = 1; j <= k + 2; j++)
			a = a + t[j] * c[j - 1];

		a = t[0] - a;
		if (fabs(a) + 1.0 == 1.0) {
			if (c) {
				delete[] c;
				c = NULL;
			}
			if (r) {
				delete[] r;
				r = NULL;
			}
			if (p) {
				delete[] p;
				p = NULL;
			}
			printf("fail\n");
			return (-1);
		}
	}

	b[0] = 1.0 / a;
	for (i = 0; i <= n - 2; i++) {
		k = i + 1;
		j = (i + 1) * n;
		b[k] = -r[i] / a;
		b[j] = -c[i] / a;
	}
	for (i = 0; i <= n - 2; i++)
		for (j = 0; j <= n - 2; j++) {
			k = (i + 1) * n + j + 1;
			b[k] = b[i * n + j] - c[i] * b[j + 1];
			b[k] = b[k] + c[n - j - 2] * b[n - i - 1];
		}

	if (c) {
		delete[] c;
		c = NULL;
	}
	if (r) {
		delete[] r;
		r = NULL;
	}
	if (p) {
		delete[] p;
		p = NULL;
	}

	return 1;
}


/**********************************************************************************************************************
 *
 *  功能：求一般行列式的值，用全选主元高斯（Gauss）消去法计算n阶方阵A所对应的行列式值
 *
 *  说明：用全选主元高斯（Gauss）消去法对方阵A进行一系列变换使之成为上三角矩阵，其对角线上的个元素乘积即为行列式值。
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[n][n]				In			存放方阵A的元素，返回时被破坏
 *		int					n					In			方阵的阶数
 *
 *  返回：double, 函数返回行列式值。
 *
**********************************************************************************************************************/
double ISL_sdet(double *a, int n)
{
	int i, j, k, is=0, js=0, l, u, v;
	double f, det, q, d;

	f = 1.0;
	det = 1.0;

	for (k = 0; k <= n - 2; k++) {
		q = 0.0;

		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++) {
				l = i * n + j;
				d = fabs(a[l]);
				if (d > q) {
					q = d;
					is = i;
					js = j;
				}
			}

		if (q + 1.0 == 1.0) {
			det = 0.0;
			return (det);
		}

		if (is != k) {
			f = -f;
			for (j = k; j <= n - 1; j++) {
				u = k * n + j;
				v = is * n + j;
				d = a[u];
				a[u] = a[v];
				a[v] = d;
			}
		}

		if (js != k) {
			f = -f;
			for (i = k; i <= n - 1; i++) {
				u = i * n + js;
				v = i * n + k;
				d = a[u];
				a[u] = a[v];
				a[v] = d;
			}
		}

		l = k * n + k;
		det = det * a[l];
		for (i = k + 1; i <= n - 1; i++) {
			d = a[i * n + k] / a[l];
			for (j = k + 1; j <= n - 1; j++) {
				u = i * n + j;
				a[u] = a[u] - d * a[k * n + j];
			}
		}
	}

	det = f * det * a[n * n - 1];

	return (det);
}


/**********************************************************************************************************************
 *
 *  功能：求矩阵的秩，用全选主元高斯（Gauss）消去法计算矩阵的秩
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[m][n]				In			存放mxn阶矩阵A的元素，返回时被破坏
 *		int					m					In			矩阵A的行数
 *		int					n					In			矩阵A的列数
 *
 *  返回：int, 函数返回矩阵的秩。
 *
**********************************************************************************************************************/
int ISL_rank(double *a, int m, int n)
{
	int i, j, k, nn, is=0, js=0, l, ll, u, v;
	double q, d;

	nn = m;
	if (m >= n)
		nn = n;

	k = 0;
	for (l = 0; l <= nn - 1; l++) {
		q = 0.0;
		for (i = l; i <= m - 1; i++)
			for (j = l; j <= n - 1; j++) {
				ll = i * n + j;
				d = fabs(a[ll]);
				if (d > q) {
					q = d;
					is = i;
					js = j;
				}
			}

		if (q + 1.0 == 1.0)
			return (k);

		k = k + 1;
		if (is != l) {
			for (j = l; j <= n - 1; j++) {
				u = l * n + j;
				v = is * n + j;
				d = a[u];
				a[u] = a[v];
				a[v] = d;
			}
		}

		if (js != l) {
			for (i = l; i <= m - 1; i++) {
				u = i * n + js;
				v = i * n + l;
				d = a[u];
				a[u] = a[v];
				a[v] = d;
			}
		}

		ll = l * n + l;
		for (i = l + 1; i <= n - 1; i++) {
			d = a[i * n + l] / a[ll];

			for (j = l + 1; j <= n - 1; j++) {
				u = i * n + j;
				a[u] = a[u] - d * a[l * n + j];
			}
		}
	}
	return (k);
}


/**********************************************************************************************************************
 *
 *  功能：对称正定矩阵的乔里斯基分解与行列式求值,用乔里斯基（Cholesky）分解法求对称正定矩阵的三角分解，并求行列式的值
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[n][n]				In			存放对称正定矩阵A，返回时其下三角部分存放分解得到的下三角阵L，其余元素均为0
 *		int					n					In			矩阵A的阶数
 *		double &			det					In/Out		指向行列式的值
 *
 *  返回：函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；
 *  				若返回标志值大于0；否则表示正常返回。
 *
**********************************************************************************************************************/
int ISL_chol(double *a, int n, double &det)
{
	int i, j, k, u, l;
	double d;
	if ((a[0] + 1.0 == 1.0) || (a[0] < 0.0)) {
		printf("fail\n");
		return (-2);
	}
	a[0] = sqrt(a[0]);
	d = a[0];

	for (i = 1; i <= n - 1; i++) {
		u = i * n;
		a[u] = a[u] / a[0];
	}

	for (j = 1; j <= n - 1; j++) {
		l = j * n + j;
		for (k = 0; k <= j - 1; k++) {
			u = j * n + k;
			a[l] = a[l] - a[u] * a[u];
		}

		if ((a[l] + 1.0 == 1.0) || (a[l] < 0.0)) {
			printf("fail\n");
			return (-2);
		}
		a[l] = sqrt(a[l]);
		d = d * a[l];

		for (i = j + 1; i <= n - 1; i++) {
			u = i * n + j;
			for (k = 0; k <= j - 1; k++)
				a[u] = a[u] - a[i * n + k] * a[j * n + k];
			a[u] = a[u] / a[l];
		}
	}
	det = d * d;

	for (i = 0; i <= n - 2; i++)
		for (j = i + 1; j <= n - 1; j++)
			a[i * n + j] = 0.0;

	return (2);
}


/**********************************************************************************************************************
 *
 *  功能：矩阵的三角分解
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[n][n]				In			存放n阶矩阵A，返回时存放Q矩阵
 *		int					n					In			矩阵阶数
 *		double *			l[n][n]				In			返回时存放下三角矩阵L
 *		double *			u[n][n]				In			返回时存放上三角矩阵U
 *
 *  返回：函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；
 *  				若返回标志值大于0；否则表示正常返回。
 *
**********************************************************************************************************************/
int ISL_lluu(double *a, int n, double *l, double *u)
{
	int i, j, k, w, v, ll;

	for (k = 0; k <= n - 2; k++) {
		ll = k * n + k;
		if (fabs(a[ll]) + 1.0 == 1.0) {
			printf("fail\n");
			return (0);
		}

		for (i = k + 1; i <= n - 1; i++) {
			w = i * n + k;
			a[w] = a[w] / a[ll];
		}
		for (i = k + 1; i <= n - 1; i++) {
			w = i * n + k;

			for (j = k + 1; j <= n - 1; j++) {
				v = i * n + j;
				a[v] = a[v] - a[w] * a[k * n + j];
			}
		}
	}
	for (i = 0; i <= n - 1; i++) {
		for (j = 0; j < i; j++) {
			w = i * n + j;
			l[w] = a[w];
			u[w] = 0.0;
		}
		w = i * n + i;
		l[w] = 1.0;
		u[w] = a[w];

		for (j = i + 1; j <= n - 1; j++) {
			w = i * n + j;
			l[w] = 0.0;
			u[w] = a[w];
		}
	}
	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：一般实矩阵的QR分解，用豪斯荷尔德（Householder）变换对一般mxn阶的实矩阵进行QR分解
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[m][n]				In			存放mxn的实矩阵A。返回时其右上三角部分存放QR分解中的上三角矩阵R
 *		int					m					In			实矩阵A的行数
 *		int					n					In			实矩阵A的列数
 *		double *			q[m][n]				In			返回QR分解中的正交矩阵Q
 *
 *  返回：函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；
 *  				若返回标志值大于0；否则表示正常返回。
 *
**********************************************************************************************************************/
int ISL_maqr(double *a, int m, int n, double *q)
{
	int i,j,k,l,nn,p,jj;
	double u,alpha,w,t;

	if (m<n){
		printf("fail\n");
		return(0);
	}

	for (i = 0; i <= m - 1; i++)
		for (j = 0; j <= m - 1; j++) {
			l = i * m + j;
			q[l] = 0.0;
			if (i == j)
				q[l] = 1.0;
		}
	nn = n;
	if (m == n)
		nn = m - 1;

	for (k = 0; k <= nn - 1; k++) {
		u = 0.0;
		l = k * n + k;
		for (i = k; i <= m - 1; i++) {
			w = fabs(a[i * n + k]);
			if (w > u)
				u = w;
		}

		alpha = 0.0;
		for (i = k; i <= m - 1; i++) {
			t = a[i * n + k] / u;
			alpha = alpha + t * t;
		}
		if (a[l] > 0.0)
			u = -u;

		alpha = u * sqrt(alpha);
		if (fabs(alpha) + 1.0 == 1.0) {
			printf("fail\n");
			return (0);
		}

		u = sqrt(2.0 * alpha * (alpha - a[l]));
		if ((u + 1.0) != 1.0) {
			a[l] = (a[l] - alpha) / u;
			for (i = k + 1; i <= m - 1; i++) {
				p = i * n + k;
				a[p] = a[p] / u;
			}

			for (j = 0; j <= m - 1; j++) {
				t = 0.0;
				for (jj = k; jj <= m - 1; jj++)
					t = t + a[jj * n + k] * q[jj * m + j];
				for (i = k; i <= m - 1; i++) {
					p = i * m + j;
					q[p] = q[p] - 2.0 * t * a[i * n + k];
				}
			}

			for (j = k + 1; j <= n - 1; j++) {
				t = 0.0;
				for (jj = k; jj <= m - 1; jj++)
					t = t + a[jj * n + k] * a[jj * n + j];
				for (i = k; i <= m - 1; i++) {
					p = i * n + j;
					a[p] = a[p] - 2.0 * t * a[i * n + k];
				}
			}

			a[l] = alpha;
			for (i = k + 1; i <= m - 1; i++)
				a[i * n + k] = 0.0;
		}
	}
	for (i = 0; i <= m - 2; i++)
		for (j = i + 1; j <= m - 1; j++) {
			p = i * m + j;
			l = j * m + i;
			t = q[p];
			q[p] = q[l];
			q[l] = t;
		}

	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：一般实矩阵的奇异值分解，用豪斯荷尔德（Householder）变换以及变形QR对一般实矩阵A进行奇异值分解
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[m][n]				In/Out		存放mxn的实矩阵A。返回时其对角线给出奇异值（以非递增次序排列）
 *		int					m					In			实矩阵A的行数
 *		int					n					In			实矩阵A的列数
 *		double *			u[m][m]				In/Out		返回左奇异向量U
 *		double *			v[n][n]				In/Out		返回右奇异向量VT
 *		double				eps					In			给定的精度要求
 *		int					ka					In			其值为max(m, n)+1
 *
 *  返回：函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；
 *  				若返回标志值大于0；否则表示正常返回。
 *
**********************************************************************************************************************/
int ISL_muav(double *a, int m, int n, double *u, double *v, double eps, int ka)
{
	int i, j, k, l, it, ll, kk, ix, iy, mm, nn, iz, m1, ks;
	double d, dd, t, sm, sm1, em1, sk, ek, b, c, shh, fg[2], cs[2];
	double *s, *e, *w;

	s = new double[ka];
	e = new double[ka];
	w = new double[ka];

	it = 60;
	k = n;
	if (m - 1 < n)
		k = m - 1;
	l = m;
	if (n - 2 < m)
		l = n - 2;
	if (l < 0)
		l = 0;

	ll = k;
	if (l > k)
		ll = l;

	if (ll >= 1) {
		for (kk = 1; kk <= ll; kk++) {
			if (kk <= k) {
				d = 0.0;
				for (i = kk; i <= m; i++) {
					ix = (i - 1) * n + kk - 1;
					d = d + a[ix] * a[ix];
				}
				s[kk - 1] = sqrt(d);
				if (s[kk - 1] != 0.0) {
					ix = (kk - 1) * n + kk - 1;
					if (a[ix] != 0.0) {
						s[kk - 1] = fabs(s[kk - 1]);
						if (a[ix] < 0.0)
							s[kk - 1] = -s[kk - 1];
					}
					for (i = kk; i <= m; i++) {
						iy = (i - 1) * n + kk - 1;
						a[iy] = a[iy] / s[kk - 1];
					}
					a[ix] = 1.0 + a[ix];
				}
				s[kk - 1] = -s[kk - 1];
			}
			if (n >= kk + 1) {
				for (j = kk + 1; j <= n; j++) {
					if ((kk <= k) && (s[kk - 1] != 0.0)) {
						d = 0.0;
						for (i = kk; i <= m; i++) {
							ix = (i - 1) * n + kk - 1;
							iy = (i - 1) * n + j - 1;
							d = d + a[ix] * a[iy];
						}
						d = -d / a[(kk - 1) * n + kk - 1];
						for (i = kk; i <= m; i++) {
							ix = (i - 1) * n + j - 1;
							iy = (i - 1) * n + kk - 1;
							a[ix] = a[ix] + d * a[iy];
						}
					}
					e[j - 1] = a[(kk - 1) * n + j - 1];
				}
			}
			if (kk <= k) {
				for (i = kk; i <= m; i++) {
					ix = (i - 1) * m + kk - 1;
					iy = (i - 1) * n + kk - 1;
					u[ix] = a[iy];
				}
			}
			if (kk <= l) {
				d = 0.0;
				for (i = kk + 1; i <= n; i++)
					d = d + e[i - 1] * e[i - 1];
				e[kk - 1] = sqrt(d);
				if (e[kk - 1] != 0.0) {
					if (e[kk] != 0.0) {
						e[kk - 1] = fabs(e[kk - 1]);
						if (e[kk] < 0.0)
							e[kk - 1] = -e[kk - 1];
					}
					for (i = kk + 1; i <= n; i++)
						e[i - 1] = e[i - 1] / e[kk - 1];
					e[kk] = 1.0 + e[kk];
				}
				e[kk - 1] = -e[kk - 1];
				if ((kk + 1 <= m) && (e[kk - 1] != 0.0)) {
					for (i = kk + 1; i <= m; i++)
						w[i - 1] = 0.0;
					for (j = kk + 1; j <= n; j++)
						for (i = kk + 1; i <= m; i++)
							w[i - 1] = w[i - 1] + e[j - 1] * a[(i - 1) * n + j
									- 1];
					for (j = kk + 1; j <= n; j++)
						for (i = kk + 1; i <= m; i++) {
							ix = (i - 1) * n + j - 1;
							a[ix] = a[ix] - w[i - 1] * e[j - 1] / e[kk];
						}
				}
				for (i = kk + 1; i <= n; i++)
					v[(i - 1) * n + kk - 1] = e[i - 1];
			}
		}
	}
	mm = n;
	if (m + 1 < n)
		mm = m + 1;
	if (k < n)
		s[k] = a[k * n + k];
	if (m < mm)
		s[mm - 1] = 0.0;
	if (l + 1 < mm)
		e[l] = a[l * n + mm - 1];
	e[mm - 1] = 0.0;
	nn = m;
	if (m > n)
		nn = n;
	if (nn >= k + 1) {
		for (j = k + 1; j <= nn; j++) {
			for (i = 1; i <= m; i++)
				u[(i - 1) * m + j - 1] = 0.0;
			u[(j - 1) * m + j - 1] = 1.0;
		}
	}
	if (k >= 1) {
		for (ll = 1; ll <= k; ll++) {
			kk = k - ll + 1;
			iz = (kk - 1) * m + kk - 1;
			if (s[kk - 1] != 0.0) {
				if (nn >= kk + 1)
					for (j = kk + 1; j <= nn; j++) {
						d = 0.0;
						for (i = kk; i <= m; i++) {
							ix = (i - 1) * m + kk - 1;
							iy = (i - 1) * m + j - 1;
							d = d + u[ix] * u[iy] / u[iz];
						}
						d = -d;
						for (i = kk; i <= m; i++) {
							ix = (i - 1) * m + j - 1;
							iy = (i - 1) * m + kk - 1;
							u[ix] = u[ix] + d * u[iy];
						}
					}
				for (i = kk; i <= m; i++) {
					ix = (i - 1) * m + kk - 1;
					u[ix] = -u[ix];
				}
				u[iz] = 1.0 + u[iz];
				if (kk - 1 >= 1)
					for (i = 1; i <= kk - 1; i++)
						u[(i - 1) * m + kk - 1] = 0.0;
			} else {
				for (i = 1; i <= m; i++)
					u[(i - 1) * m + kk - 1] = 0.0;
				u[(kk - 1) * m + kk - 1] = 1.0;
			}
		}
	}
	for (ll = 1; ll <= n; ll++) {
		kk = n - ll + 1;
		iz = kk * n + kk - 1;
		if ((kk <= l) && (e[kk - 1] != 0.0)) {
			for (j = kk + 1; j <= n; j++) {
				d = 0.0;
				for (i = kk + 1; i <= n; i++) {
					ix = (i - 1) * n + kk - 1;
					iy = (i - 1) * n + j - 1;
					d = d + v[ix] * v[iy] / v[iz];
				}
				d = -d;
				for (i = kk + 1; i <= n; i++) {
					ix = (i - 1) * n + j - 1;
					iy = (i - 1) * n + kk - 1;
					v[ix] = v[ix] + d * v[iy];
				}
			}
		}
		for (i = 1; i <= n; i++)
			v[(i - 1) * n + kk - 1] = 0.0;
		v[iz - n] = 1.0;
	}
	for (i = 1; i <= m; i++)
		for (j = 1; j <= n; j++)
			a[(i - 1) * n + j - 1] = 0.0;
	m1 = mm;
	it = 60;
	while (1 == 1) {
		if (mm == 0) {
			ISL_ppp(a, e, s, v, m, n);
			if(s){ delete []s; s = NULL; }
			if(e){ delete []e; e = NULL; }
			if(w){ delete []w; w = NULL; }
			return (1);
		}
		if (it == 0) {
			ISL_ppp(a, e, s, v, m, n);
			if(s){ delete []s; s = NULL; }
			if(e){ delete []e; e = NULL; }
			if(w){ delete []w; w = NULL; }
			return (-1);
		}
		kk = mm - 1;
		while ((kk != 0) && (fabs(e[kk - 1]) != 0.0)) {
			d = fabs(s[kk - 1]) + fabs(s[kk]);
			dd = fabs(e[kk - 1]);
			if (dd > eps * d)
				kk = kk - 1;
			else
				e[kk - 1] = 0.0;
		}
		if (kk == mm - 1) {
			kk = kk + 1;
			if (s[kk - 1] < 0.0) {
				s[kk - 1] = -s[kk - 1];
				for (i = 1; i <= n; i++) {
					ix = (i - 1) * n + kk - 1;
					v[ix] = -v[ix];
				}
			}
			while ((kk != m1) && (s[kk - 1] < s[kk])) {
				d = s[kk - 1];
				s[kk - 1] = s[kk];
				s[kk] = d;
				if (kk < n)
					for (i = 1; i <= n; i++) {
						ix = (i - 1) * n + kk - 1;
						iy = (i - 1) * n + kk;
						d = v[ix];
						v[ix] = v[iy];
						v[iy] = d;
					}
				if (kk < m)
					for (i = 1; i <= m; i++) {
						ix = (i - 1) * m + kk - 1;
						iy = (i - 1) * m + kk;
						d = u[ix];
						u[ix] = u[iy];
						u[iy] = d;
					}
				kk = kk + 1;
			}
			it = 60;
			mm = mm - 1;
		} else {
			ks = mm;
			while ((ks > kk) && (fabs(s[ks - 1]) != 0.0)) {
				d = 0.0;
				if (ks != mm)
					d = d + fabs(e[ks - 1]);
				if (ks != kk + 1)
					d = d + fabs(e[ks - 2]);
				dd = fabs(s[ks - 1]);
				if (dd > eps * d)
					ks = ks - 1;
				else
					s[ks - 1] = 0.0;
			}
			if (ks == kk) {
				kk = kk + 1;
				d = fabs(s[mm - 1]);
				t = fabs(s[mm - 2]);
				if (t > d)
					d = t;
				t = fabs(e[mm - 2]);
				if (t > d)
					d = t;
				t = fabs(s[kk - 1]);
				if (t > d)
					d = t;
				t = fabs(e[kk - 1]);
				if (t > d)
					d = t;
				sm = s[mm - 1] / d;
				sm1 = s[mm - 2] / d;
				em1 = e[mm - 2] / d;
				sk = s[kk - 1] / d;
				ek = e[kk - 1] / d;
				b = ((sm1 + sm) * (sm1 - sm) + em1 * em1) / 2.0;
				c = sm * em1;
				c = c * c;
				shh = 0.0;
				if ((b != 0.0) || (c != 0.0)) {
					shh = sqrt(b * b + c);
					if (b < 0.0)
						shh = -shh;
					shh = c / (b + shh);
				}
				fg[0] = (sk + sm) * (sk - sm) - shh;
				fg[1] = sk * ek;
				for (i = kk; i <= mm - 1; i++) {
					ISL_sss(fg, cs);
					if (i != kk)
						e[i - 2] = fg[0];
					fg[0] = cs[0] * s[i - 1] + cs[1] * e[i - 1];
					e[i - 1] = cs[0] * e[i - 1] - cs[1] * s[i - 1];
					fg[1] = cs[1] * s[i];
					s[i] = cs[0] * s[i];
					if ((cs[0] != 1.0) || (cs[1] != 0.0))
						for (j = 1; j <= n; j++) {
							ix = (j - 1) * n + i - 1;
							iy = (j - 1) * n + i;
							d = cs[0] * v[ix] + cs[1] * v[iy];
							v[iy] = -cs[1] * v[ix] + cs[0] * v[iy];
							v[ix] = d;
						}
					ISL_sss(fg, cs);
					s[i - 1] = fg[0];
					fg[0] = cs[0] * e[i - 1] + cs[1] * s[i];
					s[i] = -cs[1] * e[i - 1] + cs[0] * s[i];
					fg[1] = cs[1] * e[i];
					e[i] = cs[0] * e[i];
					if (i < m)
						if ((cs[0] != 1.0) || (cs[1] != 0.0))
							for (j = 1; j <= m; j++) {
								ix = (j - 1) * m + i - 1;
								iy = (j - 1) * m + i;
								d = cs[0] * u[ix] + cs[1] * u[iy];
								u[iy] = -cs[1] * u[ix] + cs[0] * u[iy];
								u[ix] = d;
							}
				}
				e[mm - 2] = fg[0];
				it = it - 1;
			} else {
				if (ks == mm) {
					kk = kk + 1;
					fg[1] = e[mm - 2];
					e[mm - 2] = 0.0;
					for (ll = kk; ll <= mm - 1; ll++) {
						i = mm + kk - ll - 1;
						fg[0] = s[i - 1];
						ISL_sss(fg, cs);
						s[i - 1] = fg[0];
						if (i != kk) {
							fg[1] = -cs[1] * e[i - 2];
							e[i - 2] = cs[0] * e[i - 2];
						}
						if ((cs[0] != 1.0) || (cs[1] != 0.0))
							for (j = 1; j <= n; j++) {
								ix = (j - 1) * n + i - 1;
								iy = (j - 1) * n + mm - 1;
								d = cs[0] * v[ix] + cs[1] * v[iy];
								v[iy] = -cs[1] * v[ix] + cs[0] * v[iy];
								v[ix] = d;
							}
					}
				} else {
					kk = ks + 1;
					fg[1] = e[kk - 2];
					e[kk - 2] = 0.0;
					for (i = kk; i <= mm; i++) {
						fg[0] = s[i - 1];
						ISL_sss(fg, cs);
						s[i - 1] = fg[0];
						fg[1] = -cs[1] * e[i - 1];
						e[i - 1] = cs[0] * e[i - 1];
						if ((cs[0] != 1.0) || (cs[1] != 0.0))
							for (j = 1; j <= m; j++) {
								ix = (j - 1) * m + i - 1;
								iy = (j - 1) * m + kk - 2;
								d = cs[0] * u[ix] + cs[1] * u[iy];
								u[iy] = -cs[1] * u[ix] + cs[0] * u[iy];
								u[ix] = d;
							}
					}
				}
			}
		}
	}
	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：求广义逆的奇异值分解法，利用奇异值分解求一般mxn阶实矩阵A的广义逆A+
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[m][n]				In/Out		存放mxn的实矩阵A。返回时其对角线给出奇异值（以非递增次序排列），其余元素均为0
 *		int					m					In			实矩阵A的行数
 *		int					n					In			实矩阵A的列数
 *		double *			aa[m][m]			In/Out		返回A的广义逆A+
 *		double				eps					In			给定的精度要求
 *		double *			u[m][m]				In/Out		返回右奇异向量U
 *		double *			v[n][n]				In/Out		返回右奇异向量VT
 *		int					ka					In			其值为max(m, n)+1
 *
 *  返回：函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；
 *  				若返回标志值大于0；否则表示正常返回。
 *
**********************************************************************************************************************/
int ISL_ginv(double *a, int m, int n, double *aa, double eps, double *u, double *v, int ka)
{
	int i, j, k, l, t, p, q, f;
	i = ISL_muav(a, m, n, u, v, eps, ka);
	if (i < 0)
		return (-1);
	j = n;
	if (m < n)
		j = m;
	j = j - 1;
	k = 0;
	while ((k <= j) && (a[k * n + k] != 0.0))
		k = k + 1;
	k = k - 1;
	for (i = 0; i <= n - 1; i++)
		for (j = 0; j <= m - 1; j++) {
			t = i * m + j;
			aa[t] = 0.0;
			for (l = 0; l <= k; l++) {
				f = l * n + i;
				p = j * m + l;
				q = l * n + l;
				aa[t] = aa[t] + v[f] * u[p] / a[q];
			}
		}
	return (1);
}


// ============================
// ===	矩阵特征值运算与特征向量的计算	===
// ============================
/**********************************************************************************************************************
 *
 *  功能：约化对称矩阵为对称三角阵的豪斯赫尔德变换法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[n][n]				In			存放n阶实对称矩阵A
 *		int					n					In			实对称矩阵A的阶数
 *		double *			q[n][n]				Out			返回豪斯荷尔德变换的乘积矩阵Q。在与nj_sstq()联用时，若将Q矩阵作为函数nj_sstq()
 *															中的一个参数，则可以计算一般实对称矩阵的全部特征值及相应的特征向量
 *		double *			b[n]				Out			返回对称三角阵中的主对角线元素
 *		double *			c[n]				Out			前n-1个元素返回对称三角阵中的次对角线元素
 *
 *  返回：无
 *
**********************************************************************************************************************/
void ISL_strq(double *a, int n, double *q, double *b, double *c)
{
	int i, j, k, u;
	double h, f, g, h2;

	for (i = 0; i <= n - 1; i++)
		for (j = 0; j <= n - 1; j++) {
			u = i * n + j;
			q[u] = a[u];
		}

	for (i = n - 1; i >= 1; i--) {
		h = 0.0;
		if (i > 1)
			for (k = 0; k <= i - 1; k++) {
				u = i * n + k;
				h = h + q[u] * q[u];
			}
		if (h + 1.0 == 1.0) {
			c[i] = 0.0;
			if (i == 1)
				c[i] = q[i * n + i - 1];
			b[i] = 0.0;
		} else {
			c[i] = sqrt(h);
			u = i * n + i - 1;
			if (q[u] > 0.0)
				c[i] = -c[i];
			h = h - q[u] * c[i];
			q[u] = q[u] - c[i];
			f = 0.0;
			for (j = 0; j <= i - 1; j++) {
				q[j * n + i] = q[i * n + j] / h;
				g = 0.0;
				for (k = 0; k <= j; k++)
					g = g + q[j * n + k] * q[i * n + k];
				if (j + 1 <= i - 1)
					for (k = j + 1; k <= i - 1; k++)
						g = g + q[k * n + j] * q[i * n + k];
				c[j] = g / h;
				f = f + g * q[j * n + i];
			}
			h2 = f / (h + h);
			for (j = 0; j <= i - 1; j++) {
				f = q[i * n + j];
				g = c[j] - h2 * f;
				c[j] = g;
				for (k = 0; k <= j; k++) {
					u = j * n + k;
					q[u] = q[u] - f * c[k] - g * q[i * n + k];
				}
			}
			b[i] = h;
		}
	}

	for (i = 0; i <= n - 2; i++)
		c[i] = c[i + 1];

	c[n - 1] = 0.0;
	b[0] = 0.0;

	for (i = 0; i <= n - 1; i++) {
		if ((b[i] != 0.0) && (i - 1 >= 0))
			for (j = 0; j <= i - 1; j++) {
				g = 0.0;
				for (k = 0; k <= i - 1; k++)
					g = g + q[i * n + k] * q[k * n + j];
				for (k = 0; k <= i - 1; k++) {
					u = k * n + j;
					q[u] = q[u] - g * q[k * n + i];
				}
			}
		u = i * n + i;
		b[i] = q[u];
		q[u] = 1.0;
		if (i - 1 >= 0)
			for (j = 0; j <= i - 1; j++) {
				q[i * n + j] = 0.0;
				q[j * n + i] = 0.0;
			}
	}

	return;
}


/**********************************************************************************************************************
 *
 *  功能：求对称三对角阵的全部特征值与特征向量
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					n					In			实对称三角阵的阶数
 *		double *			b[n]				In/Out		存放n阶实对称三角阵的主对角线上的元素。返回时存放全部特征值
 *		double *			c[n]				In			前n-1个元素存放对称三角阵中的次对角线上的元素
 *
 *		double *			q[n][n]				Out			若存放n阶单位矩阵，则返回实对称三对角阵T的特征向量组；若存放由函数nj_strq（）所返回的一般
 *															实对称矩阵A的豪斯荷尔德的乘积矩阵Q，则返回实对称矩阵A的特征向量组。其中q中的第j列为数组b中第j个
 *															特征值对应的特征向量
 *
 *		double				eps					In			控制精度要求
 *		int					l					Out			允许的最大迭代次数
 *
 *  返回：int，函数返回标志值。若返回的标志值小于0，则表示程序工作失败；若返回的标志值大于0，则说明程序正常返回
 *
**********************************************************************************************************************/
int ISL_sstq(int n, double *b, double *c, double *q, double eps, int l)
{
	int i, j, k, m, it, u, v;
	double d, f, h, g, p, r, e, s;

	c[n - 1] = 0.0;
	d = 0.0;
	f = 0.0;

	for (j = 0; j <= n - 1; j++) {
		it = 0;
		h = eps * (fabs(b[j]) + fabs(c[j]));
		if (h > d)
			d = h;
		m = j;
		while ((m <= n - 1) && (fabs(c[m]) > d))
			m = m + 1;
		if (m != j) {
			do {
				if (it == l) {
					printf("fail\n");
					return (-1);
				}
				it = it + 1;
				g = b[j];
				p = (b[j + 1] - g) / (2.0 * c[j]);
				r = sqrt(p * p + 1.0);
				if (p >= 0.0)
					b[j] = c[j] / (p + r);
				else
					b[j] = c[j] / (p - r);
				h = g - b[j];
				for (i = j + 1; i <= n - 1; i++)
					b[i] = b[i] - h;
				f = f + h;
				p = b[m];
				e = 1.0;
				s = 0.0;
				for (i = m - 1; i >= j; i--) {
					g = e * c[i];
					h = e * p;
					if (fabs(p) >= fabs(c[i])) {
						e = c[i] / p;
						r = sqrt(e * e + 1.0);
						c[i + 1] = s * p * r;
						s = e / r;
						e = 1.0 / r;
					} else {
						e = p / c[i];
						r = sqrt(e * e + 1.0);
						c[i + 1] = s * c[i] * r;
						s = 1.0 / r;
						e = e / r;
					}
					p = e * b[i] - s * g;
					b[i + 1] = h + s * (e * g + s * b[i]);
					for (k = 0; k <= n - 1; k++) {
						u = k * n + i + 1;
						v = u - 1;
						h = q[u];
						q[u] = s * q[v] + e * h;
						q[v] = e * q[v] - s * h;
					}
				}
				c[j] = s * p;
				b[j] = e * p;
			} while (fabs(c[j]) > d);
		}
		b[j] = b[j] + f;
	}

	for (i = 0; i <= n - 1; i++) {
		k = i;
		p = b[i];
		if (i + 1 <= n - 1) {
			j = i + 1;
			while ((j <= n - 1) && (b[j] <= p)) {
				k = j;
				p = b[j];
				j = j + 1;
			}
		}

		if (k != i) {
			b[k] = b[i];
			b[i] = p;
			for (j = 0; j <= n - 1; j++) {
				u = j * n + i;
				v = j * n + k;
				p = q[u];
				q[u] = q[v];
				q[v] = p;
			}
		}
	}

	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：约化一般实矩阵为赫申伯格矩阵的初等相似变换法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[n][n]				In			存放一般实矩阵A，返回上H矩阵
 *		int					n					In			实对称矩阵A的阶数
 *
 *  返回：无
 *
**********************************************************************************************************************/
void ISL_hhbg(double *a, int n)
{
	int i=0, j, k, u, v;
	double d, t;

	for (k = 1; k <= n - 2; k++) {
		d = 0.0;

		for (j = k; j <= n - 1; j++) {
			u = j * n + k - 1;
			t = a[u];
			if (fabs(t) > fabs(d)) {
				d = t;
				i = j;
			}
		}

		if (fabs(d) + 1.0 != 1.0) {
			if (i != k) {
				for (j = k - 1; j <= n - 1; j++) {
					u = i * n + j;
					v = k * n + j;
					t = a[u];
					a[u] = a[v];
					a[v] = t;
				}

				for (j = 0; j <= n - 1; j++) {
					u = j * n + i;
					v = j * n + k;
					t = a[u];
					a[u] = a[v];
					a[v] = t;
				}
			}

			for (i = k + 1; i <= n - 1; i++) {
				u = i * n + k - 1;
				t = a[u] / d;
				a[u] = 0.0;
				for (j = k; j <= n - 1; j++) {
					v = i * n + j;
					a[v] = a[v] - t * a[k * n + j];
				}

				for (j = 0; j <= n - 1; j++) {
					v = j * n + k;
					a[v] = a[v] + t * a[j * n + i];
				}
			}
		}
	}
	return;
}


/**********************************************************************************************************************
 *
 *  功能：求赫申伯格矩阵全部特征值的QR方法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[n][n]				In			存在上H矩阵A
 *		int					n					In			上H矩阵A的阶数
 *		double *			u[n]				Out			返回n个特征值的实部
 *		double *			v[n]				Out			返回n个特征值的虚部
 *		double				eps					In			控制精度要求
 *		int					jt					In			允许的最大迭代次数
 *
 *  返回：int，函数返回标志值。若返回的标志值小于0，则表示程序工作失败；若返回的标志值大于0，则说明程序正常返回
 *
**********************************************************************************************************************/
int ISL_hhqr(double *a, int n, double *u, double *v, double eps, int jt)
{
	int m, it, i, j, k, l, ii, jj, kk, ll;
	double b, c, w, g, xy, p, q, r, x, s, e, f, z, y;
	it = 0;
	m = n;

	while (m != 0) {
		l = m - 1;
		while ((l > 0) && (fabs(a[l * n + l - 1]) > eps * (fabs(
				a[(l - 1) * n + l - 1]) + fabs(a[l * n + l]))))
			l = l - 1;

		ii = (m - 1) * n + m - 1;
		jj = (m - 1) * n + m - 2;
		kk = (m - 2) * n + m - 1;
		ll = (m - 2) * n + m - 2;

		if (l == m - 1) {
			u[m - 1] = a[(m - 1) * n + m - 1];
			v[m - 1] = 0.0;
			m = m - 1;
			it = 0;
		} else if (l == m - 2) {
			b = -(a[ii] + a[ll]);
			c = a[ii] * a[ll] - a[jj] * a[kk];
			w = b * b - 4.0 * c;
			y = sqrt(fabs(w));
			if (w > 0.0) {
				xy = 1.0;
				if (b < 0.0)
					xy = -1.0;
				u[m - 1] = (-b - xy * y) / 2.0;
				u[m - 2] = c / u[m - 1];
				v[m - 1] = 0.0;
				v[m - 2] = 0.0;
			} else {
				u[m - 1] = -b / 2.0;
				u[m - 2] = u[m - 1];
				v[m - 1] = y / 2.0;
				v[m - 2] = -v[m - 1];
			}

			m = m - 2;
			it = 0;

		} else {
			if (it >= jt) {
				printf("fail\n");
				return (-1);
			}

			it = it + 1;
			for (j = l + 2; j <= m - 1; j++)
				a[j * n + j - 2] = 0.0;

			for (j = l + 3; j <= m - 1; j++)
				a[j * n + j - 3] = 0.0;

			for (k = l; k <= m - 2; k++) {
				if (k != l) {
					p = a[k * n + k - 1];
					q = a[(k + 1) * n + k - 1];
					r = 0.0;
					if (k != m - 2)
						r = a[(k + 2) * n + k - 1];
				} else {
					x = a[ii] + a[ll];
					y = a[ll] * a[ii] - a[kk] * a[jj];
					ii = l * n + l;
					jj = l * n + l + 1;
					kk = (l + 1) * n + l;
					ll = (l + 1) * n + l + 1;
					p = a[ii] * (a[ii] - x) + a[jj] * a[kk] + y;
					q = a[kk] * (a[ii] + a[ll] - x);
					r = a[kk] * a[(l + 2) * n + l + 1];
				}

				if ((fabs(p) + fabs(q) + fabs(r)) != 0.0) {
					xy = 1.0;
					if (p < 0.0)
						xy = -1.0;
					s = xy * sqrt(p * p + q * q + r * r);
					if (k != l)
						a[k * n + k - 1] = -s;
					e = -q / s;
					f = -r / s;
					x = -p / s;
					y = -x - f * r / (p + s);
					g = e * r / (p + s);
					z = -x - e * q / (p + s);
					for (j = k; j <= m - 1; j++) {
						ii = k * n + j;
						jj = (k + 1) * n + j;
						p = x * a[ii] + e * a[jj];
						q = e * a[ii] + y * a[jj];
						r = f * a[ii] + g * a[jj];
						if (k != m - 2) {
							kk = (k + 2) * n + j;
							p = p + f * a[kk];
							q = q + g * a[kk];
							r = r + z * a[kk];
							a[kk] = r;
						}
						a[jj] = q;
						a[ii] = p;
					}

					j = k + 3;
					if (j >= m - 1)
						j = m - 1;

					for (i = l; i <= j; i++) {
						ii = i * n + k;
						jj = i * n + k + 1;
						p = x * a[ii] + e * a[jj];
						q = e * a[ii] + y * a[jj];
						r = f * a[ii] + g * a[jj];
						if (k != m - 2) {
							kk = i * n + k + 2;
							p = p + f * a[kk];
							q = q + g * a[kk];
							r = r + z * a[kk];
							a[kk] = r;
						}
						a[jj] = q;
						a[ii] = p;
					}
				}
			}
		}
	}
	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：求实对称矩阵特征值与特征向量的雅可比法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[n][n]				In			存放n阶实对称矩阵。返回时对角线上存放n个特征值
 *		int					n					In			实对称矩阵A的阶数
 *		double *			v[n][n]				Out			返回特征向量。其中第j列为与第j个特征值对应的特征向量
 *		double				eps					In			控制精度要求
 *		int					jt					In			控制最大迭代次数
 *
 *  返回：int，函数返回标志值。若返回的标志值小于0，则表示程序工作失败；若返回的标志值大于0，则说明程序正常返回
 *
**********************************************************************************************************************/
int ISL_jcbi(double *a, int n,double * v, double eps, int jt)
{
	int i, j, p=0, q=0, u, w, t, s, l;
	double fm, cn, sn, omega, x, y, d;
	l = 1;

	for (i = 0; i <= n - 1; i++) {
		v[i * n + i] = 1.0;
		for (j = 0; j <= n - 1; j++)
			if (i != j)
				v[i * n + j] = 0.0;
	}

	while (1 == 1) {
		fm = 0.0;
		for (i = 1; i <= n - 1; i++)
			for (j = 0; j <= i - 1; j++) {
				d = fabs(a[i * n + j]);
				if ((i != j) && (d > fm)) {
					fm = d;
					p = i;
					q = j;
				}
			}

		if (fm < eps)
			return (1);

		if (l > jt)
			return (-1);

		l = l + 1;
		u = p * n + q;
		w = p * n + p;
		t = q * n + p;
		s = q * n + q;
		x = -a[u];
		y = (a[s] - a[w]) / 2.0;
		omega = x / sqrt(x * x + y * y);

		if (y < 0.0)
			omega = -omega;

		sn = 1.0 + sqrt(1.0 - omega * omega);
		sn = omega / sqrt(2.0 * sn);
		cn = sqrt(1.0 - sn * sn);
		fm = a[w];
		a[w] = fm * cn * cn + a[s] * sn * sn + a[u] * omega;
		a[s] = fm * sn * sn + a[s] * cn * cn - a[u] * omega;
		a[u] = 0.0;
		a[t] = 0.0;

		for (j = 0; j <= n - 1; j++)
			if ((j != p) && (j != q)) {
				u = p * n + j;
				w = q * n + j;
				fm = a[u];
				a[u] = fm * cn + a[w] * sn;
				a[w] = -fm * sn + a[w] * cn;
			}

		for (i = 0; i <= n - 1; i++)
			if ((i != p) && (i != q)) {
				u = i * n + p;
				w = i * n + q;
				fm = a[u];
				a[u] = fm * cn + a[w] * sn;
				a[w] = -fm * sn + a[w] * cn;
			}

		for (i = 0; i <= n - 1; i++) {
			u = i * n + p;
			w = i * n + q;
			fm = v[u];
			v[u] = fm * cn + v[w] * sn;
			v[w] = -fm * sn + v[w] * cn;
		}
	}
	return (1);
}


/**********************************************************************************************************************
 *
 *  功能：求实对称矩阵特征值与特征向量的雅可比过关法
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		double *			a[n][n]				In			存放n阶实对称矩阵。返回时对角线上存放n个特征值
 *		int					n					In			实对称矩阵A的阶数
 *		double *			v[n][n]				Out			返回特征向量。其中第j列为与第j个特征值对应的特征向量
 *		double				eps					In			控制精度要求
 *
 *  返回：int，函数返回标志值。若返回的标志值小于0，则表示程序工作失败；若返回的标志值大于0，则说明程序正常返回
 *
**********************************************************************************************************************/
void ISL_jcbj(double *a, int n, double *v, double eps)
{
	int i, j, p, q, u, w, t, s;
	double ff, fm, cn, sn, omega, x, y, d;

	for (i = 0; i <= n - 1; i++) {
		v[i * n + i] = 1.0;
		for (j = 0; j <= n - 1; j++)
			if (i != j)
				v[i * n + j] = 0.0;
	}

	ff = 0.0;
	for (i = 1; i <= n - 1; i++)
		for (j = 0; j <= i - 1; j++) {
			d = a[i * n + j];
			ff = ff + d * d;
		}

	ff = sqrt(2.0 * ff);

	loop0: ff = ff / (1.0 * n);

	loop1: for (i = 1; i <= n - 1; i++)
		for (j = 0; j <= i - 1; j++) {
			d = fabs(a[i * n + j]);
			if (d > ff) {
				p = i;
				q = j;
				goto loop;
			}
		}

	if (ff < eps)
		return;

	goto loop0;

	loop: u = p * n + q;
	w = p * n + p;
	t = q * n + p;
	s = q * n + q;
	x = -a[u];
	y = (a[s] - a[w]) / 2.0;
	omega = x / sqrt(x * x + y * y);

	if (y < 0.0)
		omega = -omega;

	sn = 1.0 + sqrt(1.0 - omega * omega);
	sn = omega / sqrt(2.0 * sn);
	cn = sqrt(1.0 - sn * sn);
	fm = a[w];
	a[w] = fm * cn * cn + a[s] * sn * sn + a[u] * omega;
	a[s] = fm * sn * sn + a[s] * cn * cn - a[u] * omega;
	a[u] = 0.0;
	a[t] = 0.0;

	for (j = 0; j <= n - 1; j++)
		if ((j != p) && (j != q)) {
			u = p * n + j;
			w = q * n + j;
			fm = a[u];
			a[u] = fm * cn + a[w] * sn;
			a[w] = -fm * sn + a[w] * cn;
		}

	for (i = 0; i <= n - 1; i++)
		if ((i != p) && (i != q)) {
			u = i * n + p;
			w = i * n + q;
			fm = a[u];
			a[u] = fm * cn + a[w] * sn;
			a[w] = -fm * sn + a[w] * cn;
		}

	for (i = 0; i <= n - 1; i++) {
		u = i * n + p;
		w = i * n + q;
		fm = v[u];
		v[u] = fm * cn + v[w] * sn;
		v[w] = -fm * sn + v[w] * cn;
	}

	goto loop1;
}

} /* End Of namespace ISLIB */


