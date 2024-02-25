/**
*	@file	ISL_NUMERICAL.h 
*	@brief	[Header file of NUMERICAL], NUMERICAL；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang, He Kai], 刘百红、杨强、宋志翔、何恺；
*	@date	2014-05-22
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_NUMERICAL_H
#define PAI_FRAME_ISEISLIB_NUMERICAL_H

#include "ISL_UserDefine.h"

namespace ISLib {


/**
* @brief	一维多项式求值
*
* @param[in]	a[n]			存放n-1次多项式的n个系数
* @param[in]	n				多项式的项数
* @param[in]	x				指定的自变量值
*
* @return	返回多项式值
*/
ISL_EXPORT double ISL_plyv(double * a, int n, double x);



/**
* @brief	一维多项式多组求值
*
* @param[in]	a[n]			存放n-1次多项式的n个系数,返回时其值将改变
* @param[in]	n				多项式的项数
* @param[in]	x[m]			存放给定的m个自变量值
* @param[in]	m				给定自变量的个数
* @param[out]	p[m]			数组大小为m， 返回时存放与给定m个自变量值对应的多项式值, 进入函数前请分配好内存空间
*
* @return	no
*/
ISL_EXPORT void ISL_plys(double *a, int n, double *x, int m, double *p);



/**
* @brief	二维多项式求值
*
* @param[in]	a[m][n]			存放二维次多项式的系数
* @param[in]	m				自变量x的最高次数为m-1
* @param[in]	n				自变量y的最高次数为n-1
* @param[in]	x				给定的x的自变量值
* @param[in]	y				给定的y的自变量值
*
* @return	返回多项式值
*/
ISL_EXPORT double ISL_bply(double *a, int m, int n, double x, double y);



/**
* @brief	复系数多项式求值
*
* @param[in]		ar[n], ai[n]	分别存放多项式系数的实部与虚部
* @param[in]		n				多项式的项数，其最高次数为n-1
* @param[in]		x, y			给定复数Z的实部与虚部
* @param[In/Out]	&u, &v			指向返回多项式值p（z）的实部与虚部
*
* @return	no
*/
ISL_EXPORT void ISL_cply(double *ar, double *ai, int n, double x, double y, double &u, double &v);
ISL_EXPORT void ISL_cmul(double a, double b, double c, double d, double &e, double &f);



/**
* @brief	多项式相乘
*
* @param[in]		p[m]			存放多项式p(x)的系数
* @param[in]		m				多项式p(x)的项数，其最高次数为m-1
* @param[in]		q[n]			存放多项式q(x)的系数
* @param[in]		n				多项式q(x)的项数，其最高次数为n-1
* @param[In/Out]	s[k]			返回乘积多项式的系数, 进入函数前请分配好内存空间
* @param[in]		k				乘积多项式s(x)的项数，其最高次数为k-1，其中k=m+n-1
*
* @return	no
*/
ISL_EXPORT void ISL_pmul(double *p, int m, double *q, int n, double *s, int k);



/**
* @brief	复系数多项式相乘
*
* @param[in]		pr[m],pi[m]		存放多项式p(x)的系数的实部与虚部
* @param[in]		m				多项式p(x)的项数，其最高次数为m-1
* @param[in]		qr[n],qi[n]		存放多项式q(x)的系数的实部与虚部
* @param[in]		n				多项式q(x)的项数，其最高次数为n-1
* @param[In/Out]	sr[k],si[k]		返回乘积多项式的系数的实部与虚部, 进入函数前请分配好内存空间
* @param[in]		k				乘积多项式s(x)的项数，其最高次数为k-1，其中k=m+n-1
*
* @return	no
*/
ISL_EXPORT void ISL_cpml(double *pr, double * pi, int m, double *qr, double *qi, int n, double *sr, double *si, int k);



/**
* @brief	多项式相除
*
* @param[in]		p[m]		存放多项式p(x)的系数。返回时其中的值将被破坏
* @param[in]		m			多项式p(x)的项数，其最高次数为m-1
* @param[in]		q[n]		存放多项式q(x)的系数
* @param[in]		n			多项式q(x)的项数，其最高次数为n-1
* @param[In/Out]	s[k]		返回商多项式s（x）的系数, 进入函数前请分配好内存空间
* @param[in]		k			商多项式s(x)的项数，其最高次数为k-1，其中k=m+n-1
* @param[In/Out]	r[l]		返回商多项式r（x）的系数, 进入函数前请分配好内存空间
* @param[in]		l			余多项式r(x)的项数，其最高次数为l-1，其中l=n-1
*
* @return	no
*/
ISL_EXPORT void ISL_pdiv(double *p, int m, double *q, int n, double *s, int k, double *r,
				int l);



/**
* @brief	复系数多项式相除
*
* @param[in]		pr[m],pi[m]		存放多项式p(x)的系数的实部与虚部。返回时其中的值将被破坏
* @param[in]		m				多项式p(x)的项数，其最高次数为m-1
* @param[in]		qr[n],qi[n]		存放多项式q(x)的系数的实部与虚部
* @param[in]		n				多项式q(x)的项数，其最高次数为n-1
* @param[In/Out]	sr[k],si[k]		返回商多项式s（x）的系数的实部与虚部, 进入函数前请分配好内存空间
* @param[in]		k				商多项式s(x)的项数，其最高次数为k-1，其中k=m+n-1
* @param[In/Out]	rr[l],ri[l]		返回商多项式r（x）的系数的实部与虚部, 进入函数前请分配好内存空间
* @param[in]		l				余多项式r(x)的项数，其最高次数为l-1，其中l=n-1
*
* @return	no
*/
ISL_EXPORT void ISL_cpdv(double *pr, double *pi, int m, double *qr, double *qi, int n,
		double *sr, double *si, int k, double *rr, double *ri, int l);
ISL_EXPORT void ISL_cdiv(double a, double b, double c, double d, double &e, double &f);



/**
* @brief	求解实系数方程组的全选主元高斯消去法
*
* @param[in]		a[n][n]			存放方程组的系数矩阵，返回时将被破坏
* @param[In/Out]	b[n]			存放方程组右端的常数向量，返回方程组的解向量
* @param[in]		n				方程组的阶数
*
* @return	函数返回整型标志位，若返回标志值为0，则表示程序工作失败（因系数矩阵奇异）；
*  				若返回标志值不为0，则表示正常返回。
*/
ISL_EXPORT int ISL_gaus(double *a, double *b, int n);



/**
* @brief	求解实系数方程组的全选主元高斯-约当消去法
*
* @param[in]		a[n][n]			存放方程组的系数矩阵，返回时将被破坏
* @param[In/Out]	b[n][m]			存放方程组右端的m组常数向量，返回方程组的m组解向量
* @param[in]		n				方程组的阶数
* @param[in]		m				方程组右端常数向量的组数
*
* @return	函数返回整型标志位，若返回标志值为0，则表示程序工作失败（因系数矩阵奇异）；
*  				若返回标志值不为0，则表示正常返回。
*/
ISL_EXPORT int ISL_gjdn(double *a, double *b, int n, int m);



/**
* @brief	求解复系数方程组的全选主元高斯消去法
*
* @param[in]		ar[n][n],ai[n][n]			存放方程组的复系数矩阵的实部和虚部，返回时将被破坏
* @param[In/Out]	br[n],bi[n]					存放方程组右端的复常数向量的实部与虚部，返回方程组的解向量的实部与虚部
* @param[in]		n							方程组的阶数
*
* @return	函数返回整型标志位，若返回标志值为0，则表示程序工作失败（因系数矩阵奇异）；
*  				若返回标志值不为0，则表示正常返回。
*/
ISL_EXPORT int ISL_cgas(double *ar, double *ai, int n, double *br, double *bi);



/**
* @brief	求解复系数方程组的全选主元高斯-约当消去法
*
* @param[in]		ar[n][n],ai[n][n]			存放方程组的复系数矩阵的实部和虚部，返回时将被破坏
* @param[In/Out]	br[n][m],bi[n][m]			存放方程组右端的m组常数向量的实部和虚部，返回方程组的m组解向量的实部和虚部
* @param[in]		n							方程组的阶数
* @param[in]		m							方程组右端常数向量的组数
*
* @return	函数返回整型标志位，若返回标志值为0，则表示程序工作失败（因系数矩阵奇异）；
*  				若返回标志值不为0，则表示正常返回。
*/
ISL_EXPORT int ISL_cjdn(double *ar, double *ai, double *br, double *bi, int n, int m);



/**
* @brief	求解三对角线方程组的追赶法
*
* @param[in]		b[m]		以行为主，存放三对角矩阵中三条对角线上的元素
* @param[in]		n			方程组的阶数
* @param[in]		m			三对角矩阵三条对角线上的元素个数，其值应为m=3n-2
* @param[In/Out]	d[n]		存放方程组右端的常数向量，返回方程组的解向量
*
* @return	函数返回整型标志位，若返回标志值小于<0，则表示m的值不正确；
*  				若返回标志值=0，则表示程序工作失败；
*  				若返回标志值>0，则表示正常返回。
*/
ISL_EXPORT int ISL_trde(double *b, int n, int m, double *d);



/**
* @brief	求解一般带型方程组
*
* @param[In/Out]	b[n][il]		存放带型矩阵A中带区内的元素。返回时将被破坏
* @param[In/Out]	d[n][m]			存放方程组右端的m组的常数向量。返回方程组的m组解向量
* @param[in]		n				方程组的阶数
* @param[in]		l				系数矩阵的半带宽h
* @param[in]		il				系数矩阵的带宽2h+1。应满足il=2l+1
* @param[in]		m				方程组右端常数向量的组数
*
* @return	函数返回整型标志位，若返回标志值小于<0，则表示m的值不正确；
*  				若返回标志值=0，则表示程序工作失败；
*  				若返回标志值>0，则表示正常返回。
*/
ISL_EXPORT int ISL_band(double *b, double *d, int n, int l, int il, int m);



/**
* @brief	求解对称方程组的分解法
*
* @param[In/Out]	a[n][n]			存放方程组的系数矩阵（应为对称矩阵)。返回时将被破坏
* @param[in]		n				方程组的阶数
* @param[in]		m				方程组右端常数向量的组数
* @param[In/Out]	c[n][m]			存放方程组右端m组常数向量。返回方程组的m组解向量
*
* @return	函数返回标志值。若返回的标志值小于0，则表示程序工作失败，若返回的标志值大于0，则表示正常返回。
*/
ISL_EXPORT int ISL_ldle(double *a, int n, int m, double *c);



/**
* @brief	求解对称正定方程组的平方根法
*
* @param[In/Out]	a[n][n]			存放对称正定的系数矩阵。返回时其上三角部分存放分解后的U矩
* @param[in]		n				方程组的阶数
* @param[in]		m				方程组右端常数向量的组数
* @param[In/Out]	d[n][m]			存放方程组右端m组常数向量。返回方程组的m组解向量
*
* @return	函数返回标志值。若返回的标志值小于0，则表示程序工作失败，若返回的标志值大于0，则表示正常返回。
*/
ISL_EXPORT int ISL_chlk(double *a, int n, int m, double *d);



/**
* @brief	求解大型稀疏方程组
*
* @param[In/Out]	a[n][n]			存放对称正定的系数矩阵。返回时其上三角部分存放分解后的U矩
* @param[in]		n				方程组的阶数
* @param[In/Out]	b[n][m]			存放方程组右端常数向量。返回方程组的解向量
*
* @return	函数返回标志值。若返回的标志值为0，则表示系数矩阵奇异;若返回的标志值不为0，则表示正常返回。
*/
ISL_EXPORT int ISL_ggje(double *a, int n, double *b);



/**
* @brief	求解托伯利兹方程组的列文逊方法
*
* @param[in]	t[n]			存放n阶T型矩阵中的元素
* @param[in]	n				方程组的阶数
* @param[in]	b[n]			存放方程组右端的常数向量
* @param[out]	x[n]			返回方程组的解向量
*
* @return	函数返回标志值。若返回的标志值小于0，则表示程序工作失败，若返回的标志值大于0，则表示正常返回。
*/
ISL_EXPORT int ISL_tlvs(double *t, int n, double *b, double *x);



/**
* @brief	高斯-赛德尔迭代法
*
* @param[in]	a[n][n]			存放方程组的系数矩阵
* @param[in]	b[n]			存放方程组右端的常数向量
* @param[in]	n				方程组的阶数
* @param[in]	eps				给定的精度要求
* @param[out]	x[n]			返回方程组的解向量
*
* @return	函数返回标志值。若返回的标志值小于0，则表示系数矩阵不具有主对角线占绝对优势，若返回的标志值大于0，则表示正常返回。
*/
ISL_EXPORT int ISL_gsdl(double *a, double *b, int n, double *x, double eps);



/**
* @brief	求解对称正定方程组的共轭梯度法
*
* @param[in]	a[n][n]			存放对称正定矩阵A
* @param[in]	b[n]			存放方程组右端的常数向量
* @param[in]	n				方程组的阶数
* @param[in]	eps				给定的精度要求
* @param[out]	x[n]			返回方程组的解向量
*
* @return	no
*/
ISL_EXPORT void ISL_grad(double *a, int n, double *b, double eps, double *x);



/**
* @brief	求解线性最小二乘问题的豪斯荷尔德变换法
*
* @param[in]	a[m][n]			存放超定方程组的系数矩阵A。返回时存放QR分解式中的R矩阵
* @param[in]	m				系数矩阵A的行数，m≥n
* @param[in]	n				系数矩阵A的列数，n≤m
* @param[in]	b[m]			存放方程组右端的常数向量。返回时前n个分量存放方程组的最小二乘解
* @param[out]	q[m][m]			返回时存放QR分解式中的正交矩阵Q
*
* @return	函数返回标志值。若返回的标志值为0，则表示程序工作失败（如A列线性相关）；
*  				若返回的标志值不为0，则表示止常返回。
*/
ISL_EXPORT int ISL_gmqr(double *a, int m, int n, double *b, double *q);



/**
* @brief	求解线性最小二乘问题的广义逆法
*
* @param[in]	a[m][n]			存放超定方程组的系数矩阵A。返回时其对角线依次给出奇异值，其余元素为0
* @param[in]	m				系数矩阵A的行数
* @param[in]	n				系数矩阵A的列数
* @param[in]	b[m]			存放方程组右端的常数向量
* @param[out]	x[n]			返回超定方程组的最小二乘解
* @param[out]	aa[n][m]		返回系数矩阵A的广义逆 A'
* @param[in]	eps				奇异值分解中的控制精度要求
* @param[out]	u[m][m]			返回系数矩阵A的奇异值分解式中的左奇异向量U
* @param[out]	v[n][n]			返回系数矩阵A的奇异值分解式中的右奇异向量Vt
* @param[in]	ka				ka = max(m,n)+1
*
* @return	函数返回标志值。函数返回标志值。若返网的标志值小于0，则表示程序工作失败；
*  				若返回的标志值大于0，则表示正常返回。
*/
ISL_EXPORT int ISL_gmiv(double *a, int m, int n, double *b, double *x, double *aa, double eps,
		double *u, double *v, int ka);



/**
* @brief	求解病态方程组
*
* @param[in]	a[m][n]			存放方程组的系数矩阵
* @param[in]	n				方程组的解数
* @param[in]	b[m]			存放方程组右端的常数向量
* @param[out]	x[n]			返回方程组的解向量
* @param[in]	eps				控制精度要求
*
* @return	函数返回标志值。函数返回标志值。若返回的标志值为O，则表示程序工作失败；
*  			若返回的标志值不为0，则表示正常返回。
*/
ISL_EXPORT int ISL_bint(double *a, int n, double *b, double eps, double *x);


}/*End of ISLib*/

#endif
