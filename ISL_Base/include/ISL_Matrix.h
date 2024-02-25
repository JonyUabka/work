/**
*	@file	ISL_Matrix.h
*	@brief	[Header file of Matrix Functions], 矩阵运算；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-06-03
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_MATRIX_H_
#define PAI_FRAME_ISEISLIB_MATRIX_H_

#include "ISL_UserDefine.h"

namespace ISLib {

/**
 * @brief	实矩阵相乘
 *
 * @param[in]	a[m][n]			存放矩阵A的元素
 * @param[in]	b[n][k]			存放矩阵B的元素
 * @param[in]	m				矩阵A与乘积矩阵C的行数
 * @param[in]	n				矩阵A的列数，矩阵B的行数
 * @param[in]	k				矩阵B与乘积矩阵C的行数
 * @param[out]	c[m][k]			返回乘积矩阵C=AB的元素，进入函数前请分配好内存空间
 *
 * @return	NO
 */
ISL_EXPORT void ISL_trmul(double *a, double *b, int m, int n, int k, double *c);


/**
 * @brief	复矩阵相乘
 *
 * @param[in]	ar[m][n]		存放矩阵A的实部元素
 * @param[in]	ai[m][n]		存放矩阵A的虚部元素
 * @param[in]	br[n][k]		存放矩阵B的实部元素
 * @param[in]	bi[n][k]		存放矩阵B的虚部元素
 * @param[in]	m				矩阵A与乘积矩阵C的行数
 * @param[in]   n				矩阵A的列数，矩阵B的行数
 * @param[in]   k				矩阵B与乘积矩阵C的行数
 * @param[out]  cr[m][k]        返回乘积矩阵C=AB的实部元素， 进入函数前请分配好内存空间
 * @param[out]	ci[m][k]		返回乘积矩阵C=AB的虚部元素， 进入函数前请分配好内存空间
 *
 * @return	NO
 */
ISL_EXPORT void ISL_tcmul(double *ar, double *ai, double *br, double *bi, int m, int n, int k,
		double *cr, double *ci);



/**
 * @brief	一般实矩阵求逆
 *
 * @param[in]	a[n][n]			存放矩阵A。返回时存放其逆矩阵A^-1
 * @param[in]	n				矩阵阶数
 *
 * @return	函数返回整型标志位，0表示A奇异；否则表示正常返回
 */
ISL_EXPORT int ISL_rinv(double *a, int n);



/**
 * @brief	一般复矩阵求逆
 *
 * @param[in]	ar[n][n]		存放矩阵A的实部。返回时存放其逆矩阵A^-1的实部
 * @param[in]	ai[n][n]		存放矩阵A的虚部。返回时存放其逆矩阵A^-1的虚部
 * @param[in]	n				矩阵阶数
 *
 * @return	函数返回整型标志位，0表示A奇异；否则表示正常返回
 */
ISL_EXPORT int ISL_cinv(double *ar, double *ai, int n);



/**
 * @brief	对称正定矩阵的求逆
 *
 * @param[in]	a[n][n]			存放 对称正定矩阵A。返回时存放其逆矩阵A^-1
 * @param[in]	n				矩阵阶数
 *
 * @return	函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；	若返回标志值大于0；否则表示正常返回。
 */
ISL_EXPORT int ISL_ssgj(double *a, int n);



/**
 * @brief	托伯利兹矩阵求逆的特兰持方法
 *
 * @param[in]	t[n]			存放 T型矩阵中的元素。t0, t1, t2,..., t(n-1)
 * @param[in]	tt[n]			后n-1个元素存放 T型矩阵中的元素
 * @param[in]	n				T型矩阵阶数
 * @param[out]	b[n][n]			返回T型矩阵的逆矩阵
 *
 * @return	函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；	若返回标志值大于0；否则表示正常返回。
 */
ISL_EXPORT int ISL_trch(double *t, double *tt, int n, double *b);



/**
 * @brief	求一般行列式的值，用全选主元高斯（Gauss）消去法计算n阶方阵A所对应的行列式值,
 * 			用全选主元高斯（Gauss）消去法对方阵A进行一系列变换使之成为上三角矩阵，其对角线上的个元素乘积即为行列式值。
 *
 * @param[in]	a[n][n]			存放方阵A的元素，返回时被破坏
 * @param[in]	n				方阵的阶数
 *
 * @return	函数返回行列式值。
 */
ISL_EXPORT double ISL_sdet(double *a, int n);



/**
 * @brief	求矩阵的秩，用全选主元高斯（Gauss）消去法计算矩阵的秩
 *
 * @param[in]	a[m][n]			存放mxn阶矩阵A的元素，返回时被破坏
 * @param[in]	m				矩阵A的行数
 * @param[in]	n				矩阵A的列数
 *
 * @return	函数返回行列式值。
 */
ISL_EXPORT int ISL_rank(double *a, int m, int n);



/**
 * @brief	对称正定矩阵的乔里斯基分解与行列式求值,用乔里斯基（Cholesky）分解法求对称正定矩阵的三角分解，并求行列式的值
 *
 * @param[in]	a[m][n]			存放对称正定矩阵A，返回时其下三角部分存放分解得到的下三角阵L，其余元素均为0
 * @param[in]	n				矩阵A的行数
 * @param[out]	det				指向行列式的值
 *
 * @return	函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；	若返回标志值大于0；否则表示正常返回。
 */
ISL_EXPORT int ISL_chol(double *a, int n, double &det);



/**
 * @brief	矩阵的三角分解
 *
 * @param[in]	a[m][n]			存放n阶矩阵A，返回时存放Q矩阵
 * @param[in]	n				矩阵阶数
 * @param[out]	l[n][n]			返回时存放下三角矩阵L
 * @param[out]	u[n][n]			返回时存放上三角矩阵U
 *
 * @return	函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；若返回标志值大于0；否则表示正常返回。
 */
ISL_EXPORT int ISL_lluu(double *a, int n, double *l, double *u);


/**
 * @brief	一般实矩阵的QR分解，用豪斯荷尔德（Householder）变换对一般mxn阶的实矩阵进行QR分解
 *
 * @param[in]	a[m][n]			存放mxn的实矩阵A。返回时其右上三角部分存放QR分解中的上三角矩阵R
 * @param[in]	m				实矩阵A的行数
 * @param[in]	n				实矩阵A的列数
 * @param[out]	q[n][n]			返回QR分解中的正交矩阵Q
 *
 * @return	函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；若返回标志值大于0；否则表示正常返回。
 */
ISL_EXPORT int ISL_maqr(double *a, int m, int n, double *q);


/**
 * @brief	一般实矩阵的奇异值分解，用豪斯荷尔德（Householder）变换以及变形QR对一般实矩阵A进行奇异值分解
 *
 * @param[in]	a[m][n]		存放mxn的实矩阵A。
 * @param[out]	a[m][n]		返回时其对角线给出奇异值（以非递增次序排列）
 * @param[in]	m			实矩阵A的行数
 * @param[in]	n			实矩阵A的列数
 * @param[out]	u[m][m]		返回左奇异向量U
 * @param[out]	v[n][n]		返回右奇异向量VT
 * @param[in]	eps			给定的精度要求
 * @param[in]	ka			其值为max(m, n)+1
 *
 * @return	函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；若返回标志值大于0；否则表示正常返回。
 */
ISL_EXPORT int ISL_muav(double *a, int m, int n, double *u, double *v, double eps, int ka);


/**
 * @brief	求广义逆的奇异值分解法，利用奇异值分解求一般mxn阶实矩阵A的广义逆A+
 *
 * @param[in]	a[m][n]		存放mxn的实矩阵A。
 * @param[out]	a[m][n]		返回时其对角线给出奇异值（以非递增次序排列），其余元素均为0。
 * @param[in]	m			实矩阵A的行数
 * @param[in]	n			实矩阵A的列数
 * @param[out]	aa[m][m]	返回A的广义逆A+
 * @param[in]	eps			给定的精度要求
 * @param[out]	u[m][m]		返回右奇异向量U
 * @param[out]	v[n][n]		返回右奇异向量VT
 * @param[out]	ka			其值为max(m, n)+1
 *
 * @return	函数返回整型标志位，若返回标志值小于0，则表示程序工作失败；若返回标志值大于0；否则表示正常返回。
 */
ISL_EXPORT int ISL_ginv(double *a, int m, int n, double *aa, double eps, double *u, double *v, int ka);

// ============================
// ===	矩阵特征值运算与特征向量的计算	===
// ============================

/**
 * @brief	约化对称矩阵为对称三角阵的豪斯赫尔德变换法
 *
 * @param[in]	a[n][n]		存放n阶实对称矩阵A
 * @param[in]	n			实对称矩阵A的阶数
 * @param[out]	q[n][n]		返回豪斯荷尔德变换的乘积矩阵Q。在与nj_sstq()联用时，若将Q矩阵作为函数nj_sstq()
 *							中的一个参数，则可以计算一般实对称矩阵的全部特征值及相应的特征向量
 * @param[in]	b[n]		返回对称三角阵中的主对角线元素
 * @param[out]	c[n]		前n-1个元素返回对称三角阵中的次对角线元素
 *
 * @return	无
 */
ISL_EXPORT void ISL_strq(double *a, int n, double *q, double *b, double *c);


/**
 * @brief	求对称三对角阵的全部特征值与特征向量
 *
 * @param[in]	n			实对称三角阵的阶数
 * @param[in]	b[n]		存放n阶实对称三角阵的主对角线上的元素。返回时存放全部特征值
 * @param[out]	b[n]		返回时存放全部特征值
 * @param[in]	c[n]		前n-1个元素存放对称三角阵中的次对角线上的元素
 * @param[out]	q[n][n]		若存放n阶单位矩阵，则返回实对称三对角阵T的特征向量组；若存放由函数nj_strq（）所返回的一般
 *							实对称矩阵A的豪斯荷尔德的乘积矩阵Q，则返回实对称矩阵A的特征向量组。其中q中的第j列为数组b中第j个
 *							特征值对应的特征向量
 * @param[in]	eps			控制精度要求
 * @param[out]	l			允许的最大迭代次数
 *
 * @return	函数返回标志值。若返回的标志值小于0，则表示程序工作失败；若返回的标志值大于0，则说明程序正常返回
 */
ISL_EXPORT int ISL_sstq(int n, double *b, double *c, double *q, double eps, int l);


/**
 * @brief	约化一般实矩阵为赫申伯格矩阵的初等相似变换法
 *
 * @param[in]	n			实对称矩阵A的阶数
 * @param[in]	a[n][n]		存放一般实矩阵A
 * @param[out]	a[n][n]		返回上H矩阵
 *
 * @return	无
 */
ISL_EXPORT void ISL_hhbg(double *a, int n);


/**
 * @brief	求赫申伯格矩阵全部特征值的QR方法
 *
 * @param[in]	a[n][n]		存在上H矩阵A
 * @param[in]	n			上H矩阵A的阶数
 * @param[out]	u[n]		返回n个特征值的实部
 * @param[out]	v[n]		返回n个特征值的虚部
 * @param[in]	eps			控制精度要求一般
 * @param[in]	jt			允许的最大迭代次数
 *
 * @return	函数返回标志值。若返回的标志值小于0，则表示程序工作失败；若返回的标志值大于0，则说明程序正常返回
 */
ISL_EXPORT int ISL_hhqr(double *a, int n, double *u, double *v, double eps, int jt);


/**
 * @brief	求实对称矩阵特征值与特征向量的雅可比法
 *
 * @param[in]	a[n][n]		存放n阶实对称矩阵。返回时对角线上存放n个特征值
 * @param[in]	n			实对称矩阵A的阶数
 * @param[out]	v[n][n]		返回特征向量。其中第j列为与第j个特征值对应的特征向量
 * @param[in]	eps			控制精度要求
 * @param[in]	jt			控制最大迭代次数
 *
 * @return	函数返回标志值。若返回的标志值小于0，则表示程序工作失败；若返回的标志值大于0，则说明程序正常返回
 */
ISL_EXPORT int ISL_jcbi(double *a, int n,double * v, double eps, int jt);


/**
 * @brief	求实对称矩阵特征值与特征向量的雅可比过关法
 *
 * @param[in]	a[n][n]		存放n阶实对称矩阵。返回时对角线上存放n个特征值
 * @param[in]	n			实对称矩阵A的阶数
 * @param[out]	v[n][n]		返回特征向量。其中第j列为与第j个特征值对应的特征向量
 * @param[in]	eps			控制精度要求
 *
 * @return	函数返回标志值。若返回的标志值小于0，则表示程序工作失败；若返回的标志值大于0，则说明程序正常返回
 */
ISL_EXPORT void ISL_jcbj(double *a, int n, double *v, double eps);

} /* End Of namespace ISLIB */
#endif /* PAI_FRAME_ISEISLIB_MATRIX_H_ */
