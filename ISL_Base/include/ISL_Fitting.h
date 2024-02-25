/**
*	@file	ISL_Fitting.h
*	@brief	[Header file of Fitting Functions], 拟合；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang, Chen Ke], 刘百红、杨强、宋志翔、陈科；
*	@date	2014-06-04
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_FITTING_H_
#define PAI_FRAME_ISEISLIB_FITTING_H_

#include "ISL_UserDefine.h"

namespace ISLib {

/**
 * @brief	最小二乘曲线拟合,用最小二乘法求给定点的拟合多项式，此函数用于模拟地震记录的振幅谱
 *
 * @param[in]	x[n]		存放n个数据点的x坐标；
 * @param[in]	y[n]		存放n个数据点的y坐标；
 * @param[in]	n			点的个数；
 * @param[Out]	cof[nterm]	返回deg-1次拟合多项式的deg个系数；
 * @param[in]	nterm		输入拟合多项式的项数
 * @param[Out]	simu[n]		返回每一个频率上拟合多项式的和
 *
 * @return	无
 */
ISL_EXPORT void ISL_lscf(float *x, float *y, int n, float *cof, int nterm, float *simu);


/**
* @brief	拟合功率谱函数, 利用公式拟合功率谱曲线
*
* @param[in]	inSpec		输入的各个点上的瞬时谱
* @param[in]	inFreq		输入的各个点上的瞬时频率
* @param[in]	num			瞬时谱的点数
* @param[in]	deltaFre	频率域采样间隔（秒）
* @param[out]	outSpec		输出拟合后的瞬时谱
*
* @return	无
*/
ISL_EXPORT void ISL_spectrumCurveFitting(float *inFreq, float *inSpec, int num, float deltaFre, float *outSpec);



/**
 * @brief	随机样本分析
 *
 * @param[in]	x[n]		存放随机变量的n个样本点值；
 * @param[in]	n			随机样本点数；
 * @param[in]	x0			直方图中随机变量的起始值；
 * @param[in]	h			直方图中随机变量等区间长度值；
 * @param[in]	m			直方图中区间总数
 * @param[in]	l			标志。若l=0，表示不需输出直方图；若l=1，表示需要输出直方图
 * @param[out]	dt[3]		dt(0)返回随机样本的算术平均值；dt（1）返回随机样本的方差；dt（2）返回随机样本的标准差； *
 * @param[out]	g[m]		返回m个区间的按高斯分布所应有的近似理论样本点数
 * @param[out]	q[m]		返回落在m个区间中每一个区间上的随机样本实际点数
 *
 * @return	无
 */
ISL_EXPORT void ISL_rhis(double *x, int n, double x0, double h, int m, int l, double *dt, int *g, int *q);



/**
 * @brief	一元线性回归分析
 *
 * @param[in]	x[n]		存放自变量x的n个取值；
 * @param[in]	y[n]		存放与自变量x的n个取值相对应的随机变量y的观测值；
 * @param[in]	n			观测点数；
 * @param[out]	a[2]		a(0)返回回归系数b，a(0)返回回归系数a；
 * @param[out]	dt[6]		dt(0)返回偏差平方和q；dt(1)返回平均标准偏差；dt(2)返回回归平方和p；
 *							dt(3)返回最大偏差umax；dt(4)返回最小偏差umin；dt(5)返回偏差平均值u；
 *
 * @return	无
 */
ISL_EXPORT void ISL_sqt1(double *x, double *y, int n, double *a, double *dt);



/**
 * @brief	多元线性回归分析
 *
 * @param[in]	x[m][n]		每列存放m个自变量的观测值；
 * @param[in]	y[n]		存放随机变量y的n个观测值；
 * @param[in]	m			自变量个数；
 * @param[in]	n			观测数据的组数；
 * @param[out]	a[m+1]		a返回回归系数a0，a1，…，am
 * @param[out]	dt[4]		dt(0)返回偏差平方和q；dt(1)返回平均标准偏差；dt(2)返回复相关系数r；dt(3)返回回归平方和u；
 * @param[out]	v[m]		返回m个自变量的偏相关系数
 *
 * @return	无
 */
ISL_EXPORT void ISL_sqt2(double *x, double *y, int m, int n, double *a, double *dt, double *v);



/**
 * @brief	逐步回归分析
 *
 * @param[in]	n			自变量x的个数；
 * @param[in]	k			观测数据的点数；
 * @param[in]	x[k][n+1]	其中前n列存放自变量因子Xi（i=0,1,…,n-1）的k次观测值，最后一列存放因变量y的k次观测值；
 * @param[in]	f1			欲选入因子时显著性检验的F分布值；
 * @param[in]	f2			欲剔除因子时显著性检验的F分布值；
 * @param[in]	eps			防止系数相关矩阵退化的判据；
 * @param[out]	xx[n+1]		前n个分量返回n个自变量因子的算术平均值Xi（i=0,1,…,n-1），最后一个分量返回闲变量y的算术平均值y~；
 * @param[out]	b[n+1]		返回回归方程中各因子的回归系数及常数项b0,b1,...,bn；
 * @param[out]	v[n+1]		前n个分量返回各因子的偏回归平方和Vi( i=0,1,...,n-1),最后一个分量返回残差平方和q；
 * @param[out]	s[n+1]		前n个分量返回各因子回归系数的标准偏差Si( i=0,1,...,n-1),最后一个分量返回估计的标准偏差s；
 * @param[out]	dt[2]		dt(0)返回复相关系数，dt(1）返回F-检验值；
 * @param[out]	ye[k]		返回对应于K个观测值的因变量条件期望值的k个估计值ei( i=0,1,...,k-1)；
 * @param[out]	yr[k]		返回因变量的K个观测值的残差ei( i=0,1,...,k-1)；
 * @param[out]	r[n+1][n+1]		返回最终的规格化的系数相关矩阵R；
 *
 * @return	无
 */
ISL_EXPORT void ISL_sqt3(int n, int k, double *x, double f1, double f2, double eps,
				double *xx, double *b, double *v, double *s, double *dt,
				double *ye, double *yr, double *r);

} /* End Of namespace ISLIB */
#endif /* PAI_FRAME_ISEISLIB_FITTING_H_ */
