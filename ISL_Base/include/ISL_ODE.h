/**
*	@file	ISL_ODE.h
*	@brief	[Header file of ODE Functions], 常微分方程数值解；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-08-18
*	@refer
*/

#ifndef PAI_FRAME_ISEISLIB_ODE_H_
#define PAI_FRAME_ISEISLIB_ODE_H_

#include "ISL_UserDefine.h"

namespace ISLib {


/** @class ISL_gelr1
    @brief 用改进的Euler法对一阶微分方程组进行定步长全区间积分
*/
class ISL_EXPORT ISL_gelr1
{
public:
	/** @brief 构造函数 */
	ISL_gelr1(){}

	/** @brief 析构函数（虚函数） */
	virtual ~ISL_gelr1(){}

	/** @brief 计算方程组左端函数值方（虚函数），本类需要继承使用，并重新实现gelr1函数,
	 *  @param[in] d[n] 输入d序列
	 *  @param[in] y[n] 输入y序列
	 *  @param[in] n
	 *  @param[in] t
	 *  @param[out] d[n] 输入d序列
	 *  @return 无
	 *  @note 无
	 *  @see
	 */
	virtual void gelr1f(double t, double *y, int n, double *&d)
	{
		y = NULL;
		d = NULL;
		n = 0;
		t = 0;
	}

	/** @brief 计算函数主体
	 *  @param[in] t 		双精度实型变量，积分的起始点
	 *  @param[in] y[n] 	双精度实型数组，存放起始点的函数值
	 *  @param[in] n		整型变量，微分方程组个数，也是未知数的个数
	 *  @param[in] h		双精度实型变量，积分步长
	 *  @param[in] k		整型变量，积分步数
	 *  @param[out] z[n][k]		双精度实型数组，存放k个积分点上的函数值
	 *  @return 无
	 *  @note 无
	 *  @see 测试程序（ODE_test）
	 */
	void run(double t, double *y, int n, double h, int k, double *z);
};



/** @class ISL_gmrsn
    @brief 用变步长Merson法对一阶微分方程组进行全区间积分
*/
class ISL_EXPORT ISL_gmrsn
{
public:
	/** @brief 构造函数 */
	ISL_gmrsn(){}

	/** @brief 析构函数（虚函数） */
	virtual ~ISL_gmrsn(){}

	/** @brief 计算方程组左端函数值方（虚函数），本类需要继承使用，并重新实现ggil2f函数,
	 *  @param[in] d[n] 输入d序列
	 *  @param[in] y[n] 输入y序列
	 *  @param[in] n
	 *  @param[in] t
	 *  @param[out] d[n] 输入d序列
	 *  @return 无
	 *  @note 无
	 *  @see
	 */
	virtual void gmrsnf(double t, double *y, int n, double *&d)
	{
		y = NULL;
		d = NULL;
		n = 0;
		t = 0;
	}

	/** @brief 计算函数主体
	 *  @param[in] t 		双精度实型变量，积分的起始点
	 *  @param[in] h		双精度实型变量，积分步长
	 *  @param[in] y[n] 	双精度实型数组，存放起始点的函数值
	 *  @param[in] n		整型变量，微分方程组个数，也是未知数的个数
	 *  @param[in] eps		双精度实型变量，控制积分一步的精度
	 *  @param[in] k		整型变量，积分步数
	 *  @param[out] z[n][k]		双精度实型数组，存放k个积分点上的函数值
	 *  @return 无
	 *  @note 无
	 *  @see 测试程序（ODE_test）
	 */
	void run(double t, double h, int n, double *y, double eps, int k, double *z);
};





/** @class ISL_grkt1
    @brief 定步长四阶Runge-Kutta法对一阶微分方程组进行全区间积分
*/
class ISL_EXPORT ISL_grkt1
{
public:
	/** @brief 构造函数 */
	ISL_grkt1(){}

	/** @brief 析构函数（虚函数） */
	virtual ~ISL_grkt1(){}

	/** @brief 计算方程组左端函数值方（虚函数），本类需要继承使用，并重新实现ggil2f函数,
	 *  @param[in] d[n] 输入d序列
	 *  @param[in] y[n] 输入y序列
	 *  @param[in] n
	 *  @param[in] t
	 *  @param[out] d[n] 输入d序列
	 *  @return 无
	 *  @note 无
	 *  @see
	 */
	virtual void grkt1f(double t, double *y, int n, double *&d)
	{
		y = NULL;
		d = NULL;
		n = 0;
		t = 0;
	}

	/** @brief 计算函数主体
	 *  @param[in] t 		双精度实型变量，积分的起始点
	 *  @param[in] h		双精度实型变量，积分步长
	 *  @param[in] y[n] 	双精度实型数组，存放起始点的函数值
	 *  @param[in] n		整型变量，微分方程组个数，也是未知数的个数
	 *  @param[in] k		整型变量，积分步数
	 *  @param[out] z[n][k]		双精度实型数组，存放k个积分点上的函数值
	 *  @return 无
	 *  @note 无
	 *  @see 测试程序（ODE_test）
	 */
	void run(double t, double *y, int n, double h, int k, double *z);
};


}/* End Of namespace ISLIB */

#endif /* PAI_FRAME_ISEISLIB_ODE_H_ */
