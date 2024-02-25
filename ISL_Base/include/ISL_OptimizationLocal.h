/**
*	@file	ISL_OptimizationLocal.h
*	@brief	[Header file of ISL_OptimizationLocal], 局部最优化；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-07-1
*	@refer
*/


#ifndef PAI_FRAME_ISEISLIB_OptimizationLocal_H_
#define PAI_FRAME_ISEISLIB_OptimizationLocal_H_

#include "ISL_UserDefine.h"

namespace ISLib {

/** @class ISL_dsnse
    @brief 求非线性方程组的梯度法.

	用最速下降法求非线性方程组的一组实数解 \n
*/
class ISL_EXPORT ISL_dsnse
{
public:
	/** @brief 构造函数 */
	ISL_dsnse(){}

	/** @brief 析构函数（虚函数） */
	virtual ~ISL_dsnse(){}

	/** @brief 计算方程组左端函数值方（虚函数），本类需要继承使用，并重新实现dsnsef函数,
	 *  @param[in] x[n] 输入x序列
	 *  @param[in] y[n] 输入y序列
	 *  @param[in] n 	输入y序列的长度
	 *  @return 返回方程的平方和
	 *  @note 无
	 *  @see
	 */
	virtual double dsnsef(double *x, double *y, int n){
		x = NULL;
		y = NULL;
		n = 0;
		return 0;
	}

	/** @brief 计算函数主体
	 *  @param[in] x[n] 	输入x序列
	 *  @param[in] n 		输入y序列的长度
	 *  @param[in] eps		控制精度
	 *  @param[in] js		迭代次数
	 *  @param[out] x[n]	返回一组实根
	 *  @return
	 *  @note 无
	 *  @see 测试程序（optimizationLocal_test）
	 */
	int run(double *x, int n, double eps, int js);
};




/** @class ISL_dnetn
    @brief 求非线性方程组的拟牛顿法.

	用最速下降法求非线性方程组的一组实数解 \n
*/
class ISL_EXPORT ISL_dnetn {
private:
	int agaus(double *a, double *b, int n);

public:
	/** @brief 构造函数 */
	ISL_dnetn(){}

	/** @brief 析构函数（虚函数） */
	virtual ~ISL_dnetn(){}


	/** @brief 计算方程组左端函数值（虚函数），本类需要继承使用，并重新实现dnetnf函数
	 *  @param[in] x[n] 输入x序列
	 *  @param[in] y[n] 输入y序列
	 *  @param[in] n 	输入y序列的长度
	 *  @return 返回方程的平方和
	 *  @note 无
	 *  @see 测试程序（optimizationLocal_test）
	 */
	virtual void dnetnf(double *x, double *y, int n){
		x = NULL;
		y = NULL;
		n = 0;
	}



	/** @brief 计算函数主体
	 *  @param[in] x[n] 	输入x序列
	 *  @param[in] n 		输入y序列的长度
	 *  @param[in] eps		控制精度
	 *  @param[in] t		控制h的大小，0<t<1
	 *  @param[in] h		增量初值
	 *  @param[in] k		迭代次数
	 *  @param[out] x[n]	返回方程组的一组实数解
	 *  @return
	 *  @note 无
	 *  @see 测试程序（optimizationLocal_test）
	 */
	int run(int n, double eps, double t, double h, double *x, int k);
};

}/*End of ISLib*/

#endif /* ISL_OptimizationLocal_H_ */
