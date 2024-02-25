/****************************************************************************
**
**				ISeis Lib [Base Algorithm] - 插值与拟合
**
**				Writer By	Mouri Song
**							Baihong Liu
**							Qiang Yang
**
**				(Soft Center IGP)
**
**				DATA : 2014-03-19
**
****************************************************************************/

/**
*	@file	ISL_Interpolation.h
*	@brief	[Header file of Interpolation Functions], 插值与拟合；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-03-20
*	@refer	SU CWP
*/

#ifndef PAI_FRAME_ISEISLIB_INTERPOLATION_H
#define PAI_FRAME_ISEISLIB_INTERPOLATION_H

#include "ISL_FFT.h"
#include "ISL_SU_Matrix.h"

namespace ISLib {

/* number of tabulated interpolation coefficients */
#define ICMAX 99 /* must be odd, so that ICMAC-ic!=ic, for ic=0 to ICMAX/2! */
#define NTABLE (ICMAX+1)

#define LTABLE 8
#define NTABLE2 513

#define O2 0.5000000
#define O6 0.1666667

class ISL_EXPORT Interpolation
{
public:
	Interpolation(){ ; }
	~Interpolation(){ ; }

	/**
	* @brief	一维插值( 对序列y(x[0]), y(x[1])进行插值 )
	* @param[in]	nin			输入数据的长度
	* @param[in]	xin			输入数据的X
	* @param[in]	yin			输入数据的Y
	* @param[in]	yinl		y(x)的左极值（设置数值计算范围）
	* @param[in]	yinr		y(x)的右极值（设置数值计算范围）
	* @param[in]	nout		输出数据长度
	* @param[in]	xout		输出的X数据（等间隔数据可由函数ISL_intlin_getx来求得）
	* @param[in]	flag		插值方法
	* 								1 = 线性插值；
	* 								2 = 一元三点插值；
	* 								3 = 拉格朗日（Lagrange）插值；
	* 								4 = 连分式插值;
	* 								5 = 埃特金插值;
	*
	*
	* @param[in]	eps			埃特金插值方法中的调整系数(如不选埃特金插值方法，则此参数不参与计算)
	*
	* @param[out]	yout		计算结果数据的数组【浮点型】
	*
	* @return	无
	*/
	static void ISL_intlin (int nin, float xin[], float yin[],
							float yinl, float yinr,
							int nout, float xout[], float yout[],
							int flag = 1, float eps = 0.5);


	/**
	* @brief	正弦一维插值(浮点版本,波形显示的优化插值)
	* @param[in]	nyin		输入数据的长度
	* @param[in]	dxin		x的间隔
	* @param[in]	fxin		x的起始值
	*
	* @param[in]	yin			输入数据的Y
	* @param[in]	yinl		y(x)的左极值（设置数值计算范围）
	* @param[in]	yinr		y(x)的右极值（设置数值计算范围）
	* @param[in]	nxout		输出数据长度
	* @param[in]	xout		输出的X数据
	*
	* @param[out]	yout		计算结果数据的数组【浮点型】
	*
	* @return	无
	*/
	static void ISL_ints8r (int nxin, float dxin, float fxin, float yin[],
							float yinl, float yinr, int nxout, float xout[], float yout[]);


/***********************************************************************************************
 *
 * 以上函数是主要调用函数
 *
 **************************************************************************************************/

public:
	/* 辅助函数  */
	/**
	* @brief	等间隔一维线性插值得输出X位置( ISL_intlin 的辅助函数 )
	*
	* @param[in]	nin			输入数据的长度
	* @param[in]	xin			输入数据的X
	* @param[in]	nout		输出数据长度
	*
	* @param[out]	xout		输出的X数据的数组【浮点型】
	*
	* @return	无
	*/
	static void ISL_intlin_getx (int nin, float *xin, int nout, float *xout);

	/**
	* @brief	8系数一维插值(浮点版本, 正弦插值的辅助函数)
	* @param[in]	ntable		插值算子的数目，ntable>=2
	* @param[in]	table[][8]	8点插值算子的序列
	*
	* @param[in]	nyin		输入数据的长度
	* @param[in]	dxin		x的间隔
	* @param[in]	fxin		x的起始值
	*
	* @param[in]	yin			输入数据的Y
	* @param[in]	yinl		y(x)的左极值（设置数值计算范围）
	* @param[in]	yinr		y(x)的右极值（设置数值计算范围）
	* @param[in]	nxout		输出数据长度
	* @param[in]	xout		输出的X数据
	*
	* @param[out]	yout		计算结果数据的数组【浮点型】
	*
	* @return	无
	*/
	static void ISL_intt8r (int ntable, float table[][8],
							int nxin, float dxin, float fxin, float yin[], float yinl, float yinr,
							int nxout, float xout[], float yout[]);

	/**
	* @brief	compute cubic spline interpolation coefficients for interpolation with continuous second derivatives
	* 			(计算连续的二阶导数的三次样条插值系数, 	和ISL_intcub配合使用)
	*
	* @param[in]	n			输入数据的长度
	* @param[in]	x[n]		输入数据的X
	* @param[in]	y[n]		array[nin][4] of y(x), y'(x), y''(x), and y'''(x)
	*
	* @param[out]	y[n][4]		计算结果数据的数组【浮点型】
	*
	* @return	无
	*/
	static void ISL_csplin ( int n, float x[], float y[], float yd[][4] );


public:
	/*单点一维插值方法*/

	/**
	* @brief	单点线性插值
	*
	* @param[in]	x[n]		输入数据的X
	* @param[in]	y[n]		输入数据的Y
	* @param[in]	n			输入数据的长度
	* @param[in]	t			添加的X定位点
	*
	* @return	double,  函数返回指定插值点t处的函数近似值 z=f（t）
	*/
	static float ISL_line(float *x, float *y, int n, float t);

	/**
	* @brief	一元三点插值,用抛物线插值公式计算指定插值点t处的函数近似值z=f（t）
	*
	* @param[in]	x[n]		输入数据的X
	* @param[in]	y[n]		输入数据的Y
	* @param[in]	n			输入数据的长度
	* @param[in]	t			添加的X定位点
	*
	* @return	double,  函数返回指定插值点t处的函数近似值 z=f（t）
	*/
	static float ISL_lg3(float *x, float *y, int n, float t);

	/**
	* @brief	拉格朗日（Lagrange）插值
	*
	* @param[in]	x[n]		输入数据的X
	* @param[in]	y[n]		输入数据的Y
	* @param[in]	n			输入数据的长度
	* @param[in]	t			添加的X定位点
	*
	* @return	double,  函数返回指定插值点t处的函数近似值 z=f（t）
	*/
	static float ISL_lgr(float *x, float *y, int n, float t);

	/**
	* @brief	连分式插值，用连分式插值公式计算指定插值点t处的函数近似值z=f（t）
	*
	* @param[in]	x[n]		输入数据的X
	* @param[in]	y[n]		输入数据的Y
	* @param[in]	n			输入数据的长度
	* @param[in]	t			添加的X定位点
	*
	* @return	double,  函数返回指定插值点t处的函数近似值 z=f（t）
	*/
	static float ISL_pqs(float *x, float *y, int n, float t);

	/**
	* @brief	埃特金插值，用埃特金（Aitken）插值公式计算指定插值点t处的函数近似值z=f（t）
	*
	* @param[in]	x[n]		输入数据的X
	* @param[in]	y[n]		输入数据的Y
	* @param[in]	n			输入数据的长度
	* @param[in]	t			添加的X定位点
	*
	* @return	double,  函数返回指定插值点t处的函数近似值 z=f（t）
	*/
	static float ISL_atk(float *x, float *y, int n, float t, float eps);

	// =============================================================================================
	//
	// =============================================================================================

	/**
	* @brief	埃尔米特插值，
	* 			给定n个结点Xi上的函数Yi=f（Xi）以及一阶导数值Yi'=f'（Xi），
	* 			用埃尔米特（Hermite）插值公式计算指定插值点t处的函数近似值z=f（t）
	*
	* @param[in]	x[n]		输入数据的X
	* @param[in]	y[n]		输入数据的Y
	* @param[in]	dy[n]		存放n个给定结点上的一阶导数值,y0’,y1’,...,yn-1’
	* @param[in]	n			输入数据的长度
	* @param[in]	t			添加的X定位点
	*
	* @return	double,  函数返回指定插值点t处的函数近似值 z=f（t）
	*/
	static float ISL_hmt(float *x, float *y, float *dy, int n, float t);


	/**
	* @brief	evaluate y(x), y'(x), y''(x), ... via piecewise cubic interpolation(浮点版本, 分段插值)
	* @param[in]	ntable		=0 if y(x) desired; =1 if y'(x) desired
	*
	* @param[in]	nin			输入数据的长度
	* @param[in]	xin			输入数据的X
	* @param[in]	ydin		array[nin][4] of y(x), y'(x), y''(x), and y'''(x) (参数从ISL_csplin函数获得)
	*
	* @param[in]	nout		输出数据长度
	* @param[in]	xout		输出的X数据
	*
	* @param[out]	yout		计算结果数据的数组【浮点型】
	*
	* @return	无
	*/
	static void ISL_intcub (	int ideriv,
								int nin, float xin[], float ydin[][4],
								int nout, float xout[], float yout[] );


};





}/* End Of namespace ISLIB */
#endif /* ISL_INTERPOLATION_H_ */
