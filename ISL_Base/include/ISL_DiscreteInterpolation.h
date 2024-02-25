/**
*	@file	ISL_DiscreteInterpolation.h
*	@brief	[Header file of Discrete Interpolation Functions], 离散插值；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-07-02
*	@refer
*/

#ifndef PAI_FRAME_ISEISLIB_DISCRETEINTERPOLATION_H_
#define PAI_FRAME_ISEISLIB_DISCRETEINTERPOLATION_H_


#include "ISL_UserDefine.h"

namespace ISLib {


/**
* @brief	二维反距离加权插值法\n
* 			near_num <= nin \n
* 			nin < noutx*nouty
*
* @param[in]	nin					给定样点的个数；
* @param[in]	xi[nin]				给定样点的X坐标；
* @param[in]	yi[nin]				给定样点的Y坐标；
* @param[in]	zi[nin]				给定样点的X,Y坐标位置的样点值；
* @param[in]	near_num			就近计算的样点个数（near_num <= nin）；
* @param[out]	noutx				输出数据X轴的个数；
* @param[out]	nouty				输出数据Y轴的个数；
* @param[out]	xout[noutx]			输出数据的X坐标；
* @param[out]	yout[nouty]			输出数据的X坐标；
* @param[out]	zout[noutx*nouty]	输出数据的样点值；
*
* @return	无
*/
ISL_EXPORT void ISL_dintIDW2D (	int nin,
						float *xi,
						float *yi,
						float *zi,
						int near_num,
						int noutx,
						int nouty,
						float *xout,
						float *yout,
						float *&zout );



/**
* @brief	三维反距离加权插值法\n
* 			near_num <= nin \n
* 			nin < noutx*nouty*noutz
*
* @param[in]	nin						给定样点的个数；
* @param[in]	xi[nin]					给定样点的X坐标；
* @param[in]	yi[nin]					给定样点的Y坐标；
* @param[in]	zi[nin]					给定样点的Z坐标；
* @param[in]	si[nin]					给定样点的X,Y,Z坐标位置的样点值；
* @param[in]	radius					搜索半径，小于noutx，nouty，noutz中的最小值 ；
* @param[in]	near_num				就近计算的样点个数（near_num <= nin）；
* @param[out]	noutx					输出数据X轴的个数；
* @param[out]	nouty					输出数据Y轴的个数；
* @param[out]	noutz					输出数据Z轴的个数；
* @param[out]	xout[noutx]				输出数据的X坐标；
* @param[out]	yout[nouty]				输出数据的Y坐标；
* @param[out]	zout[noutz]				输出数据的Z坐标；
* @param[out]	sout[noutx*nouty*noutz]	输出数据的样点值；
* @param[in]	mode					插值方式，0=常规，1=快速；
*
* @return	无
*/
ISL_EXPORT void ISL_dintIDW3D (	int nin,
						float *xi,
						float *yi,
						float *zi,
						float *si,
						float radius,
						int near_num,
						int noutx,
						int nouty,
						int noutz,
						float *xout,
						float *yout,
						float *zout,
						float *&sout,
						int mode = 1);





/**
* @brief	二维克里金插值法\n
* 			near_num <= nin \n
* 			nin < noutx*nouty
*
* @param[in]	*range				插值区域范围(矩形)，数组大小为4，array[4]；
* 									range[0] = 左上角 X
* 									range[1] = 左上角 Y
* 									range[2] = 右下角 X
* 									range[3] = 右下角 Y
* @param[in]	mode				计算半方差矩阵的 3三种模式
* 									1 = 球状模型
* 									2 = 指数模型
* 									3 = 高斯模型
* @param[in]	nin					给定样点的个数；
* @param[in]	nin					给定样点的个数；
* @param[in]	xi[nin]				给定样点的X坐标；
* @param[in]	yi[nin]				给定样点的Y坐标；
* @param[in]	zi[nin]				给定样点的X,Y坐标位置的样点值；
* @param[in]	c0					块金值(Nugget)；
* @param[in]	c1					基抬值(sill)；
* @param[in]	a					变程值(a)；
* @param[out]	noutx				输出数据X轴的个数；
* @param[out]	nouty				输出数据Y轴的个数；（noutx*nouty > nin）
* @param[out]	zout[noutx*nouty]	输出数据的样点值；
*
* @return	无
*/
ISL_EXPORT void ISL_dintKriging2D (	int *range,
						int mode,
						int nin,
						float *xi,
						float *yi,
						float *zi, /* 输入随即样点 */
						int noutx,
						int nouty,
						float *&zout/* array[noutx * nouty] */,
						float c0,
						float c1,
						float a );


}/* End Of namespace ISLIB */
#endif /* PAI_FRAME_ISEISLIB_DISCRETEINTERPOLATION_H_ */
