/**
*	@file	ISL_Geometric3D.h
*	@brief	[Header file of Geometric 3D Functions], 三维空间计算几何；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-07-02
*	@refer
*/

#ifndef PAI_FRAME_ISEISLIB_GEOMETRIC3D_H_
#define PAI_FRAME_ISEISLIB_GEOMETRIC3D_H_


#include "ISL_UserDefine.h"

namespace ISLib {

/**
*	@brief	结构类型	三维维空间点结构
*/
struct ISL_EXPORT POINT3D
{
	float x;	/**< 点的X坐标 */
	float y;	/**< 点的Y坐标 */
	float z;	/**< 点的Y坐标 */

	POINT3D(float a = 0, float b = 0, float c = 0) {
		x = a;
		y = b;
		z = c;
	}
};


/**
* @brief	三维空间两点之间的距离
*
* @param[in]	p1			点1；
* @param[in]	p2			点2；
*
* @return	返回两点之间的距离（float）
*/
ISL_EXPORT float ISL_p2pDist3D(POINT3D p1, POINT3D p2);


}/* End Of namespace ISLIB */




#endif /* PAI_FRAME_ISEISLIB_GEOMETRIC3D_H_ */
