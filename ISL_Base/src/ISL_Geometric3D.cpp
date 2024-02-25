/**
*	@file	ISL_Geometric3D.cpp
*	@brief	[Source file of Geometric 3D Functions], 三维空间计算几何；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-07-02
*	@refer
*/



#include "ISL_UserDefine.h"
#include "ISL_Geometric3D.h"

namespace ISLib {

float ISL_p2pDist3D(POINT3D p1, POINT3D p2)
{
	return (sqrt(
				(p1.x - p2.x) * (p1.x - p2.x) +
				(p1.y - p2.y) * (p1.y - p2.y) +
				(p1.z - p2.z) * (p1.z - p2.z)
				)
			);
}



}/* End Of namespace ISLIB */
