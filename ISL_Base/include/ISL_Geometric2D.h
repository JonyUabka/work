/**
*	@file	ISL_Geometric2D.h
*	@brief	[Header file of Geometric Functions], 计算几何；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-06-24
*	@refer	NEWS
*/

#ifndef PAI_FRAME_ISEISLIB_GEOMETRIC2D_H_
#define PAI_FRAME_ISEISLIB_GEOMETRIC2D_H_


#include "ISL_UserDefine.h"

namespace ISLib {


/**
*	@brief	结构类型	二维空间点结构
*/
struct ISL_EXPORT POINT
{
	float x;	/**< 点的X坐标 */
	float y;	/**< 点的Y坐标 */

	POINT(float a = 0, float b = 0) {
		x = a;
		y = b;
	}
};


/**
*	@brief	结构类型	二维空间线段（由两点组成的一条线）
*/
struct ISL_EXPORT LINESEG
{
	POINT s;	/**< 起始点 */
	POINT e;	/**< 终止点 */

	LINESEG(POINT a, POINT b){
		s = a;
		e = b;
	}
	LINESEG(){}
};


/**
*	@brief	结构类型	二维空间直线方程参数
*
*	a*x+b*y+c=0 为统一表示，约定 a>= 0
*/
struct ISL_EXPORT LINE
{
	float a;	/**< 参数a */
	float b;	/**< 参数b */
	float c;	/**< 参数c */

	LINE(float d1 = 1, float d2 = -1, float d3 = 0){
		a = d1;
		b = d2;
		c = d3;
	}
};


/**
* @brief	两点之间的距离
*
* @param[in]	p1			点1；
* @param[in]	p2			点2；
*
* @return	返回两点之间的距离（float）
*/
ISL_EXPORT float ISL_p2pDist(POINT p1, POINT p2);


/**
* @brief	判断两个点是否重合
*
* @param[in]	p1			点1；
* @param[in]	p2			点2；
*
* @return	返回true为重合，false为不重合（bool）
*/
ISL_EXPORT bool ISL_equalPoint(POINT p1, POINT p2);


/**
* @brief	矢量叉乘
*
* @param[in]	sp			起始点的X,Y坐标；
* @param[in]	ep			终了点的X,Y坐标；
* @param[in]	op			原点的X,Y坐标；
*
* @return			返回矢量叉乘的值r（float）
*  					r>0:sp在矢量op ep的顺时针方向；
*  					r=0:op sp ep三点共线；
*  					r<0:sp在矢量op ep的逆时针方向 。
*/
ISL_EXPORT float ISL_multiply(POINT sp, POINT ep, POINT op);


/**
* @brief	绝对值矢量叉乘
*
* @param[in]	sp			起始点的X,Y坐标；
* @param[in]	ep			终了点的X,Y坐标；
* @param[in]	op			原点的X,Y坐标；
*
* @return		返回返回矢量叉乘的绝对值（float）
*/
ISL_EXPORT float ISL_absMultiply(POINT sp, POINT ep, POINT op);


/**
* @brief	矢量(p1-op)和(p2-op)的点积
*
* @param[in]	p1			点p1的X,Y坐标；
* @param[in]	p2			点p2的X,Y坐标；
* @param[in]	p0			原点的X,Y坐标；
*
* @return			返回矢量点积的值r（float）
* 					如果两个矢量都非零矢量
*  					r < 0: 两矢量夹角为锐角；
*  					r = 0: 两矢量夹角为直角；
*  					r > 0: 两矢量夹角为钝角 。
*/
ISL_EXPORT float ISL_dotMultiply(POINT p1, POINT p2, POINT p0);


/**
* @brief	判断点是否在矩形内(起始点、矩形宽、高)
*
* @param[in]	pLT			矩形左上角的X,Y坐标；
* @param[in]	w			矩形的宽度；
* @param[in]	h			矩形的高度；
* @param[in]	p0			给定点的X,Y坐标；
*
* @return		bool，true为在矩形内，false为不在矩形内
*/
ISL_EXPORT bool ISL_onRect(POINT pLT, float w, float h, POINT p0);


/**
* @brief	判断点是否在矩形内(矩形左上角点、矩形右下角点)
*
* @param[in]	pLT			矩形左上角的X,Y坐标；
* @param[in]	pRB			矩形右下角的X,Y坐标；
* @param[in]	p0			给定点的X,Y坐标；
*
* @return		bool，true为在矩形内，false为不在矩形内
*/
ISL_EXPORT bool ISL_onRect(POINT pLT, POINT pRB, POINT p0);


/**
* @brief	判断点是否在矩形内(线段l为对角线的矩形)
*
* @param[in]	l			线段l为对角线的矩形；
* @param[in]	p0			给定点的X,Y坐标；
*
* @return		返回bool型，true为在矩形内，false为不在矩形内
*/
ISL_EXPORT bool ISL_onRect(LINESEG l, POINT p);


/**
* @brief	判断点p是否在线段l上
*
* 			条件：(p在线段l所在的直线上) && (点p在以线段l为对角线的矩形内)
*
* @param[in]	l			一条线段的坐标；
* @param[in]	p0			给定点的X,Y坐标；
*
* @return		返回bool型，true为在线段上，false为不在线段上
*/
ISL_EXPORT bool ISL_onLine(LINESEG l, POINT p);



/**
* @brief	判断点q是否在多边形polygon内。
* 			利用射线法判断点q与多边形polygon的位置关系, 要求polygon为简单多边形，顶点逆时针排列
*
* @param[in]	polygon		多边形的结点数组；
* @param[in]	num			多边形结点的个数；
* @param[in]	q			点q的坐标；
*
* @return		int，如果点在多边形内： 返回0
* 					如果点在多边形边上：返回1
* 					如果点在多边形外： 返回2。
*/
ISL_EXPORT int ISL_inPolygon(POINT * polygon, int num, POINT q);



/**
* @brief	多边形面积； 输入顶点按逆时针排列时，返回正值；否则返回负值。
*
* @param[in]	polygon		多边形的结点数组；
* @param[in]	num			多边形结点的个数；
*
* @return		float,返回面积大小，输入顶点按逆时针排列时，返回正值；否则返回负值。
*/
ISL_EXPORT float ISL_polygonArea(POINT * polygon, int num);



/**
* @brief	点p以点o为圆心逆时针旋转alpha(单位：弧度)后所在的位置
*
* @param[in]	o			圆心点的坐标；
* @param[in]	p			初始指针末端坐标；
* @param[in]	alpha		旋转的角度；
*
* @return		返回点p以点o为圆心逆时针旋转alpha(单位：弧度)后所在的位置
*/
ISL_EXPORT POINT ISL_rotate(POINT o, double alpha, POINT p);



/**
* @brief	返回顶角在o点，起始边为os，终止边为oe的夹角(单位：弧度),可以用于求线段之间的夹角
*
* @param[in]	o			圆心点的坐标；
* @param[in]	os			起始边s指针末端坐标；
* @param[in]	oe			终止边e指针末端坐标；
*
* @return		角度小于pi，返回正值；角度小于pi，返回正值；
*/
ISL_EXPORT float ISL_angle(POINT o, POINT os, POINT oe);



/**
* @brief	判断点C在线段AB所在的直线l上垂足P的与线段AB的关系。
*
* 			本函数是根据下面的公式写的，P是点C到线段AB所在直线的垂足
*			AC dot AB
*			r = ----------------------
*			||AB||^2
*			(Cx-Ax)(Bx-Ax) + (Cy-Ay)(By-Ay)
*			= ----------------------------------------------------
*			L^2
*
*			r has the following meaning:
*			r = 0, P = A
*			r = 1, P = B
*			r < 0, P is on the backward extension of AB
*			r > 1, P is on the forward extension of AB
*			0 < r < 1, P is interior to AB
*
* @param[in]	c			给定点的坐标；
* @param[in]	l			线段l的坐标；
*
* @return		r值(float)
*/
ISL_EXPORT float ISL_relation(POINT c, LINESEG l);



/**
* @brief	求点C到线段AB所在直线的垂足P
*
* @param[in]	p			给定点的坐标；
* @param[in]	l			线段l的坐标；
*
* @return		返回垂足的坐标；
*/
ISL_EXPORT POINT ISL_perpendicularFoot(POINT p, LINESEG l);



/**
* @brief	求点p到线段l的最短距离, 返回线段上距该点最近的点np 注意：np是线段l上到点p最近的点，不一定是垂足。
*
* @param[in]	p			给定点的坐标；
* @param[in]	l			线段l的坐标；
* @param[out]	np			返回线段上距该点最近的点的坐标；
*
* @return		返回p到线段l的最短距离；
*/
ISL_EXPORT float ISL_p2lDist(POINT p, LINESEG l, POINT &np);



/**
* @brief	求点p到线段l所在直线的距离,
* 			请注意本函数与[ISL_p2lDist(POINT p, LINESEG l, POINT &np)]函数的区别。
*
* @param[in]	p			给定点的坐标；
* @param[in]	l			线段l的坐标；
*
* @return		返回p到线段l所在直线的距离；
*/
ISL_EXPORT float ISL_p2lDist(POINT p, LINESEG l);



/**
* @brief	计算点到折线集的最近距离,并返回最近点, 注意：调用的是ISL_p2lDist(POINT p, LINESEG l, POINT &np)函数。
*
* @param[in]	curve[num]		折线样点；
* @param[in]	num				折线样点个数；
* @param[in]	p				给定点的坐标；
* @param[in]	q				返回最近的点的坐标；
*
* @return		返回p到线段l所在直线的距离；
*/
ISL_EXPORT float ISL_p2curve(POINT * curve, int num, POINT p, POINT &q);



/**
* @brief	判断圆是否在多边形内。
*
* @param[in]	center			远点的坐标；
* @param[in]	radius			圆的半径；
* @param[in]	polygon[num]	多边形结点数组；
* @param[in]	num				多边形结点的个数；
*
* @return		返回布尔型，true：在多边形内，false：在多边形外；
*/
ISL_EXPORT bool ISL_circleInPolygon(POINT center, float radius, POINT polygon[], int num);



/**
* @brief	返回两个矢量l1和l2的夹角的余弦 (-1 ~ 1), 注意：如果想从余弦求夹角的话，注意反余弦函数的值域是从 0到PI。
*
* @param[in]	l1			线段l1的坐标；
* @param[in]	l2			线段l2的坐标；
*
* @return		float，两个矢量l1和l2的夹角的余弦 (-1 ~ 1)；
*/
ISL_EXPORT float ISL_cosine(LINESEG l1, LINESEG l2);



/**
* @brief	返回线段l1与l2之间的夹角,单位：弧度 范围(-pi，pi)。
*
* @param[in]	l1			线段l1的坐标；
* @param[in]	l2			线段l2的坐标；
*
* @return		float，返回线段l1与l2之间的夹角,单位：弧度 范围(-pi，pi)；
*/
ISL_EXPORT float ISL_lsAngle(LINESEG l1, LINESEG l2);


/**
* @brief	判断线段u和v相交(包括相交在端点处)。
*
* @param[in]	u			线段l1的坐标；
* @param[in]	v			线段l2的坐标；
*
* @return		bool，值为true，表示相交；值为false，表示不相交。
*/
ISL_EXPORT bool ISL_interSect(LINESEG u, LINESEG v);


/**
* @brief	判断线段u和v相交（不包括双方的端点）。
*
* @param[in]	u			线段l1的坐标；
* @param[in]	v			线段l2的坐标；
*
* @return		bool，值为true，表示相交；值为false，表示不相交。
*/
ISL_EXPORT bool ISL_interSect_A(LINESEG u, LINESEG v);


/**
* @brief	判断线段v所在直线与线段u相交,方法：判断线段u是否跨立线段v。
*
* @param[in]	u			线段l1的坐标；
* @param[in]	v			线段l2的坐标；
*
* @return		bool，值为true，表示相交；值为false，表示不相交。
*/
ISL_EXPORT bool ISL_interSect_l(LINESEG u, LINESEG v);


/**
* @brief	根据已知两点坐标，求过这两点的直线解析方程： a*x+b*y+c = 0 (a >= 0)。
*
* @param[in]	p1			p1点的坐标；
* @param[in]	p2			p2点的坐标；
*
* @return		LINE，返回直线方程。
*/
ISL_EXPORT LINE ISL_makeLine(POINT p1, POINT p2);


/**
* @brief	根据直线解析方程返回直线的斜率k, 水平线返回 0, 竖直线返回 1e200。
*
* @param[in]	l			直线方程；
*
* @return		float，返回直线的斜率k, 水平线返回 0, 竖直线返回 1e200。
*/
ISL_EXPORT float ISL_slope(LINE l);



/**
* @brief	返回直线的倾斜角alpha ( 0 - pi),注意：atan()返回的是 -PI/2 ~ PI/2。
*
* @param[in]	l			直线方程；
*
* @return		float，返回直线的倾斜角alpha ( 0 - pi)。
*/
ISL_EXPORT float ISL_alpha(LINE l);



/**
* @brief	求点p关于直线l的对称点。
*
* @param[in]	l			直线方程；
* @param[in]	p			给定的点的坐标；
*
* @return		POINT， 返回p关于直线l的对称点。
*/
ISL_EXPORT POINT ISL_symmetry(LINE l, POINT p);


/**
* @brief	直线相交
* 			如果两条直线 l1(a1*x+b1*y+c1 = 0),
* 			l2(a2*x+b2*y+c2 = 0)相交，
* 			返回true，且返回交点p。
*
* @param[in]	l1			直线方程l1；
* @param[in]	l2			直线方程l2；
* @param[out]	p			返回的交点的坐标；
*
* @return		bool，值为true，表示相交；值为false，表示不相交。
*/
ISL_EXPORT bool ISL_lineInterSect(LINE l1, LINE l2, POINT &p);



/**
* @brief	线段相交
* 			如果线段l1和l2相交，
* 			返回true且交点由(inter)返回，否则返回false。
*
* @param[in]	l1			线段l1坐标；
* @param[in]	l2			线段l2坐标；
* @param[out]	p			返回的交点的坐标；
*
* @return		bool,返回true且交点由(p)返回，否则返回false。
*/
ISL_EXPORT bool ISL_lineSegInterSect(LINESEG l1, LINESEG l2, POINT &p);



/**
* @brief	判断顶点是否按逆时针排列,如果输入顶点按逆时针排列，返回true。
*
* @param[in]	polygon		多边形的结点数组；
* @param[in]	num			多边形结点的个数；
*
* @return		bool，输入顶点按逆时针排列，返回true，顺时针排序，返回false。
*/
ISL_EXPORT float ISL_isConterClock(POINT * polygon, int num);


/**
* @brief	判断点q在凸多边形polygon内。
* 			点q是凸多边形polygon内[包括边上]时，返回true
* 			注意：多边形polygon一定要是凸多边形
*
* @param[in]	polygon		多边形的结点数组；
* @param[in]	num			多边形结点的个数；
* @param[in]	q			点q的坐标；
*
* @return		bool，内[包括边上]时，返回true；不在其中，返回false。
*/
ISL_EXPORT bool ISL_inConvexPolygon(POINT * polygon, int num, POINT q);



/**
* @brief	寻找凸包的graham 扫描法。
* 			点q是凸多边形polygon内[包括边上]时，返回true
* 			注意：多边形polygon一定要是凸多边形
*
* @param[in]	points		输入的点集数组；
* @param[in]	pnum		输入的点集数组的个数；
* @param[out]	ch			输出的凸包上的点集，按照逆时针方向排列；
* @param[out]	cnum		输出的凸包上的点的个数；
*
* @return		无。
*/
ISL_EXPORT void ISL_grahamScan(POINT * points, int pnum, POINT *ch, int &cnum);




/**
* @brief	卷包裹法求点集凸壳。
*
* @param[in]	points		输入的点集数组；
* @param[in]	pnum		输入的点集数组的个数；
* @param[out]	ch			输出的凸包上的点集，按照逆时针方向排列；
* @param[out]	cnum		输出的凸包上的点的个数；
*
* @return		无。
*/
ISL_EXPORT void ISL_convexClosure(POINT *points, int pnum, POINT *ch, int &cnum);



/**
* @brief	求凸多边形的重心,要求输入多边形按逆时针排序。
*
* @param[in]	polygon		输入的点集数组；
* @param[in]	num			输入的点集数组的个数；
*
* @return		POINT，重心点的坐标。
*/
ISL_EXPORT POINT ISL_gravityCenter(POINT *polygon, int num);


} /* End Of namespace ISLIB */
#endif /* PAI_FRAME_ISEISLIB_GEOMETRIC2D_H_ */
