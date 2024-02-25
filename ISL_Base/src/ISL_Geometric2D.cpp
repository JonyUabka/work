//*
//*	@file	ISL_Geometric2D.cpp
//*	@brief	[Header file of ISL_Geometric Functions], 计算几何；
//*	@see	ISeisLib Manual
//*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
//*	@date	2014-06-24
//*	@refer	Baidu
//*/


//#include "ISL_Geometric2D.h"

//namespace ISLib {

///* 两点之间的距离  */
//float ISL_p2pDist(POINT p1, POINT p2)
//{
//	return (sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y)));
//}

///* 判断两个点是否重合  */
//bool ISL_equalPoint(POINT p1, POINT p2)
//{
//	return ((fabs(p1.x - p2.x) < EP) && (fabs(p1.y - p2.y) < EP));
//}

///* 矢量叉乘 */
//float ISL_multiply(POINT sp, POINT ep, POINT op)
//{
//	return ((sp.x - op.x) * (ep.y - op.y) - (ep.x - op.x) * (sp.y - op.y));
//}

///* 矢量叉乘的绝对值  */
//float ISL_absMultiply(POINT sp, POINT ep, POINT op)
//{
//	return fabs((sp.x - op.x) * (ep.y - op.y) - (ep.x - op.x) * (sp.y - op.y));
//}

///* 矢量点积  */
//float ISL_dotMultiply(POINT p1, POINT p2, POINT p0)
//{
//	return ((p1.x - p0.x) * (p2.x - p0.x) + (p1.y - p0.y) * (p2.y - p0.y));
//}

///*判断点是否在矩形内(起始点、矩形宽、高)*/
//bool ISL_onRect(POINT pLT, float w, float h, POINT p0)
//{
//	w += pLT.x;
//	h += pLT.y;

//	if (p0.x < pLT.x || p0.x > w)
//		return false;
//	if (p0.y < pLT.y || p0.y > h)
//		return false;

//	return true;
//}

///*判断点是否在矩形内(矩形左上角点、矩形右下角点)*/
//bool ISL_onRect(POINT pLT, POINT pRB, POINT p0)
//{
//	float w = pRB.x - pLT.x;
//	float h = pRB.y - pLT.y;

//	return ISL_onRect(pLT, w, h, p0);
//}

///*判断点是否在矩形内(线段l为对角线的矩形)*/
//bool ISL_onRect(LINESEG l, POINT p)
//{
//	POINT pLT(l.s.x, l.s.y), pRB(l.e.x, l.e.y);
//	return ISL_onRect(pLT, pRB, p);
//}

///*判断点p是否在线段l上*/
//bool ISL_onLine(LINESEG l, POINT p)
//{
//	return ( (ISL_multiply(l.e, p, l.s) == 0 )
//			&&
//			ISL_onRect(l, p)
//			);
//}

///*点p以点o为圆心逆时针旋转alpha(单位：弧度)后所在的位置*/
//POINT ISL_rotate(POINT o, double alpha, POINT p)
//{
//	POINT tp;
//	p.x -= o.x;
//	p.y -= o.y;
//	tp.x = p.x * cos(alpha) - p.y * sin(alpha) + o.x;
//	tp.y = p.y * cos(alpha) + p.x * sin(alpha) + o.y;

//	return tp;
//}

///*返回顶角在o点，起始边为os，终止边为oe的夹角(单位：弧度)，可以用于求线段之间的夹角*/
//float ISL_angle(POINT o, POINT os, POINT oe)
//{
//	float cosfi, fi, norm;
//	float dsx = os.x - o.x;
//	float dsy = os.y - o.y;
//	float dex = oe.x - o.x;
//	float dey = oe.y - o.y;

//	cosfi = dsx * dex + dsy * dey;
//	norm = (dsx * dsx + dey * dey) * (dex * dex + dey * dey);
//	cosfi /= sqrt(norm);
//	if (cosfi >= 1.0)
//		return 0;
//	if (cosfi <= -1.0)
//		return -3.1415926;

//	fi = acos(cosfi);
//	if (dsx * dey - dsy * dex > 0)
//		return fi; // 说明矢量os 在矢量 oe的顺时针方向

//	return -fi;
//}

///*判断点C在线段AB所在的直线l上垂足P的与线段AB的关系*/
//float ISL_relation(POINT c, LINESEG l)
//{
//	LINESEG tl;
//	tl.s = l.s;
//	tl.e = c;
//	return ISL_dotMultiply(tl.e, l.e, l.s)
//			/ (ISL_p2pDist(l.s, l.e) * ISL_p2pDist(l.s, l.e));
//}

///*求点C到线段AB所在直线的垂足P*/
//POINT ISL_perpendicularFoot(POINT p, LINESEG l)
//{
//	double r = ISL_relation(p, l);

//	POINT tp;
//	tp.x = l.s.x + r * (l.e.x - l.s.x);
//	tp.y = l.s.y + r * (l.e.y - l.s.y);

//	return tp;
//}


///*求点p到线段l的最短距离, 返回线段上距该点最近的点np 注意：np是线段l上到点p最近的点，不一定是垂足*/
//float ISL_p2lDist(POINT p, LINESEG l, POINT &np)
//{
//	float r = ISL_relation(p, l);

//	if (r < 0) {
//		np = l.s;
//		return ISL_p2pDist(p, l.s);
//	}

//	if (r > 1) {
//		np = l.e;
//		return ISL_p2pDist(p, l.e);
//	}

//	np = ISL_perpendicularFoot(p, l);

//	return ISL_p2pDist(p, np);
//}

///*求点p到线段l所在直线的距离,请注意本函数与上个函数的区别*/
//float ISL_p2lDist(POINT p, LINESEG l)
//{
//	return fabs(ISL_multiply(p, l.e, l.s)) / ISL_p2pDist(l.s, l.e);
//}

///*计算点到折线集的最近距离,并返回最近点*/
//float ISL_p2curve(POINT * curve, int num, POINT p, POINT &q)
//{
//	int i;

//	float cd = (float) (INF), td;
//	LINESEG l;
//	POINT tq, cq;

//	for (i = 0; i < num - 1; i++) {
//		l.s = curve[i];
//		l.e = curve[i + 1];
//		td = ISL_p2lDist(p, l, tq);

//		if (td < cd) {
//			cd = td;
//			cq = tq;
//		}
//	}
//	q = cq;

//	return cd;
//}


///*判断圆是否在多边形内*/
//bool ISL_circleInPolygon(POINT center, float radius, POINT * polygon, int num)
//{
//	POINT q;
//	float d;

//	q.x = 0;
//	q.y = 0;

//	d = ISL_p2curve(polygon, num, center, q);

//	if (d < radius || fabs(d - radius) < EP)
//		return true;
//	else
//		return false;
//}


///*返回两个矢量l1和l2的夹角的余弦 (-1 ~ 1), 注意：如果想从余弦求夹角的话，注意反余弦函数的值域是从 0到PI*/
//float ISL_cosine(LINESEG l1, LINESEG l2)
//{
//	float A1 = (l1.e.x - l1.s.x) * (l2.e.x - l2.s.x);
//	float A2 = (l1.e.y - l1.s.y) * (l2.e.y - l2.s.y);
//	float D1 = ISL_p2pDist(l1.e, l1.s);
//	float D2 = ISL_p2pDist(l2.e, l2.s);

//	return (A1 + A2) / (D1 * D2);
//}

///*线段l1与l2之间的夹角,单位：弧度 范围(-pi，pi)*/
//float ISL_lsAngle(LINESEG l1, LINESEG l2)
//{
//	POINT o, s, e;

//	o.x = o.y = 0;
//	s.x = l1.e.x - l1.s.x;
//	s.y = l1.e.y - l1.s.y;
//	e.x = l2.e.x - l2.s.x;
//	e.y = l2.e.y - l2.s.y;

//	return ISL_angle(o, s, e);
//}

///*判断线段u和v相交(包括相交在端点处)*/
//bool ISL_interSect(LINESEG u, LINESEG v)
//{
//	return ((MAX(u.s.x, u.e.x) >= MIN(v.s.x, v.e.x)) && //排斥实验
//			(MAX(v.s.x, v.e.x) >= MIN(u.s.x, u.e.x)) && (MAX(u.s.y, u.e.y)
//			>= MIN(v.s.y, v.e.y)) && (MAX(v.s.y, v.e.y) >= MIN(u.s.y, u.e.y))
//			&& (ISL_multiply(v.s, u.e, u.s) * ISL_multiply(u.e, v.e, u.s) >= 0) && //跨立实验
//			(ISL_multiply(u.s, v.e, v.s) * ISL_multiply(v.e, u.e, v.s) >= 0));
//}


///* 判断线段u和v相交（不包括双方的端点）*/
//bool ISL_interSect_A(LINESEG u, LINESEG v)
//{
//	return ((ISL_interSect(u, v)) && (!ISL_onLine(u, v.s))
//			&& (!ISL_onLine(u, v.e)) && (!ISL_onLine(v, u.e))
//			&& (!ISL_onLine(v, u.s)));
//}

///* 判断线段v所在直线与线段u相交,判断线段v所在直线与线段u相交(延长线相交)*/
//bool ISL_interSect_l(LINESEG u, LINESEG v)
//{
//	return ISL_multiply(u.s, v.e, v.s) * ISL_multiply(v.e, u.e, v.s) >= 0;
//}

///*根据已知两点坐标，求过这两点的直线解析方程： a*x+b*y+c = 0 (a >= 0)*/
//LINE ISL_makeLine(POINT p1, POINT p2)
//{
//	LINE tl;
//	int sign = 1;

//	tl.a = p2.y - p1.y;

//	if (tl.a < 0) {
//		sign = -1;
//		tl.a = sign * tl.a;
//	}
//	tl.b = sign * (p1.x - p2.x);
//	tl.c = sign * (p1.y * p2.x - p1.x * p2.y);

//	return tl;
//}

///*根据直线解析方程返回直线的斜率k, 水平线返回 0, 竖直线返回 1e200*/
//float ISL_slope(LINE l)
//{
//	if (fabs(l.a) < 1e-20)
//		return 0;
//	if (fabs(l.b) < 1e-20)
//		return INF;

//	return -(l.a / l.b);
//}

///*返回直线的倾斜角alpha ( 0 - pi),注意：atan()返回的是 -PI/2 ~ PI/2*/
//float ISL_alpha(LINE l)
//{
//	if (fabs(l.a) < EP)
//		return 0;
//	if (fabs(l.b) < EP)
//		return PI / 2;

//	float k = ISL_slope(l);

//	if (k > 0)
//		return atan(k);
//	else
//		return PI + atan(k);
//}

///*求点p关于直线l的对称点*/
//POINT ISL_symmetry(LINE l, POINT p)
//{
//	POINT tp;
//	tp.x = ((l.b * l.b - l.a * l.a) * p.x - 2 * l.a * l.b * p.y - 2 * l.a * l.c)
//			/ (l.a * l.a + l.b * l.b);

//	tp.y = ((l.a * l.a - l.b * l.b) * p.y - 2 * l.a * l.b * p.x - 2 * l.b * l.c)
//			/ (l.a * l.a + l.b * l.b);

//	return tp;
//}

///*如果两条直线 l1(a1*x+b1*y+c1 = 0), l2(a2*x+b2*y+c2 = 0)相交，返回true，且返回交点p*/
//bool ISL_lineInterSect(LINE l1, LINE l2, POINT &p)
//{
//	float d = l1.a * l2.b - l2.a * l1.b;

//	if (fabs(d) < EP) // 不相交
//		return false;

//	p.x = (l2.c * l1.b - l1.c * l2.b) / d;
//	p.y = (l2.a * l1.c - l1.a * l2.c) / d;

//	return true;
//}

///*如果线段l1和l2相交，返回true且交点由(p)返回，否则返回false*/
//bool ISL_lineSegInterSect(LINESEG l1, LINESEG l2, POINT &p)
//{
//	LINE ll1, ll2;

//	ll1 = ISL_makeLine(l1.s, l1.e);
//	ll2 = ISL_makeLine(l2.s, l2.e);

//	if (ISL_lineInterSect(ll1, ll2, p))
//		return ISL_onLine(l1, p);
//	else
//		return false;
//}

///*多边形面积； 输入顶点按逆时针排列时，返回正值；否则返回负值*/
//float ISL_polygonArea(POINT polygon[], int num)
//{
//	int i;
//	float s;

//	if (num < 3)
//		return 0;

//	s = polygon[0].y * (polygon[num - 1].x - polygon[1].x);

//	for (i = 1; i < num; i++) {
//		s += polygon[i].y * (polygon[(i - 1)].x - polygon[(i + 1) % num].x);
//	}

//	return s / 2;
//}

///*判断顶点是否按逆时针排列,如果输入顶点按逆时针排列，返回true*/
//float ISL_isConterClock(POINT * polygon, int num)
//{
//	return ISL_polygonArea(polygon, num) > 0;
//}

///*判断点q是否在多边形polygon内*/
//int ISL_inPolygon(POINT * polygon, int num, POINT q)
//{
//	int c = 0, i, n;
//	LINESEG l1, l2;

//	l1.s = q;
//	l1.e = q;
//	l1.e.x = float(INF);

//	n = num;
//	for (i = 0; i < num; i++) {
//		l2.s = polygon[i];
//		l2.e = polygon[(i + 1) % num];
////		float ee = polygon[(i + 2) % num].x;
////		float ss = polygon[(i + 3) % num].y;

//		if (ISL_onLine(l2, q))
//			return 1;

//		if (ISL_interSect_A(l1, l2))
//			c++; // 相交且不在端点

//		if (ISL_onLine(l1, l2.e) && !ISL_onLine(l1, l2.s) && l2.e.y > l2.e.y)
//			c++;// l2的一个端点在l1上且该端点是两端点中纵坐标较大的那个

//		if (!ISL_onLine(l1, l2.e) && ISL_onLine(l1, l2.s) && l2.e.y < l2.e.y)
//			c++;// 忽略平行边
//	}
//	if (c % 2 == 1)
//		return 0;
//	else
//		return 2;
//}


///*判断点q在凸多边形polygon内*/
//bool ISL_inConvexPolygon(POINT * polygon, int num, POINT q)
//{
//	POINT p;
//	LINESEG l;
//	int i;

//	p.x = 0;
//	p.y = 0;
//	for (i = 0; i < num; i++) {
//		// 寻找一个肯定在多边形polygon内的点p：多边形顶点平均值
//		p.x += polygon[i].x;
//		p.y += polygon[i].y;
//	}
//	p.x /= num;
//	p.y /= num;

//	for (i = 0; i < num; i++) {
//		l.s = polygon[i];
//		l.e = polygon[(i + 1) % num];
//		if (ISL_multiply(p, l.e, l.s) * ISL_multiply(q, l.e, l.s) < 0)/* 点p和点q在边l的两侧，说明点q肯定在多边形外 */
//			return false;
//	}

//	return true;
//}


///*寻找凸包的graham 扫描法*/
//void ISL_grahamScan(POINT * points, int pnum, POINT ch[], int &cnum)
//{
//	int i, j, k = 0, top = 2;
//	POINT tmp;
//	// 选取PointSet中y坐标最小的点PointSet[k]，如果这样的点有多个，则取最左边的一个
//	for (i = 1; i < pnum; i++)
//		if (points[i].y < points[k].y || ((points[i].y == points[k].y)
//				&& (points[i].x < points[k].x)))
//			k = i;

//	tmp = points[0];
//	points[0] = points[k];
//	points[k] = tmp; // 现在PointSet中y坐标最小的点在PointSet[0]

//	/*
//	 * 对顶点按照相对PointSet[0]的极角从小到大进行排序，
//	 * 极角相同 的按照距离PointSet[0]从近到远进行排序
//	 * */
//	for (i = 1; i < pnum - 1; i++) {
//		k = i;
//		for (j = i + 1; j < pnum; j++)
//			if (ISL_multiply(points[j], points[k], points[0]) > 0
//					|| /*极角更小 */((ISL_multiply(points[j], points[k], points[0]) == 0)
//							&& /*极角相等，距离更短 */ISL_p2pDist(points[0], points[j]))< ISL_p2pDist(points[0], points[k]))
//				k = j;

//		tmp = points[i];
//		points[i] = points[k];
//		points[k] = tmp;
//	}

//	ch[0] = points[0];
//	ch[1] = points[1];
//	ch[2] = points[2];
//	for (i = 3; i < pnum; i++) {
//		while (ISL_multiply(points[i], ch[top], ch[top - 1]) >= 0)
//			top--;
//		ch[++top] = points[i];
//	}

//	cnum = top + 1;
//}


///*卷包裹法求点集凸壳*/
//void ISL_convexClosure(POINT * points, int pnum, POINT ch[], int &cnum)
//{
//	int top = 0, i, index, first;
//	float curmax, curcos, curdis;
//	POINT tmp;
//	LINESEG l1, l2;
//	bool use[MAXV];
//	tmp = points[0];
//	index = 0;

//	// 选取y最小点，如果多于一个，则选取最左点
//	for (i = 1; i < pnum; i++) {
//		if ((points[i].y < tmp.y || points[i].y) == (tmp.y && points[i].x < tmp.x)) {
//			index = i;
//		}
//		use[i] = false;
//	}
//	tmp = points[index];
//	first = index;
//	use[index] = true;
//	index = -1;
//	ch[top++] = tmp;
//	tmp.x -= 100;
//	l1.s = tmp;
//	l1.e = ch[0];
//	l2.s = ch[0];

//	while (index != first) {
//		curmax = -100;
//		curdis = 0;

//		// 选取与最后一条确定边夹角最小的点，即余弦值最大者
//		for (i = 0; i < pnum; i++) {
//			if (use[i])
//				continue;
//			l2.e = points[i];
//			curcos = ISL_cosine(l1, l2); // 根据cos值求夹角余弦，范围在 （-1 -- 1 ）
//			if ((curcos > curmax) || ((fabs(curcos - curmax) < 1e-6)
//					&& (ISL_p2pDist(l2.s, l2.e) > curdis))) {
//				curmax = curcos;
//				index = i;
//				curdis = ISL_p2pDist(l2.s, l2.e);
//			}
//		}

//		use[first] = false; //清空第first个顶点标志，使最后能形成封闭的hull
//		use[index] = true;
//		ch[top++] = points[index];
//		l1.s = ch[top - 2];
//		l1.e = ch[top - 1];
//		l2.s = ch[top - 1];
//	}

//	cnum = top - 1;
//}


///*求凸多边形的重心,要求输入多边形按逆时针排序*/
//POINT ISL_gravityCenter(POINT * polygon, int num)
//{
//	POINT tp;
//	float x, y, s, x0, y0, cs, k;

//	x = 0;
//	y = 0;
//	s = 0;

//	for (int i = 1; i < num - 1; i++) {
//		x0 = (polygon[0].x + polygon[i].x + polygon[i + 1].x) / 3;
//		y0 = (polygon[0].y + polygon[i].y + polygon[i + 1].y) / 3; //求当前三角形的重心
//		cs = ISL_multiply(polygon[i], polygon[i + 1], polygon[0]) / 2;

//		//三角形面积可以直接利用该公式求解
//		if (fabs(s) < 1e-20) {
//			x = x0;
//			y = y0;
//			s += cs;
//			continue;
//		}

//		k = cs / s; //求面积比例
//		x = (x + k * x0) / (1 + k);
//		y = (y + k * y0) / (1 + k);
//		s += cs;
//	}
//	tp.x = x;
//	tp.y = y;

//	return tp;
//}

//} /* End Of namespace ISLIB
