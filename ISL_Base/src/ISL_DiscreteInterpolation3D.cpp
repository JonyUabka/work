///**
//*	@file	ISL_DiscreteInterpolation3D.cpp
//*	@brief	[Source file of Discrete Interpolation Functions], 离散插值3D；
//*	@see	ISeisLib Manual
//*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
//*	@date	2014-07-02
//*	@refer
//*/


//#include "ISL_DiscreteInterpolation.h"
//#include "ISL_Geometric2D.h"
//#include "ISL_Geometric3D.h"
//#include "ISL_MathBase.h"
//#include "ISL_Numerical.h"


//namespace ISLib {

///* 计算任意数据点（Xi, Yi, Zi）到 离插值点A的平面投影距离, d为计算的距离  */
//void ISL_IDW_distance2p3D(float xi, float yi, float zi, float xa, float ya, float za, float &d)
//{
//	POINT3D i, a;
//	i.x = xi;
//	i.y = yi;
//	i.z = zi;

//	a.x = xa;
//	a.y = ya;
//	a.z = za;

//	d = ISL_p2pDist3D(i, a);
//}


///* 将两点距离最近的一系列坐标数组的下标返回,
// *  nout = 给定的返回数组个数；
// *  out[] = 返回的数组，记录输入数组的下标
// *  注意：nout<=n*/
//void ISL_IDW_nearpoints3D(	float * xi, float * yi, float * zi, int n, float radius,
//							float xa, float ya, float za,
//							int nout, int *& out, int mode)
//{
//	struct distance
//	{
//		float d;
//		int idx;
//	};

//	float di = 0;
//	distance * d = NULL;


//	if(mode == 1){
//		/** 快速插值 */

//		d = new distance[n];
//		int iCount = 0;
//		for (int i = 0; i < n; i++) {
//			ISL_IDW_distance2p3D(xi[i], yi[i], zi[i], xa, ya, za, di);

//			d[i].d = di;
//			d[i].idx = i;

//			if(di < radius){
//				out[iCount] = d[i].idx;
//				iCount++;
//			}

//			if(iCount>=nout)
//				break;
//		}

//		if(iCount<nout){
//			for (int i = iCount; i < nout; i++) {
//				out[i] = i-iCount;
//			}
//		}


//	}else{
//		/** 常规插值 */

//		d = new distance[n];
//		for (int i = 0; i < n; i++) {
//			ISL_IDW_distance2p3D(xi[i], yi[i], zi[i], xa, ya, za, di);
//			d[i].d = di;
//			d[i].idx = i;
//		}

//		for (int m = n - 1; m > 0; m--) {
//			distance temp;
//			for (int i = 0; i < m; i++) {
//				if (d[i].d > d[i + 1].d) {
//					temp = d[i];
//					d[i] = d[i + 1];
//					d[i + 1] = temp;
//				}
//			}
//		}

//		for (int i = 0; i < nout; i++) {
//			out[i] = d[i].idx;
//		}
//	}

//	if(d){
//		delete []d;
//		d = NULL;
//	}
//}


///** 利用筛选的多个点与插值点计算内插值，与ISL_IDW_nearpoints3D配合使用  */
//void ISL_IDW_intz3D(float *xi, float *yi, float *zi, float *si, int n,
//					int *id, int id_num,
//					float xa, float ya, float za, float &sa)
//{
//	if(n<0) return;
//	float up = 0, down = 0, d = 0;

//	for(int i=0; i<id_num; i++){
//		ISL_IDW_distance2p3D(xi[id[i]], yi[id[i]], zi[id[i]], xa, ya, za, d);
//		if(ISL_floatEqual(d, (float)0) == 1){
//			sa = si[id[i]];
//			return;
//		}
//		up = up + si[id[i]]/d;
//		down = down + 1/d;
//	}
//	sa = up / down;
//}


///**********************************************************************************************************************
// *
// *  功能：三维反距离加权插值法
// *
// *  说明:	near_num <= nin
// *  	nin < noutx*nouty*noutz
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		int					nin					In			给定样点的个数
// *		float				xi[]				In			给定样点的X坐标	 array[nin]
// *		float				yi[]				In			给定样点的Y坐标	 array[nin]
// *		float				zi[]				In			给定样点的Z坐标	 array[nin]
// *		float				si[]				In			给定样点的X,Y,Z坐标位置的样点值 array[nin]
// *
// *		float				radius				In			搜索半径，小于noutx，nouty，noutz中的最小值 ；
// *		int					near_num			In			就近计算的样点个数（near_num <= nin）
// *
// *		int					noutx				In			输出数据X轴的个数
// *		int					nouty				In			输出数据y轴的个数（noutx*nouty*noutz > nin）
// *		int					noutz				In			输出数据z轴的个数（noutx*nouty*noutz > nin）
// *		float				xout[]				In			输出数据的X坐标  array[noutx]
// *		float				yout[]				In			输出数据的Y坐标  array[nouty]
// *		float				zout[]				In			输出数据的Y坐标  array[noutz]
// *
// *		float				sout[]				Out			输出数据的样点值  array[noutx*nouty*noutz]
// *		int					mode				In			插值方式，0=常规，1=快速；
// *
// *  返回：无
// *
//**********************************************************************************************************************/
//void ISL_dintIDW3D(	int nin, float *xi, float *yi, float *zi, float *si, float radius, int near_num,
//					int noutx, int nouty, int noutz, float *xout, float *yout, float *zout, float *&sout, int mode)
//{
//	if(nin < near_num){
//		printf("input failed\n");
//		return;
//	}

//	if( nin <= 5){
//		near_num = 5;
//	}

//	int * near_idx = new int[near_num];

//	for (int i = 0; i < noutx; i++) {
//		cout<<"i = "<<i<<endl;

//		for (int j = 0; j < nouty; j++) {

//			for (int k = 0; k < noutz; k++) {

//				/* 获得最近的已知值 */
//				ISL_IDW_nearpoints3D(xi, yi, zi, nin, radius, xout[i], yout[j], zout[k], near_num, near_idx, mode);

//				/* 计算插值位置的内插值 */
//				ISL_IDW_intz3D(xi, yi, zi, si, nin, near_idx, near_num, xout[i], yout[j], zout[k],
//								sout[i * nouty * noutz + j * noutz + k]);
//			}
//		}
//	}
//	if(near_idx){
//		delete []near_idx;
//		near_idx = NULL;
//	}
//}

//}/* End Of namespace ISLIB */




