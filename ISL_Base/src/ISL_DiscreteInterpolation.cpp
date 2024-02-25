//*
//*	@file	ISL_DiscreteInterpolation.cpp
//*	@brief	[Source file of Discrete Interpolation Functions], 离散插值；
//*	@see	ISeisLib Manual
//*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
//*	@date	2014-07-02
//*	@refer
//*/


//#include "ISL_DiscreteInterpolation.h"
//#include "ISL_Geometric2D.h"
//#include "ISL_MathBase.h"
//#include "ISL_Numerical.h"
//#include "ISL_Fitting.h"


//namespace ISLib {

///* 计算任意数据点（Xi, Yi, Zi）到 离插值点A的平面投影距离, d为计算的距离  */
//void ISL_IDW_distance2p(float xi, float yi, float xa, float ya, float &d)
//{
//	POINT i, a;
//	i.x = xi;
//	i.y = yi;
//	a.x = xa;
//	a.y = ya;

//	d = ISL_p2pDist(i, a);
//}


///* 将两点距离最近的一系列坐标数组的下标返回,
// *  nout = 给定的返回数组个数；
// *  out[] = 返回的数组，记录输入数组的下标
// *  注意：nout<=n*/
//void ISL_IDW_nearpoints(float * xi, float * yi, int n,
//						float xa, float ya,
//						int nout, int *& out)
//{
//	struct distance
//	{
//		float d;
//		int idx;
//	};

//	float di = 0;
//	distance * d = new distance[n];

//	for (int i = 0; i < n; i++) {
//		ISL_IDW_distance2p(xi[i], yi[i], xa, ya, di);
//		d[i].d = di;
//		d[i].idx = i;
//	}

//	for (int m = n - 1; m > 0; m--) {
//		distance temp;
//		for (int i = 0; i < m; i++) {
//			if (d[i].d > d[i + 1].d) {
//				temp = d[i];
//				d[i] = d[i + 1];
//				d[i + 1] = temp;
//			}
//		}
//	}

//	for (int i = 0; i < nout; i++) {
//		out[i] = d[i].idx;
//	}
//	if(d){
//		delete []d;
//		d = NULL;
//	}
//}


///* 计算与多点的内插值  */
//void ISL_IDW_intz(float *xi, float *yi, float *zi, int n,
//				float xa, float ya, float &za)
//{
//	float up = 0, down = 0, d = 0;

//	for(int i=0; i<n; i++){
//		ISL_IDW_distance2p(xi[i], yi[i], xa, ya, d);
//		up = up + zi[i]/d;
//		down = down + 1/d;
//	}
//	za = up / down;
//}

///* 利用筛选的多个点与插值点计算内插值，与ISL_IDW_nearpoints配合使用  */
//void ISL_IDW_intz(float *xi, float *yi, float *zi, int n,
//					int *id, int id_num,
//					float xa, float ya, float &za)
//{
//	if(n<0) return;
//	float up = 0, down = 0, d = 0;

//	for(int i=0; i<id_num; i++){
//		ISL_IDW_distance2p(xi[id[i]], yi[id[i]], xa, ya, d);
//		if(ISL_floatEqual(d, (float)0) == 1){
//			za = zi[id[i]];
//			return;
//		}
//		up = up + zi[id[i]]/d;
//		down = down + 1/d;
//	}
//	za = up / down;

//}


///**********************************************************************************************************************
// *
// *  功能：二维反距离加权插值法
// *
// *  说明:	near_num <= nin
// *  	nin < noutx*nouty
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		int					nin					In			给定样点的个数
// *		float				xi[]				In			给定样点的X坐标	 array[nin]
// *		float				yi[]				In			给定样点的Y坐标	 array[nin]
// *		float				zi[]				In			给定样点的X,Y坐标位置的样点值 array[nin]
// *
// *		int					near_num			In			就近计算的样点个数（near_num <= nin）
// *
// *		int					noutx				In			输出数据X轴的个数
// *		int					nouty				In			输出数据y轴的个数（noutx*nouty > nin）
// *		float				xout[]				In			输出数据的X坐标  array[noutx]
// *		float				yout[]				In			输出数据的Y坐标  array[nouty]
// *
// *		float				yout[]				Out			输出数据的样点值  array[noutx*nouty]
// *
// *  返回：无
// *
//**********************************************************************************************************************/
//void ISL_dintIDW2D(	int nin, float *xi, float *yi, float *zi, int near_num,
//					int noutx, int nouty, float *xout, float *yout, float *&zout)
//{
//	int * near_idx = new int[near_num];

//	for (int i = 0; i < noutx; i++) {

//		for (int j = 0; j < nouty; j++) {
//			/* 获得最近的已知值 */
//			ISL_IDW_nearpoints(xi, yi, nin, xout[i], yout[j], near_num, near_idx);

//			/* 计算插值位置的内插值 */
//			ISL_IDW_intz(xi, yi, zi, nin, near_idx, near_num, xout[i], yout[j], zout[i * nouty + j]);
//		}
//	}
//	if(near_idx){
//		delete []near_idx;
//		near_idx = NULL;
//	}
//}

//// ==============================================================================================================
//// ==============================================================================================================
//// ==============================================================================================================

///************************************************************
//	range 所有特征点的Xmin、 Ymin、  Xmax、 Ymax
//	int mode， 计算半方差矩阵的模式 1 2 3三种情况
//	int item， 特征点数量
//	float *Z_s， item 个点的采样值
//	float *pos : X1 Y1  X2  Y2  …  Xitem  Yitem 所有特征点的坐标

//	int *resol,x方向和y方向的分辨率。
//		即，	resol_x为x方向的间隔数；
//			resol_y为y方向的间隔数，
//			最后总点数为：（resol_x+1）*（resol_y+1）


//	=== 三个常量参数 ===
//	块金值(Nugget), float c0
//	基抬值(sill), float c1
//	变程(a), float a

//	result 输出结果:（resol_x+1）*（resol_y+1）个点的插值结果
//	x带表行位置，y代表列位置
//***************************************************************/
//void ISL_krigingCore(int * range/* array[4] */, int mode,		/* 数据范围， 计算方式 */
//					int item, float * Z_s, float * pos,			/* 输入的随机样点数 */
//					float c0, float c1, float a,				/* 计算常量 */
//					int * resol/* array[2] */, float * result)	/* 输出数据 */
//{
//	int dim, i, j, k, l;
//	float i_f, j_f, cnt_x, cnt_y;
//	float begin_row, begin_col, end_row, end_col, *Cd, test_t;
//	int resol_x, resol_y;
//	double *D = NULL, *V = NULL, *temp = NULL;

////	FILE *fp; //定义文件用于存放结果
////	fp = fopen("test_out.txt", "w");//打开文件

//	/* initialize values */
//	//得到范围坐标
//	begin_row = (float) *range;
//	begin_col = (float) *(range + 1);
//	end_row = (float) *(range + 2);
//	end_col = (float) *(range + 3);

//	//为特征点数+1
//	dim = item + 1;

//	//得到分辩率
//	resol_x = *(resol);
//	resol_y = *(resol + 1);

//	/* allocate V D array */
//	//V为(n+1)*(n+1)矩阵，V*W=D
//	//D为(n+1)*1矩阵，V*W=D
//	//temp临时存放V，大小与V同
//	//W放在D中，故不为W分配空间
//	V = (double *) malloc(sizeof(double) * dim * dim);
//	D = (double *) malloc(sizeof(double) * dim);
//	temp = (double *) malloc(sizeof(double) * dim * dim);
//	/* allocate Cd array */
//	// Cd为(n+1)*(n+1)矩阵，用于存储距离
//	// n 即item
//	//     D0 0        D0 1 …      D0 item-1         1
//	//     D1 0        D1 1 …      D1 item-1         1
//	//     Ditem-1 0   Ditem-1 1 … Ditem-1 item-1    1
//	//     1           1            1                0
//	Cd = (float *) malloc(sizeof(float) * dim * dim);
//	//计算上面的距离，用的欧氏方法
//	/* caculate the distance between sample datas put into Cd array*/
//	for (i = 0; i < dim - 1; i++)
//		for (j = i; j < dim - 1; j++) {
//			test_t = (pos[i * 2] - pos[j * 2]) * (pos[i * 2] - pos[j * 2])
//					+ (pos[i * 2 + 1] - pos[j * 2 + 1]) * (pos[i * 2 + 1] - pos[j * 2 + 1]);
//			Cd[i * dim + j] = (float) sqrt(test_t);
//		}
//	//补1
//	for (i = 0; i < dim - 1; i++) {
//		V[i * dim + dim - 1] = 1;
//		V[(dim - 1) * (dim) + i] = 1;
//	}
//	//添0
//	V[(dim - 1) * (dim) + i] = 0;

//	/* caculate the variogram of sample datas and put into  V array */
//	//依据Dij计算采样点的半方差矩阵并放入V
//	for (i = 0; i < dim - 1; i++){
//		for (j = i; j < dim - 1; j++) {	//由于对称，只计算一半
//			switch (mode)
//			{	//模式不同，公式不同
//				case 1: /* Spher mode */
//					if (Cd[i * dim + j] < a)
//						V[i * dim + j] = V[j * dim + i] = (float) (c0
//								+ c1 * (1.5 * Cd[i * dim + j] / a - 0.5 * (Cd[i * dim + j] / a)
//																	* (Cd[i * dim + j] / a)
//																	* (Cd[i * dim + j] / a)));
//					else
//						V[i * dim + j] = V[j * dim + i] = c0 + c1;
//				break;
//				case 2: /* Expon mode */
//					V[i * dim + j] = V[j * dim + i] = (float) (c0
//							+ c1 * (1 - exp(-3 * Cd[i * dim + j] / a)));
//				break;
//				case 3: /* Gauss mode */
//					V[i * dim + j] = V[j * dim + i] = (float) (c0
//							+ c1 * (1 - exp(-3 * Cd[i * dim + j] * Cd[i * dim + j] / a / a)));
//				break;
//				default:
//					V[i * dim + j] = V[j * dim + i] = (float) (1 * Cd[i * dim + j]);
//				break;
//			}
//		}
//	}

//	/* release Cd array */
//	free(Cd);
//	//Cd完成任务，释放

//	//先暂存V,因GAUSS算法要改V
//	for (k = 0; k < dim * dim; k++)
//		temp[k] = V[k];

//	//计算所有范围内，以分辨率为间隔的所有点
//	//例：range:
//	//2,4,   10,20
//	//resol:4,4 两方向均取4格，
//	//则待算行列方向步长为：cnt_x=（10-2）/4=2    cnt_y=(20-4)/4=4
//	//待算点为：
//	//2  4    2  8    2  12    2  16    2  20
//	//4  4    4  8    4  12    4  16    4  20
//	//6  4    6  8    6  12    6  16    6  20
//	//8  4    8  8    8  12    8  16    8  20
//	//10 4    10 8    10 12    10 16    10 20
//	/* for loop for each point of the estimated block */
//	cnt_x = (end_row - begin_row) / (float) resol_x;	//x方向步长
//	cnt_y = (end_col - begin_col) / (float) resol_y;	// y方向步长

//	l = 0;
//	for (i = 0; i <= resol_x; i++) {
////		cout<<"i = "<<i<<endl;
//		i_f = cnt_x * i + begin_row; //2 4 6  8  10
//		for (j = 0; j <= resol_y; j++) {
//			j_f = cnt_y * j + begin_col; //4 8 12 16 20

//			//           得到待计算点坐标：	行：i_f
//			//           				列：j_f
//			//           针对该点计算D的item个值,并补一个1
//			for (k = 0; k < dim - 1; k++) {
//				test_t = (i_f - pos[2 * k]) * (i_f - pos[2 * k])
//						+ (j_f - pos[2 * k + 1]) * (j_f - pos[2 * k + 1]);
//				test_t = (float) (sqrt(test_t));
//				switch (mode)
//				{			//模式不同，公式不同
//					case 1: /* Spher mode */
//						if (test_t < a)
//							D[k] = (float) (c0
//									+ c1 * (1.5 * test_t / a - 0.5 * (test_t / a) * (test_t / a)
//																* (test_t / a)));
//						else
//							D[k] = c0 + c1;
//					break;
//					case 2: /* Expon mode */
//						D[k] = (float) (c0 + c1 * (1 - exp(-3 * test_t / a)));
//					break;
//					case 3: /* Gauss mode */
//						D[k] = (float) (c0 + c1 * (1 - exp(-3 * test_t * test_t / a / a)));
//					break;
//					default:
//						D[k] = (float) (1 * test_t);
//					break;
//				}

//			}
//			D[dim - 1] = 1;			//并补一个1
//			//现在V有了,D有了,可解方程了

//			//先取回V
//			for (k = 0; k < dim * dim; k++)
//				V[k] = temp[k];

//			ISL_gaus(V, D, dim);
//			//解出的结果:W1,W2,...,Witem,λ在D中，替换D原来的值

//			//下面计算该点的预测值,依据公式:
//			//W1Y1+W2Y2+...+WitemYitem有:

//			test_t = 0;			//用于累加
//			for (k = 0; k < dim - 1; k++)
//				test_t += D[k] * Z_s[k];

//			*(result + l++) = test_t;
//			//得到该点的最终结果

//			//////////////////////////
//			//(int)(i_f+0.5)     坐标 x
//			//(int)(j_f+0.5)     坐标 y
//			//(int)(test_t+0.5)  值   z
//			//写入文件
////			fprintf(fp, "%d,%d,%d\n", (int) (i_f + 0.5), (int) (j_f + 0.5), (int) (test_t + 0.5));
//			//////////////////////////
//		}
//	}
//	//关闭文件
////	fclose(fp);
//	free(V);
//	free(D);
//	free(temp);
//}

///**********************************************************************************************************************
// *
// *  功能：二维普通克里金插值法
// *
// *  说明:
// *
// *  参数：
// *		Type				Name				In/Out		Description
// *		----				----				------		-----------
// *		int	*				range				In			插值区域范围(矩形)，数组大小为4，array[4]
// *															range[0] = 左上角 X
// *															range[1] = 左上角 Y
// *															range[2] = 右下角 X
// *															range[3] = 右下角 Y
// *
// *		int					mode				In			 计算半方差矩阵的 3三种模式
// *																1 = 球状模型
// *																2 = 指数模型
// *																3 = 高斯模型
// *
// *		int					nin					In			给定的样点个数
// *		float				xi[]				In			给定样点的X坐标	 array[nin]
// *		float				yi[]				In			给定样点的Y坐标	 array[nin]
// *		float				zi[]				In			给定样点的X,Y坐标位置的样点值 array[nin]
// *
// *		float				c0					In			块金值(Nugget)
// *		float				c1					In			基抬值(sill)
// *		float				a					In			变程值(a)
// *
// *		int					noutx				In			输出数据X轴的个数
// *		int					nouty				In			输出数据y轴的个数（noutx*nouty > nin）
// *
// *		float				yout[]				Out			输出数据的样点值  array[noutx*nouty]
// *
// *  返回：无
// *
//**********************************************************************************************************************/
//void ISL_dintKriging2D (	int *range, int mode,
//						int nin, float *xi, float *yi, float *zi, /* 输入随即样点 */
//						int noutx, int nouty, float *&zout/* array[noutx * nouty] */,
//						float c0, float c1, float a)
//{
//	float * pos = new float[2*nin];
//	for(int i=0; i<nin; i++){
//		pos[2*i] = xi[i];
//		pos[2*i+1] = yi[i];
//	}

//	int resol[2];
//	resol[0] = noutx - 1;
//	resol[1] = nouty - 1;

//	ISL_krigingCore(range, mode,
//					nin, zi, pos,
//					c0, c1, a,
//					resol, zout);

//	if(pos){ delete []pos; pos = NULL; }
//}


///** @brief 一般趋势拟合值Q，求多项式p，最少6个样点坐标
// * 	@param[in]	* x		给定样点的X坐标；
// * 	@param[in]	* y		给定样点的Y坐标；
// * 	@param[in]	* z		给定样点的样点值；
// * 	@param[in]	n		给定样点的个数；
// * 	@return 一般趋势拟合值Q（double）
// */
//double ISL_QC_LSDWI(float * x, float * y, float * z, int n)
//{
//	float cof[6];
//	float * simu = new float[n];

//	ISL_lscf(x, y, n, cof, 6, simu);

//	double QC = 0;
//	for(int i=0; i<n; i++){
//		double p = cof[0] + cof[1]*x[i] + cof[2]*y[i] + cof[3]*x[i]*x[i] + cof[4]*x[i]*y[i] + cof[5]*y[i]*y[i];
//		QC = QC + (p - z[i])*(p - z[i]);
//	}

//	if(simu){
//		delete []simu;
//		simu = NULL;
//	}

//	return QC;
//}

///** @brief 加权趋势拟合值Q，最少6个样点坐标
// * 	@param[in]	* x		给定样点的X坐标；
// * 	@param[in]	* y		给定样点的Y坐标；
// * 	@param[in]	* z		给定样点的样点值；
// * 	@param[in]	n		给定样点的个数；
// * 	@param[in]	a		网格点X坐标；
// * 	@param[in]	b		网格点Y坐标；
// * 	@return 加权趋势拟合值Q（double）
// */
//double ISL_QW_LSDWI(float * x, float * y, float * z, int n, float a, float b)
//{
//	if(x == NULL || y == NULL || z== NULL || n<6)
//			return NullValue;

//	double QW = 0;
//	for(int i=0; i<n; i++){
//		double d = 1.0 / pow( ( (x[i]-a)*(x[i]-a) + (y[i]-b)*(y[i]-b) + EP ), 2 );
//		QW = QW + d;
//	}

//	double QC = ISL_QC_LSDWI(x, y, z, n);
//	QW = QW * QC;

//	return QW;
//}


//}/* End Of namespace ISLIB
