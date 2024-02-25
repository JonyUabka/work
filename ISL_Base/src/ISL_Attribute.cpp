/**
*	@file	ISL_Attribute.cpp
*	@brief	[Source file of Attribute Functions], 地震属性；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-06-26
*	@refer
*/

#include "ISL_Attribute.h"
#include "ISL_Transform.h"
#include "ISL_Absorption.h"


namespace ISLib {

/**
* @brief	计算三瞬属性, 三瞬属性计算，包括瞬时振幅、瞬时相位和瞬时频率
*
* @param[in]	mode		输入的计算类型，0表示计算瞬时振幅，1表示计算瞬时相位，2表示计算瞬时频率
* @param[in]	in			输入的数组
* @param[in]	num			输入的数组的长度
* @param[in]	delta		输入的地震数据的采样间隔，以秒为单位
* @param[out]	out			输出经希尔伯特变换得到的复信号的虚部(外部需分配好内存空间)
*
* @return	无
*/
void ISL_insAttributes(int mode, float *in, int num, float delta, float *out)
{
	float * complexIm = new float[num];

	ISL_hilbert(num, in, complexIm);
	/////计算瞬时振幅/////
	if (mode == 0) {
		for (int i = 0; i < num; i++)
			out[i] = sqrt(in[i] * in[i] + complexIm[i] * complexIm[i]);
	}

	/////计算瞬时相位/////
	if (mode == 1) {
		for (int i = 0; i < num; i++)
			ISL_calcPh(in[i], complexIm[i], out[i]);
	}

	/////计算瞬时频率/////
	if (mode == 2) {
		float Tmp;
		float phase = 0;
		for (int i = 0; i < num; i++) {
			ISL_calcPh(in[i], complexIm[i], Tmp);
			ISL_calcFreq(phase, Tmp, delta, out[i]);
			phase = Tmp;
		}
	}

	if(complexIm){
		delete[] complexIm;
		complexIm = NULL;
	}
}



/**********************************************************************************************************************
 *
 *  功能：快速傅里叶变换(相位转换专用)
 *
 *  说明：本函数用于计算旋转地震数据的相位
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		float*				dateRe				In/Out		当l=0时，存放nPoints个采样的实部，返回离散傅里叶变换的模；
 *															当l=1时，存放傅里叶变换的nPoints个实部，返回傅里叶变换的模
 *		float*				dateIm				In/Out		当l=0时，存放nPoints个采样的虚部，返回离散傅里叶变换的幅角；
 *															当l=1时，存放傅里叶变换的nPoints个虚部，返回傅里叶变换的幅角。
 *															其中幅角的单位为度
 *		int					nPoints				In			输入数据的采样点数
 *		int					nExp				In			满足nPoints=2^nExp   由于nPoints和nExp存在这样的关系，应用时受到一定的限制
 *		float*				fftRe				In/Out		当l=0时，返回傅里叶变换的实部；
 *															当l=1时，返回傅里叶变换的实部
 *		float*				fftIm				In/Out		当l=0时，返回傅里叶变换的虚部；
 *															当l=1时，返回傅里叶变换的虚部
 *		int					l					In			当l=0时，表示要求本函数计算傅里叶变换；
 *															当l=1时，表示要求本函数计算逆傅里叶变换
 *		int					il					In			当il=0时，表示不要求本函数计算傅里叶变换或逆变换的模与幅角；
 *															当il=1时，表示要求本函数计算傅里叶变换或逆变换的模与幅角
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void AlterFFT(float *pr, float *pi, int n, int k, float *fr, float *fi, int l, int il)
{
	int it, m, is, i, j, nv, l0;
	float ave = 0.;
	float p, q, s, vr, vi, poddr, poddi;

	for (it = 0; it <= n - 1; it++) {
		m = it;
		is = 0;
		for (i = 0; i <= k - 1; i++) {
			j = m / 2;
			is = 2 * is + (m - 2 * j);
			m = j;
		}
		fr[it] = pr[is];
		fi[it] = pi[is];
	}

	pr[0] = 1.0;
	pi[0] = 0.0;
	p = float(6.28318530717959 / (1.0 * n));

	pr[1] = float(cos(p));
	pi[1] = float(-sin(p));
	if (l != 0)
		pi[1] = -pi[1];

	for (i = 2; i <= n - 1; i++) {
		p = pr[i - 1] * pr[1];
		q = pi[i - 1] * pi[1];
		s = (pr[i - 1] + pi[i - 1]) * (pr[1] + pi[1]);
		pr[i] = p - q;
		pi[i] = s - p - q;
	}

	for (it = 0; it <= n - 2; it = it + 2) {
		vr = fr[it];
		vi = fi[it];
		fr[it] = vr + fr[it + 1];
		fi[it] = vi + fi[it + 1];
		fr[it + 1] = vr - fr[it + 1];
		fi[it + 1] = vi - fi[it + 1];
	}

	m = n / 2;
	nv = 2;

	for (l0 = k - 2; l0 >= 0; l0--) {
		m = m / 2;
		nv = 2 * nv;
		for (it = 0; it <= (m - 1) * nv; it = it + nv) {
			for (j = 0; j <= (nv / 2) - 1; j++) {
				p = pr[m * j] * fr[it + j + nv / 2];
				q = pi[m * j] * fi[it + j + nv / 2];
				s = pr[m * j] + pi[m * j];
				s = s * (fr[it + j + nv / 2] + fi[it + j + nv / 2]);
				poddr = p - q;
				poddi = s - p - q;
				fr[it + j + nv / 2] = fr[it + j] - poddr;
				fi[it + j + nv / 2] = fi[it + j] - poddi;
				fr[it + j] = fr[it + j] + poddr;
				fi[it + j] = fi[it + j] + poddi;
			}
		}
	}

	if (l != 0) {
		for (i = 0; i <= n - 1; i++) {
			fr[i] = float(fr[i] / (1. * n));
			fi[i] = float(fi[i] / (1. * n));
		}
	}
	if (il != 0) {
		for (i = 0; i <= n - 1; i++)
			if (ave < fabs(fr[i]))
				ave = float(fabs(fr[i]));

		for (i = 0; i <= n - 1; i++) {
			pr[i] = sqrt(fr[i] * fr[i] + fi[i] * fi[i]);
			if (fr[i] < 0) {
				pr[i] *= -1;
				fr[i] *= -1;
				fi[i] *= -1;
			}
			if (fabs(fr[i]) <= 0.000001 * fabs(fi[i]) || fr[i] == 0.0) {
				if ((fi[i] * fr[i]) >= 0)
					pi[i] = 90.0;
				else
					pi[i] = -90.0;
			} else
				pi[i] = float(atan(fi[i] / fr[i]) * 360.0 / 6.28318530717959);
		}
	}
	return;
}


/* 相位旋转 */
void ISL_alterPh(float *in, int num, int angle, float *out)
{
	int i, N, n;
	float pai = 3.141592653589f, *pr = NULL, *pi = NULL, *fr = NULL, *fi = NULL;

	for (N = 1, i = 1;; i++) {
		N *= 2;
		n = i;
		if (N >= num)
			break;
	}

	pr = new float[N];
	pi = new float[N];
	fr = new float[N];
	fi = new float[N];

	memset(pr, 0, sizeof(float) * N);
	memset(pi, 0, sizeof(float) * N);
	memset(fi, 0, sizeof(float) * N);
	memset(fr, 0, sizeof(float) * N);
	memcpy(pr, in, sizeof(float) * num);

	AlterFFT(pr, pi, N, n, fr, fi, 0, 1);

	for (i = 0; i < N; i++) {
		if (i < N / 2) {
			fr[i] = pr[i] * cos((pi[i] + angle) * pai / 180);
			fi[i] = pr[i] * sin((pi[i] + angle) * pai / 180);
		} else {
			fr[i] = pr[i] * cos((pi[i] - angle) * pai / 180);
			fi[i] = pr[i] * sin((pi[i] - angle) * pai / 180);
		}
	}

	ISL_kfft(fr, fi, N, n, pr, pi, 1, 0);
	memcpy(out, pr, sizeof(float) * num);

	delete[] pr;
	pr = NULL;
	delete[] pi;
	pi = NULL;
	delete[] fr;
	fr = NULL;
	delete[] fi;
	fi = NULL;

	return;
}

} /* End Of namespace ISLIB */
