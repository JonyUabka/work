///**
//*	@file	ISL_Transform.cpp
//*	@brief	[Source file of Transform Functions], 数值变换；
//*	@see	ISeisLib Manual
//*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang, Chen Ke], 刘百红、杨强、宋志翔、陈科；
//*	@date	2014-06-03
//*	@refer	SU CWP
//*/

//#include "ISL_Transform.h"
//#include "ISL_Alloc.h"

//namespace ISLib {

///* 希尔伯特变换 */
//void ISL_hilbert(int num, float * in, float * out, hilbert_type type)
//{
//	if(type == NEWS_hlb){
//		float * outRe = new float[num];
//		ISL_hilbert(num, in, outRe, out);
//		if(outRe){
//			delete []outRe;
//			outRe = NULL;
//		}
//	}
//	if(type == SU_hlb){
//		static int madeh = 0;
//		static float h[LH];
//		int i;
//		float taper;

//		/* if not made, make Hilbert transform filter; use Hamming window */
//		if (!madeh) {
//			h[LHHALF] = 0.0;
//			for (i = 1; i <= LHHALF; i++) {
//				taper = 0.54 + 0.46 * cos(PI * (float) i / (float) (LHHALF));
//				h[LHHALF + i] = taper * (-(float) (i % 2) * 2.0
//						/ (PI * (float) (i)));
//				h[LHHALF - i] = -h[LHHALF + i];
//			}
//			madeh = 1;
//		}

//		/* convolve Hilbert transform with input array */
//		ISL_conv(LH, -LHHALF, h, num, 0, in, num, 0, out);
//	}
//}

///* 希尔伯特变换, 返回实部与虚部 */
//void ISL_hilbert(int num, float *in, float *outRe, float *outIm)
//{
//	int nNumber = 1, nexp, half;

//	for (int i = 1;; i++) {
//		nNumber *= 2;
//		nexp = i;
//		if (nNumber >= num)
//			break;
//	}

//	float *dataRe = NULL, *dataIm = NULL, *fftRe = NULL, *fftIm = NULL;
//	half = nNumber / 2;
//	dataRe = new float[nNumber];
//	dataIm = new float[nNumber];
//	fftRe = new float[nNumber];
//	fftIm = new float[nNumber];

//	memset(dataRe, 0, sizeof(float) * nNumber);
//	memset(dataIm, 0, sizeof(float) * nNumber);
//	memset(fftRe, 0, sizeof(float) * nNumber);
//	memset(fftIm, 0, sizeof(float) * nNumber);
//	memcpy(dataRe, in, sizeof(float) * num);
//	ISL_kfft(dataRe, dataIm, nNumber, nexp, fftRe, fftIm, 0, 0);

//	dataRe[0] = fftRe[0];
//	dataIm[0] = fftIm[0];
//	for (int i = 1; i < nNumber; i++) {
//		if (i <= half) {
//			dataRe[i] = 2 * fftRe[i];
//			dataIm[i] = 2 * fftIm[i];
//		} else {
//			dataRe[i] = 0.;
//			dataIm[i] = 0.;
//		}
//	}

//	ISL_kfft(dataRe, dataIm, nNumber, nexp, fftRe, fftIm, 1, 0);
//	for (int i = 0; i < num; i++){
//		outRe[i] = fftRe[i];
//		outIm[i] = fftIm[i];
//	}

//	delete[] dataRe;
//	dataRe = NULL;
//	delete[] dataIm;
//	dataIm = NULL;
//	delete[] fftRe;
//	fftRe = NULL;
//	delete[] fftIm;
//	fftIm = NULL;
//}

///* 沃什（Walsh）变换 */
//void ISL_kfwt(double *p, int n, int k, double *x)
//{
//	int m, l, it, ii, i, j, is;
//	double q;
//	m = 1;
//	l = n;
//	it = 2;
//	x[0] = 1;
//	ii = n / 2;
//	x[ii] = 2;

//	for (i = 1; i <= k - 1; i++) {
//		m = m + m;
//		l = l / 2;
//		it = it + it;
//		for (j = 0; j <= m - 1; j++)
//			x[j * l + l / 2] = it + 1 - x[j * l];
//	}

//	for (i = 0; i <= n - 1; i++) {
//		ii = x[i] - 1;
//		x[i] = p[ii];
//	}

//	l = 1;
//	for (i = 1; i <= k; i++) {
//		m = n / (2 * l) - 1;
//		for (j = 0; j <= m; j++) {
//			it = 2 * l * j;
//			for (is = 0; is <= l - 1; is++) {
//				q = x[it + is] + x[it + is + l];
//				x[it + is + l] = x[it + is] - x[it + is + l];
//				x[it + is] = q;
//			}
//		}
//		l = 2 * l;
//	}
//	return;
//}


///* 单频小波变换，输入一个信号，输出信号的小波变换的单频谱 */
//void ISL_sfCWT(float *in, int num, float *wavelet, int nw, float *ampSpec, float *phSpec)
//{
//	float *seismic = new float[num];
//	float *complexIm = new float[num];

//	int counter = 0;
//	memset(ampSpec, 0, sizeof(float) * num);
//	memset(phSpec, 0, sizeof(float) * num);

//	ISL_synConv(in, num, wavelet, nw, seismic);
//	ISL_hilbert(num, seismic, complexIm);

//	for (int j = 0; j < num; j++) {
//		ampSpec[j] += sqrt(seismic[j] * seismic[j] + complexIm[j] * complexIm[j]);
//		if (seismic[j] * seismic[j] + complexIm[j] * complexIm[j])
//			phSpec[j] += atan2(complexIm[j], seismic[j]);
//		else
//			phSpec[j] += 0.;
//	}
//	counter++;

//	for (int j = 0; j < num; j++) {
//		ampSpec[j] /= counter;
//		phSpec[j] /= counter;
//	}

//	delete[] seismic;
//	delete[] complexIm;
//}


///*在Z平面单位圆上计算有限长序列x(n)的Z变换的采样值*/
//void ISL_czt(double *xr, double *xi, int n, int m, double f1, double f2)
//{
//	int i, j, n1, n2, len;
//	double e, t, ar, ai, ph, pi, tr, ti, *wr = NULL, *wr1 = NULL, *wi = NULL, *wi1 = NULL;
//	len = n + m - 1;
//	for (j = 1, i = 1; i < 16; i++) {
//		j = j * 2;
//		if (j >= len) {
//			len = j;
//			break;
//		}
//	}

//	wr 	= alloc1double( len );
//	wi 	= alloc1double( len );
//	wr1 = alloc1double( len );
//	wi1 = alloc1double( len );
//	pi 	= 3.14159265358979;

//	ph = 2.0 * pi * (f2 - f1) / (m - 1);

//	n1 = (n >= m) ? n : m;
//	for (i = 0; i < n1; i++) {
//		e = ph * i * i / 2.0;
//		wr[i] = cos(e);
//		wi[i] = sin(e);
//		wr1[i] = wr[i];
//		wi1[i] = -wi[i];
//	}

//	n2 = len - n + 1;
//	for (i = m; i < n2; i++) {
//		wr[i] = 0.0;
//		wi[i] = 0.0;
//	}
//	for (i = n2; i < len; i++) {
//		j = len - i;
//		wr[i] = wr[j];
//		wi[i] = wi[j];
//	}

//	ISL_kfft_d(wr, wi, len, 1);
//	ph = -2.0 * pi * f1;
//	for (i = 0; i < n; i++) {
//		e = ph * i;
//		ar = cos(e);
//		ai = sin(e);
//		tr = ar * wr[i] - ai * wi[i];
//		ti = ai * wr1[i] + ar * wi1[i];
//		t = xr[i] * tr - xi[i] * ti;
//		xi[i] = xr[i] * ti + xi[j] * tr;
//		xr[i] = t;
//	}
//	for (i = n; i < len; i++) {
//		xr[i] = 0.0;
//		xi[i] = 0.0;
//	}

//	ISL_kfft_d(xr, xi, len, 1);
//	for (i = 0; i < len; i++) {
//		tr = xr[i] * wr[i] - xi[i] * wi[i];
//		xi[i] = xr[i] * wi[i] + xi[i] * wr[i];
//		xr[i] = tr;
//	}

//	ISL_kfft_d(xr, xi, len, -1);
//	for (i = 0; i < m; i++) {
//		tr = xr[i] * wr1[i] - xi[i] * wi1[i];
//		xi[i] = xr[i] * wi1[i] + xi[i] * wr1[i];
//		xr[i] = tr;
//	}
//	free1double(wr);
//	free1double(wi);
//	free1double(wr1);
//	free1double(wi1);
//}

//} /* End Of namespace ISLIB */

