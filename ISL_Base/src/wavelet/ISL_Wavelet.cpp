/**
 *	@file	ISL_Wavelet.h
 *	@brief	[Header file of Wavelet], 子波；
 *	@see	ISeisLib Manual
 *	@author [Zhang Yang], 张杨；
 *	@date	2014-05-27
 *	@refer	SU CWP
 */

#include "ISL_Wavelet.h"
#include "ISL_Normalize.h"
#include "ISL_Absorption.h"

namespace ISLib {



/**********************************************************************************************************************
 *
 *  功能：雷克子波
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					nt					In			雷克子波的长度, number samples in output wavelet
 *		float				dt					In			时间采样间隔(秒), time step
 *		float				fpeak				In			输入的主频(HZ), peak frequency of the Ricker wavelet
 *		float				wavelet[nt]			Out			array[n] of computed wavelet
 *
 *  返回：float * ，返回数据的指针，得到后请自己负责释放内存空间,数组的大小为n
 *
 **********************************************************************************************************************/
void ISL_classicRicker(int nt, float dt, float fpeak, float *wavelet, int mode) {
	if (nt <= 0 || dt <= 0 || fpeak < 0)
		return;

	if (mode == 1) {
		float A = 1;
		float tmp = 0;

		for (int i = -nt/2; i<nt / 2; i++) {
			tmp = float(pow(((float)PI * fpeak * i * dt), 2));
			wavelet[i + nt / 2] = (float) (A * (1 - 2 * tmp) * exp(-tmp));
		}
	}

	if(mode == 2){
		int n1 = 0;
		bool flag = ISL_isOdd(nt);

		if (flag == false)
			n1 = nt - 1;
		else
			n1 = nt;

		float t = 0, T = 0;
		T = dt * (nt - 1) / 2;

		for (int i = 0; i < n1; i++) {
			t = -T + i * dt;
			wavelet[i] = (1 - 2 * (PI * fpeak * t) * (PI * fpeak * t)) * exp(-(PI
					* fpeak * t) * (PI * fpeak * t));
		}

		if (flag == false)
			wavelet[nt - 1] = wavelet[nt - 2];
	}

	if(mode == 3){
		int it;
		float t1, t0;

		t0 = 1.0 / fpeak;

		for (it = 0; it < nt; it++) {
			t1 = it * dt;
			wavelet[it] = exp(-PI * PI * fpeak * fpeak * (t1 - t0) * (t1 - t0))
					* (1.0 - 2. * PI * PI * fpeak * fpeak * (t1 - t0) * (t1 - t0));
		}
	}
}


/**********************************************************************************************************************
 *
 *  功能：生成宽带雷克子波
 *
 *  说明：本函数用于生成宽带雷克子波，零相位的
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		float				dt					In			采样间隔
 *		int					fLowCut				In			低截频率
 *		int					fHighCut			In			高截频率
 *		float				nt					In			采样点数
 *		float*				wavelet				Out			生成的宽带雷克子波
 *
 *  返回：无
 *
**********************************************************************************************************************/
void ISL_broadBandRicker(int nt, float dt, float lowCut, float highCut, float *wavelet)
{
	float A = 1;
	float dTempHigh;
	float dTempLow;

	for (int i = -nt / 2; i <= nt / 2 - 1; i++) {
		dTempHigh = PI * highCut * i * dt * PI * highCut * i * dt;
		dTempLow = PI * lowCut * i * dt * PI * lowCut * i * dt;
		wavelet[i + nt / 2] = float(
				A * (highCut * exp(-dTempHigh) - lowCut * exp(-dTempLow)) / (highCut - lowCut));
	}
}


/**********************************************************************************************************************
 *
 *  功能：雷克子波
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					hlw					In			half-length of the wavelet including center (samples)
 *		float				dt					In			sampling interval (s)
 *		float				period				In			wavelet period (s)
 *		float				ampl				In			wavelet amplitude
 *		float				distort				In			wavelet distortion factor
 *
 *		float				wavelet				Out			Ricker wavelet
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_ricker3(int hlw, float dt, float period, float ampl, float distort, float *wavelet)
{
	float t = 0., tnorm, xsq, wsym;
	int j;

	for (j = 0; j < hlw; j++) {
		tnorm = t / period;
		xsq = (2.5 * tnorm) * (2.5 * tnorm);
		wsym = ampl * (1.0 - 2.0 * xsq) * exp(-xsq);
		wavelet[hlw + j - 1] = wsym * (1.0 - 2.0 * distort * tnorm);
		wavelet[hlw - j - 1] = wsym * (1.0 + 2.0 * distort * tnorm);
		t = t + dt;
	}
}


/*计算三参数小波,此函数用于单独计算三参数小波*/
float ISL_triParWave(float fTal, float fSigma, int nWav, float *fWavRe, float *fWavIm)
{
	int T0 = -4;      //母小波的初始坐标
	float dt = 0.0156f; //母波时间采样间隔
	float Bata = 0.;    //能量延迟因子
	float pi = 3.1415926f;

	memset(fWavRe, 0, 4 * nWav);
	memset(fWavIm, 0, 4 * nWav);

	float Tmp1, Tmp2;
	float aa, bb, cc, p1, q1, p, q, rk, ik, c1, c2, c3, d, CentFre, fTmp;
	cc = fSigma * fSigma / fTal;
	aa = float(sin(Bata * fSigma));
	bb = float(cos(Bata * fSigma));
	c1 = float(exp(-0.5 * cc));
	c2 = float(exp(-0.25 * cc));
	c3 = float(exp(-0.375 * cc));
	d = float(pow(2 * fTal / pi, 0.25f));

	p1 = float(4 * (c1 - c3) * bb * bb + 1 - c1);
	p = d * pow(p1, -0.5f);

	q1 = float(4 * (c1 - c3) * aa * aa + 1 - c1);
	q = d * pow(q1, -0.5f);

	rk = c2 * bb;
	ik = c2 * aa * q / p;
	CentFre = float(0.5 * sqrt(pi) * fSigma * p * q * (1 - c3) / sqrt(2 * fTal) / pi);

	/////计算实部和虚部/////
	float Pi25 = pow(pi, -0.25f);
	float TmpSum1 = 0., TmpSum2 = 0.;
	for (int i = 0; i < nWav; i++) {
		fTmp = float(T0 + dt * i);
		Tmp1 = exp(-fTal * (fTmp - Bata) * (fTmp - Bata));
		Tmp2 = fSigma * fTmp;

		TmpSum1 += float(Tmp1 * (p * (cos(Tmp2) - rk)) * Pi25);
		TmpSum2 += float(Tmp1 * (q * sin(Tmp2) + p * ik) * Pi25);

		fWavRe[i] = float(TmpSum1 * dt);
		fWavIm[i] = float(TmpSum2 * dt);
	}

	return CentFre;
}


/**********************************************************************************************************************
 *
 *  功能：akb子波
 *
 *  说明：Compute the time response of a source function as a wavelet based on a wavelet.
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					nt					In			number of samples in output wavelet
 *		float				dt					In			time step
 *		float				fpeak				In			peak frequency of the wavelet
 *
 *		float				wavelet				Out			array[nt] of computed wavelet
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_akb_wavelet(int nt, float dt, float fpeak, float *wavelet) {
	int it;
	float t1;
	float t0 = 1.8 / fpeak;

	for (it = 0; it < nt; it++) {
		t1 = it * dt;
		wavelet[it] = -60.0 * (t1 - t0) * exp(-2.0 * fpeak * fpeak * (t1 - t0)
				* (t1 - t0));
	}
}

/**********************************************************************************************************************
 *
 *  功能：尖脉冲子波
 *
 *  说明：Compute the time response of a source function as a spike.
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					nt					In			number of time step
 *		int					tindex				In			time index to locate the spike
 *
 *		float				wavelet				Out			array[nt] of computed wavelet
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_spike_wavelet(int nt, int tindex, float *wavelet) {
	int it;

	for (it = 0; it < nt; it++) {
		wavelet[it] = 0.0;
	}

	wavelet[tindex] = 1.0;
}

/**********************************************************************************************************************
 *
 *  功能：unit_wavelet
 *
 *  说明：Compute the	time response of a source function as a constant unit shift.
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					nt					In			number of samples in output wavelet
 *
 *		float				wavelet				Out			array[nt] of computed wavelet
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_unit_wavelet(int nt, float *wavelet) {
	int it;

	for (it = 0; it < nt; it++) {
		wavelet[it] = 1.0;
	}

}

/**********************************************************************************************************************
 *
 *  功能：零子波
 *
 *  说明：Compute the time response of a source function as zero everywhere.
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					nt					In			number of samples in output wavelet
 *
 *		float				wavelet				Out			array[nt] of computed wavelet
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_zero_wavelet(int nt, float *wavelet) {
	int it;

	for (it = 0; it < nt; it++) {
		wavelet[it] = 0.0;
	}
}

/**********************************************************************************************************************
 *
 *  功能：贝尔拉格子波
 *
 *  说明：Compute the time response of a source function as a
 *  Berlage wavelet with peak frequency "fpeak" Hz, exponential decay
 *  factor "decay", time exponent "tn", and initial phase angle "ipa".
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					nt					In			number of samples in output wavelet
 *		float				dt					In			time step
 *		float				fpeak				In			peak frequency of the Berlage wavelet
 *		float				ampl				In			wavelet amplitude
 *		float				tn					In			non-negative time exponent (typically an integer number)
 *		float				decay				In			non-negative exponential decay factor
 *		float				ipa					In			initial phase angle in radians
 *
 *		float				wavelet				Out			array[nt] of computed wavelet,Berlage wavelet
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_berlage_wavelet(int nt, float dt, float fpeak, float ampl, float tn,
							float decay, float ipa, float *wavelet)
{
	register int it;
	float t;

	for (it = 0; it < nt; it++) {
		t = dt * (float) it;
		wavelet[it] = ampl * pow(t, tn) * exp(-decay * t) * cos(2.0 * PI
				* fpeak * t + ipa);
	}
}

/**********************************************************************************************************************
 *
 *  功能：高斯子波
 *
 *  说明：Compute the time response of a source function as a
 *  	Gaussian wavelet with peak frequency "fpeak" in Hz.
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					nt					In			number of samples in output wavelet
 *		float				dt					In			time step
 *		float				fpeak				In			peak frequency of the Gaussian wavelet
 *
 *		float				wavelet				Out			array[nt] of computed wavelet, Gaussian wavelet
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_gaussian_wavelet(int nt, float dt, float fpeak, float *wavelet) {
	float t, t0, s;

	t0 = 1.0 / fpeak;
	s = 1.0 / (sqrt(2.0) * PI * fpeak);

	for (int it = 0; it < nt; it++) {
		t = dt * (float) it;
		wavelet[it] = 1.0 / (s * sqrt(2.0 * PI)) * exp(-(t - t0) * (t - t0)
				/ (2.0 * s * s));
	}
}

/**********************************************************************************************************************
 *
 *  功能：高斯一阶导数子波
 *
 *  说明：Compute the time response of a source function as a
 *  	Gaussian first derivative wavelet with peak frequency "fpeak" in Hz.
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					nt					In			number of samples in output wavelet
 *		float				dt					In			time step
 *		float				fpeak				In			peak frequency of the Gaussian wavelet
 *
 *		float				wavelet				Out			array[nt] of computed Gaussian wavelet, first derivative
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_gaussderiv_wavelet(int nt, float dt, float fpeak, float *wavelet) {
	register int it;
	float t, t0, s;

	t0 = 1.0 / fpeak;
	s = 1.0 / (sqrt(2.0) * PI * fpeak);

	for (it = 0; it < nt; it++) {
		t = dt * (float) it;
		wavelet[it] = -1.0 / (s * s * s * sqrt(2.0 * PI)) * (t - t0) * exp(-(t
				- t0) * (t - t0) / (2.0 * s * s));
	}
}

/**********************************************************************************************************************
 *
 *  功能：高斯一阶导数子波
 *
 *  说明：Compute n-th order derivative of a Gaussian in double precision.
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					nt					In			length of waveform in samples
 *		float				dt					In			sampling interval
 *		double				t0					In			time shift for (pseudo-) causality
 *		float				fpeak				In			maximum frequency
 *		int					n					In			order of derivative
 *		int					sign				In			multiplier for polarity of waveform
 *		int					verbose				In			flag for diagnostic messages
 *
 *		double				*w					Out			array of size nt containing the waveform
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_deriv_n_gauss(double dt, int nt, double t0, float fpeak, int n,
						double *w, int sign, int verbose)
{
	int i; /* loop variable			*/
	double sigma; /* temporal variance of Gaussian	*/
	double C; /* normalization constant		*/
	double *h = NULL; /* Hermite polynomial			*/
	double *h0 = NULL; /* temp array for H_{n-1}		*/
	double *h1 = NULL; /* temp array for H_{n}			*/
	double *t = NULL; /* time vector				*/

	/* allocate & initialize memory */
	t = alloc1double(nt);
	h = alloc1double(nt);
	h0 = alloc1double(nt);
	h1 = alloc1double(nt);

	memset((void *) t, 0, DSIZE * nt);
	memset((void *) h, 0, DSIZE * nt);
	memset((void *) h0, 0, DSIZE * nt);
	memset((void *) h1, 0, DSIZE * nt);
	if (verbose)
		fprintf(stderr, "memory allocated and initialized/n");

	/* initialize time vector */
	for (i = 0; i < nt; ++i)
		t[i] = i * dt - t0;
	if (verbose)
		fprintf(stderr, "t[] initialized/n");

	/* compute Gaussian */
	sigma = n / (4 * PI * PI * fpeak * fpeak);
	if (verbose)
		fprintf(stderr, "sigma=%f", sigma);
	for (i = 0; i < nt; ++i)
		w[i] = exp(-t[i] * t[i] / (2 * sigma));
	if (verbose)
		fprintf(stderr, "Gaussian computed/n");

	/* compute Hermite polynomial */
	for (i = 0; i < nt; ++i) {
		h0[i] = 1.0;
		h1[i] = t[i] / sigma;
	}
	if (n == 1)
		memcpy((void *) h, (const void *) h1, DSIZE * nt);
	if (n > 1)
		ISL_hermite_n_polynomial(h, h0, h1, t, nt, n, sigma);
	if (verbose)
		fprintf(stderr, "Hermite polynomial H_%d computed/n", n);

	/* compute waveform */
	for (i = 0; i < nt; ++i)
		w[i] = h[i] * w[i];
	if (verbose)
		fprintf(stderr, "waveform computed/n");

	/* find normalization constant */
	C = fabs(w[0]);
	for (i = 1; i < nt; ++i)
		if (fabs(w[i]) > C)
			C = fabs(w[i]);
	if (ISODD(n))
		C = -C; /* to account for (-1)^n */
	if (verbose)
		fprintf(stderr, "C=%f/n", C);

	/* and finally normalize */
	for (i = 0; i < nt; ++i)
		w[i] = sign * w[i] / C;
	if (verbose)
		fprintf(stderr, "waveform normalized/n");

	/* check amplitude a t=0 */
	if (verbose)
		fprintf(stderr, "w[o]=%.12f", w[0]);

	/* free memory */
	free1double(h);
	free1double(h0);
	free1double(h1);
	free1double(t);
	if (verbose)
		fprintf(stderr, "memory freed/n");
}

/**********************************************************************************************************************
 *
 *  功能：地震子波提取 1
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		float *				buf					In			输入的地震数据，数组形式
 *		int					nt					In			输入的地震数据长度，不能小于182
 *		int					k					In/Out		返回子波的长度
 *
 *  返回：float * ，返回数据的指针，得到后请自己负责释放内存空间,数组的大小为n
 *
 **********************************************************************************************************************/
float * ISL_extractSeisWavelet(float *buf, int nt, int &k) //求自相关---傅立叶变换----傅立叶逆变换
{
	int ntfft; // nt after padding for FFT
	int nw; // number of frequencies
	float *pfft = NULL; // float FFT workspace
	complex *cpfft = NULL; // complexex FFT workspace

	if (buf == NULL || nt < 182)
		return NULL;

	if (k < 91)
		k = 91;
	if (k > nt / 2)
		k = nt / 2;

	ntfft = su_fft::ISL_npfar(nt);
	nw = ntfft / 2 + 1;

	cpfft = new complex[nw];
	if (cpfft == NULL)
		return NULL;

	pfft = new float[ntfft];
	if (pfft == NULL) {
		delete[] cpfft;
		cpfft = NULL;
		return NULL;
	}

	ISL_xcor(nt, 0, buf, nt, 0, buf, nt, 0, pfft);

	for (int i = nt; i < ntfft; i++)
		pfft[i] = 0.0;

	su_fft::ISL_pfarc(1, ntfft, pfft, cpfft);

	ISL_normalizeFFT(cpfft, nw, 1);

	su_fft::ISL_pfacr(-1, ntfft, cpfft, pfft);

	//	estimate seismic wavelet
	ntfft = su_fft::ISL_npfar(k);
	nw = ntfft / 2 + 1;

	for (int i = k; i < ntfft; i++)
		pfft[i] = 0.0;

	su_fft::ISL_pfarc(1, ntfft, pfft, cpfft);

	ISL_arrangeWavelet(nw, cpfft, ntfft, pfft, k);

	return pfft;
}


/**********************************************************************************************************************
 *
 *  功能：求均方根差
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		float *				vals				In			输入的数组
 *		int					num					In			数组的长度
 *
 *  返回：float ，返回求均方根差
 *
 **********************************************************************************************************************/
float ISL_rms(float *vals, int num)/*求均方根差*/
{
	if (vals == NULL || num <= 0)
		return -9999.99;

	double rr = 0.0;
	for (int i = 0; i < num; i++) {
		rr = rr + (double) (vals[i] * vals[i]);
	}
	rr = rr / (double) num;

	float rms = sqrt(rr);

	return rms;
}


/**********************************************************************************************************************
 *
 *  功能：地震子波提取 2
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		float *				buf					In			输入的地震数据，数组形式
 *		float *				refs				In			输入的反射系数，数组形式
 *		int					nt					In			输入的地震数据和反射系数的长度，不能小于130
 *		int					k					In/Out		返回子波的长度
 *
 *  返回：float * ，返回数据的指针，得到后请自己负责释放内存空间,数组的大小为n
 *
 **********************************************************************************************************************/
float * ISL_extractSeisWavelet(float *buf, float *refs, int nt, int &k)//根据反射系数计算地震子波
{
	int ntfft; // nt after padding for FFT
	int nw; // number of frequencies
	float * pfft = NULL; // float FFT workspace
	complex *cpfft1 = NULL, *cpfft2 = NULL, *cpfft3 = NULL; // complexex FFT workspace
	complex cpaa, cpbb;

	if (buf == NULL || refs == NULL || nt < 130)
		return NULL;

	if (k < 65)
		k = 65;
	if (k > nt / 2)
		k = nt / 2;

	float rms1 = ISL_rms(refs, nt);
	float rms2 = ISL_rms(buf, nt);

	if (rms1 == 0.0 || rms2 == 0.0)
		return NULL;

	ntfft = su_fft::ISL_npfar(nt);
	nw = ntfft / 2 + 1;

	cpfft1 = new complex[nw];
	if (cpfft1 == NULL)
		return NULL;

	cpfft2 = new complex[nw];
	if (cpfft2 == NULL) {
		delete[] cpfft1;
		return NULL;
	}

	cpfft3 = new complex[nw];
	if (cpfft3 == NULL) {
		delete[] cpfft1;
		delete[] cpfft2;
		return NULL;
	}

	pfft = new float[ntfft];
	if (pfft == NULL) {
		delete[] cpfft1;
		delete[] cpfft2;
		delete[] cpfft3;
		return NULL;
	}

	for (int i = 0; i < nt; i++)
		pfft[i] = rms1 * (buf[i] / rms2);
	for (int i = nt; i < ntfft; i++)
		pfft[i] = 0.0;

	su_fft::ISL_pfarc(1, ntfft, pfft, cpfft1);

	for (int i = 0; i < nt; i++)
		pfft[i] = refs[i];
	for (int i = nt; i < ntfft; i++)
		pfft[i] = 0.0;

	su_fft::ISL_pfarc(1, ntfft, pfft, cpfft2);

	// calculate wavelet values

	/*
	 // 1. 频率域中直接使用除法

	 for(int i=0;i<nw;i++){
	 cpfft3[i] = cdiv(cpfft1[i],cpfft2[i]);
	 }

	 delete []cpfft1;
	 delete []cpfft2;

	 su_fft::ISL_pfacr(-1, ntfft, cpfft3, pfft);

	 float *vvs = ExtractSeisWavelet(pfft,nt,k);

	 delete []cpfft3;
	 delete pfft;

	 return vvs;
	 */

	// 2. 频率域中计算对数谱，再应用减法
	ntfft = su_fft::ISL_npfar(k);
	nw = ntfft / 2 + 1;

	for (int i = 0; i < nw; i++) {
		cpaa = ISL_cln(cpfft1[i]);
		cpbb = ISL_cln(cpfft2[i]);
		cpfft3[i] = ISL_csub(cpaa, cpbb);
	}

	delete[] cpfft1;
	delete[] cpfft2;

	ISL_arrangeWavelet(nw, cpfft3, ntfft, pfft, k);

	delete[] cpfft3;

	Smooth::ISL_smooth_WA(pfft, k, 5);
	Smooth::ISL_smooth_WA(pfft, k, 5);

	return pfft;
}


/**********************************************************************************************************************
 *
 *  功能：合成记录
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					wavelet_type		In			=1 for a spike
 *															=2 for Tong Fei's Ricker wavelet
 *															=3 for Larner's Ricker wavelet
 *															=4 for Akima wavelet
 *
 *		int					nx					In			number of horizontal samples (traces)
 *		int					nt					In			number of vertical samples (samples/trace)
 *		float				dt					In			time sampling interval
 *		float				fpeak				In			frequency peak for Ricker or Akima wavelets
 *
 *		float **			wfield				In/Out		[In]	array[nx][nt] of reflectivities
 *															[Out]	array[nx][nt] of seismic traces
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_convolve_wavelet(int wavelet_type, int nx, int nt, float dt, float fpeak, float **wfield)
{
	int ix, it; /* loop counters */
	float *wavelet; /* wavelet array */
	float *temp; /* scratch array */

	/* allocate working space */
	wavelet = alloc1float(nt);
	temp = alloc1float(nt);

	/* compute wavelet */
	if (wavelet_type == 1) {
		ISL_spike_wavelet(nt, 1, wavelet);
	} else if (wavelet_type == 2) {
		ISL_classicRicker(nt, dt, fpeak, wavelet, 3);
	} else if (wavelet_type == 3) {
		ISL_ricker3(nt / 2, dt, 1. / fpeak, 1., 0., wavelet);
	} else if (wavelet_type == 4) {
		ISL_akb_wavelet(nt, dt, fpeak, wavelet);
	} else {
		printf("wavelet_type has to be 1,2,3 or 4\n");
		return;
	}

	/* loop over traces to convolve with wavelet */
	for (ix = 0; ix < nx; ix++) {

		/* save input reflectivity trace */
		for (it = 0; it < nt; it++)
			temp[it] = wfield[ix][it];

		/* convolve with wavelet */
		ISL_conv(nt, 0, temp, nt, 0, wavelet, nt, 0, wfield[ix]);
	}

	/* free allocated space */
	free1float(temp);
	free1float(wavelet);
}



/**********************************************************************************************************************
 *
 *  功能：Bath 震源爆炸子波（单点）
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		float				A0					In
 *		float				r					In			爆破半径  m
 *		float				Vp					In			纵波速度   m/s
 *		float				z					In			传播距离  m
 *		float				t					In			传播时间  s
 *
 *  返回：float  单个杨点样点值
 *
 **********************************************************************************************************************/
float ISL_bathWavelet_SinglePoint(float A0, float r, float Vp, float z, float t)
{
	float B = (2 * Vp * t) / (3 * r);
	float C = sqrt((float) 2) * B;

	float expB = exp((float) -1 * B);
	float sinC = sqrt((float) 2) * sin(C);
	float cosC = r * cos(C);

	float func_1 = A0 * expB * ((z - r / (float) 2) * sinC - cosC);

	return func_1;
}



/**********************************************************************************************************************
 *
 *  功能：Bath 震源爆炸子波（数组）
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		float				A0					In
 *		float				r					In			爆破半径 (m)
 *		float				Vp					In			纵波速度  (m/s)
 *		float				z					In			传播距离  (m)
 *		int					nin					In			输入数组个数
 *		float				tin[]				In			输入的传播时间（ s）数组，个数为nin
 *		float				tout[]				out			输出的计算结果数组，个数为nin
 *
 *  返回：无
 *
 **********************************************************************************************************************/
void ISL_bathWavelet(float A0, float r, float Vp, float z, int nin, float tin[], float tout[])
{
	for (int i = 0; i < nin; i++) {
		tout[i] = ISL_bathWavelet_SinglePoint(A0, r, Vp, z, tin[i]);
	}
}

// =========================================================================================
/* 寻找子波的主频和相位 的 辅助函数 */
float ISL_lineSpace(float d1, float d2, int n, float *y)
{
	int n1 = floor(n) - 1;
	int nn = n - 1;
	float * vec = new float[nn];
	for (int i = 0; i < nn; i++)
		vec[i] = i;

	float d2md1 = d2 - d1;

	for (int i = 0; i < nn; i++)
		y[i] = d1 + vec[i] * d2md1 / n1;

	y[nn] = d2;

	if (vec) {
		delete[] vec;
		vec = NULL;
	}
	return 0;
}


/**********************************************************************************************************************
 *
 *  功能：计算子波的主屏与相位
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		float *				wdata				In			子波数据
 *		int					wnum				In			子波长度
 *		float				dt					In			时间采样间隔（秒）
 *		float				&freq				In			得到的主频(HZ)
 *		float				&ph					In			得到的相位
 *
 *  返回：0 ，表示一切正常；-1，表示程序异常
 *
 **********************************************************************************************************************/
int ISL_waveletFreqAndPh(float * wdata, int wnum, float dt/*秒*/, float &freq, float &ph)
{
	if (wdata == NULL || wnum <= 0)
		return -1;

	int nNumber = 1, nexp;
	for (int i = 1;; i++) {
		nNumber *= 2;
		nexp = i;
		if (nNumber >= wnum)
			break;
	}
	float *dataRe = new float[nNumber];
	float *dataIm = new float[nNumber];
	float *fftRe = new float[nNumber];
	float *fftIm = new float [nNumber];

	memset(dataRe, 0, sizeof(float) * nNumber);
	memset(dataIm, 0, sizeof(float) * nNumber);
	memset(fftRe, 0, sizeof(float) * nNumber);
	memset(fftIm, 0, sizeof(float) * nNumber);
	memcpy(dataRe, wdata, sizeof(float) * wnum);

	ISL_kfft(dataRe, dataIm, nNumber, nexp, fftRe, fftIm, 0, 0);

	int nf = nNumber / 2 + 1;
	float *freqs = new float[nf];

	float fnyq = 1.0 / (dt * 2.0);
	ISL_lineSpace(0, fnyq, nf, freqs);

	complex * ct = new complex[nf];
	for (int i = 0; i < nf; i++) {
		ct[i].r = fftRe[i];
		ct[i].i = fftIm[i];
	}

	int indexi = 0;
	for (int i = 1; i < nf; i++) {
		if (ISL_rcabs(ct[i]) > ISL_rcabs(ct[indexi]))
			indexi = i;
	}

	freq = freqs[indexi];
	//	ph = atan2(ct[indexi].i,ct[indexi].r);
	ph = fabs(atan2(ct[indexi].i, ct[indexi].r) / 3.1415926 * 180.0);
	if (ph < 1)
		ph = 0;
	if (ph > 360)
		ph = 360;

	if (ct) { delete[] ct;	ct = NULL; }
	if (dataRe) { delete[] dataRe;	dataRe = NULL; }
	if (dataIm) { delete[] dataIm;	dataIm = NULL; }
	if (fftRe) { delete[] fftRe;	fftRe = NULL; }
	if (fftIm) { delete[] fftIm;	fftIm = NULL; }
	if (freqs) { delete[] freqs;	freqs = NULL; }

	return 0;
}


/*高斯窗 Gauss Window*/
void ISL_gaussWindow(float * out, int n)
{
	float norm = 0;
	for (int i = 0; i < n; i++) {
		out[i] = float(
				exp(-0.5 * ((i - 0.5 * (n - 1)) / (0.5 * 0.5 * (n - 1)))
					* ((i - 0.5 * (n - 1)) / (0.5 * 0.5 * (n - 1)))));
		norm += out[i] * out[i];
	}

	norm = sqrt(norm);
	for (int i = 0; i < n; i++)
		out[i] /= norm;
}



/*汉宁窗 Hanning Window*/
void ISL_hanningWindow(float * out, int n)
{
	float norm = 0;
	for (int i = 0; i < n; i++) {
		out[i] = float(0.5 * (1 - cos(2 * (float)PI * i / (n - 1))));
		norm += out[i] * out[i];
	}

	norm = sqrt(norm);
	for (int i = 0; i < n; i++)
		out[i] /= norm;
}



/*海明窗  Hamming Window*/
void ISL_hammingWindow(float * out, int n)
{
	float norm = 0;
	for (int i = 0; i < n; i++) {
		out[i] = float(0.53836 - 0.46164 * cos(2 * (float) PI * i / (n - 1)));
		norm += out[i] * out[i];
	}

	norm = sqrt(norm);
	for (int i = 0; i < n; i++)
		out[i] /= norm;
}



/*Blackman*/
void ISL_blackmanWindow(float * out, int n)
{
	float norm = 0;
	for (int i = 0; i < n; i++) {
		out[i] = float(
				0.46 - 0.5 * cos(2 * (float)PI * i / (n - 1)) + 0.08
						* cos(4 * (float)PI * i / (n - 1)));
		norm += out[i] * out[i];
	}

	norm = sqrt(norm);
	for (int i = 0; i < n; i++)
		out[i] /= norm;
}

}/*End of ISLib*/
