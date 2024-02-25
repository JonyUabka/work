/**
*	@file	ISL_Filter.cpp
*	@brief	[Source file of Filter Functions], 滤波；
*	@see	ISeisLib Manual
*	@author [Liu Baihong, Yang Qiang, Song ZhiXiang], 刘百红、杨强、宋志翔；
*	@date	2014-06-03
*	@refer	SU CWP
*/

#include "ISL_Filter.h"
#include "ISL_Matrix.h"
#include "ISL_Statistics.h"
#include "ISL_FFT.h"
#include "ISL_Alloc.h"
#include "ISL_Interpolation.h"
#include "ISL_Convolution.h"

namespace ISLib {

/*提取滤波因子*/
int ISL_filterFactor(float f1, float f2, float f3, float f4, float sim,
							int il, float *&filbuf, float scl)
{
	if (filbuf) {
		delete[] filbuf;
		filbuf = NULL;
	}

	if (f1 >= f2)
		return (-1);

	filbuf = new float[il];
	if (filbuf == NULL)
		return -1;

	int i, i2;
	float a, b, c, d, t, fa, fb, fc, fd, df1, df2;
	float xdou, ydou;

	d = sim * scl;
	i2 = il / 2;
	filbuf[i2] = 1.0;
	i2 = i2 + 1;
	a = f2 - f1;
	b = f2 + f1;
	fa = 0.5 * a;
	fb = 4.0 / (a * a);
	t = d;

	for (i = 0; i < i2; i++) {
		a = PI * t;
		xdou = a * fa;
		ydou = (float) sin(xdou);
		c = ydou / a;
		c = c * c;
		a = a * b;
		df1 = fb * c / a;
		xdou = a;
		ydou = (float) sin(xdou);
		filbuf[i2 - (i + 1)] = df1 * ydou;
		t = t + d;
	}

	if (f3 > 0.0 && f4 > 0.0 && f3 > f2) {
		df2 = b;
		a = f4 - f3;
		b = f4 + f3;
		fa = 0.5 * a;
		fb = 4.0 / (a * a);
		t = d;
		df1 = b - df2;

		for (i = 0; i < i2; i++) {
			a = PI * t;
			xdou = a * fa;
			ydou = (float) sin(xdou);
			c = ydou / a;
			c = c * c;
			xdou = a * b;
			ydou = (float) sin(xdou);
			fc = ydou / a;
			fd = fb * fc * c;
			filbuf[i2 - (i + 1)] = (fd - filbuf[i2 - (i + 1)] * df2) / df1;
			t = t + d;
		}
	}

	d = 0.0;
	for (i = 0; i < i2; i++) {
		filbuf[il - (i + 1)] = filbuf[i];
		d = d + 2.0 * filbuf[i];
	}

	d = 1.0 / d;
	for (i = 0; i < il; i++) {
		filbuf[i] = filbuf[i] * d;
	}

	return 0;
}


/*滤波或褶积*/
int ISL_filter(float *filter, int len, float *orig, int num, float *out, int flag)
{
	if (filter == NULL || len <= 0)
		return -1;
	if (orig == NULL || num <= 0)
		return -1;

	float *ftmp = NULL;
	if (!out) {
		ftmp = new float[num];
		if (ftmp == NULL)
			return -1;
	}

	int i, j;
	int n = len / 2;
	float aa;
	for (i = 0; i < num; i++) {
		aa = 0.0;
		if (out)
			out[i] = 0.0;
		else
			ftmp[i] = 0.0;
		for (j = -n; j < n; j++) {
			if (i - j < 0 || i - j >= num)
				continue;
			if (orig[i - j] == 0.0 || filter[j + n] == 0.0)
				continue;
			if (out)
				out[i] = out[i] + orig[i - j] * filter[j + n];
			else
				ftmp[i] = ftmp[i] + orig[i - j] * filter[j + n];
			aa = aa + filter[j + n];
		}
		if (flag == 1) {
			if (out)
				out[i] = out[i] / aa;
			else
				ftmp[i] = ftmp[i] / aa;
		}
	}

	if (!out) {
		for (i = 0; i < num; i++)
			orig[i] = ftmp[i];
		delete ftmp;
	}

	return 0;
}


/* 带通滤波 */
float * ISL_bandpassFilter(float *vals, int num, float dtime, float f1, float f2, float f3, float f4, int len, float scl)
{
	if (vals == NULL || num <= 0)
		return NULL;
	if (dtime <= 0.0 || len <= 0 || scl <= 0.0)
		return NULL;
	if (f1 >= f2)
		return NULL;

	float *filbuf = NULL;
	int status = ISL_filterFactor(f1, f2, f3, f4, dtime, len, filbuf, scl);
	if (status != 0)
		return NULL;

	float *outs = new float[num];
	if (outs == NULL) {
		delete[] filbuf;
		filbuf = NULL;
		return NULL;
	}

	status = ISL_filter(filbuf, len, vals, num, outs);
	delete filbuf;

	if (status != 0) {
		delete[] outs;
		outs = NULL;
		return NULL;
	}

	return outs;
}

/**********************************************************************************************************************
 *
 *  功能：巴特沃斯滤波器 (Butterworth filter)
 *
 *  说明： compute number of poles and -3 db frequency
 *  for a low-pass or high-pass filter, given a frequency response
 *  constrained at two frequencies.
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		float				fpass				In			frequency in pass band at which amplitude is >= apass
 *		float				apass				In			amplitude in pass band corresponding to frequency fpass
 *		float				fstop				In			frequency in stop band at which amplitude is <= astop
 *		float				astop				In			amplitude in stop band corresponding to frequency fstop
 *
 *		float *				npoles				Out			number of poles
 *		float *				f3db				Out			frequency at which amplitude is sqrt(0.5) (-3 db)
 *
 *  返回：无
 *  (1) Nyquist frequency equals 0.5
 *
 *  (2) The following conditions must be true:
 *  	(0.0<fpass && fpass<0.5) &&
 *  	(0.0<fstop && fstop<0.5) &&
 *  	(fpass!=fstop) &&
 *  	(0.0<astop && astop<apass && apass<1.0)
 *
 *  (3) if (fpass<fstop)
 *  		a low-pass filter is assumed
 *  	else
 *  		a high-pass filter is assumed
 *
**********************************************************************************************************************/
void ISL_bfdesign (float fpass, float apass, float fstop, float astop, int &npoles, float &f3db)
{
	float wpass,wstop,fnpoles,w3db;

	/* warp frequencies according to bilinear transform */
	wpass = 2.0*tan(PI*fpass);
	wstop = 2.0*tan(PI*fstop);

	/* if lowpass filter, then */
	if (fstop>fpass) {
		fnpoles = log((1.0/(apass*apass)-1.0)/(1.0/(astop*astop)-1.0))
			/ log(pow(wpass/wstop,2.0));
		w3db = wpass/pow((1.0/(apass*apass)-1.0),0.5/fnpoles);

	/* else, if highpass filter, then */
	} else {
		fnpoles = log((1.0/(apass*apass)-1.0)/(1.0/(astop*astop)-1.0))
			/ log(pow(wstop/wpass,2.0));
		w3db = wpass*pow((1.0/(apass*apass)-1.0),0.5/fnpoles);
	}

	/* determine integer number of poles */
	npoles = 1+(int)fnpoles;

	/* determine (unwarped) -3 db frequency */
	f3db = atan(0.5*w3db)/PI;
}


/*巴特沃斯-高通*/
void ISL_bfhighpass (int npoles, float f3db, int n, float p[], float q[])
{
	int jpair,j;
	float r,scale,theta,a,b1,b2,pj,pjm1,pjm2,qjm1,qjm2;

	r = 2.0*tan(PI*fabs(f3db));
	if (npoles%2!=0) {
		scale = r+2.0;
		a = 2.0/scale;
		b1 = (r-2.0)/scale;
		pj = 0.0;
		qjm1 = 0.0;
		for (j=0; j<n; j++) {
			pjm1 = pj;
			pj = p[j];
			q[j] = a*(pj-pjm1)-b1*qjm1;
			qjm1 = q[j];
		}
	} else {
		for (j=0; j<n; j++)
			q[j] = p[j];
	}
	for (jpair=0; jpair<npoles/2; jpair++) {
		theta = PI*(2*jpair+1)/(2*npoles);
		scale = 4.0+4.0*r*sin(theta)+r*r;
		a = 4.0/scale;
		b1 = (2.0*r*r-8.0)/scale;
		b2 = (4.0-4.0*r*sin(theta)+r*r)/scale;
		pjm1 = 0.0;
		pj = 0.0;
		qjm2 = 0.0;
		qjm1 = 0.0;
		for (j=0; j<n; j++) {
			pjm2 = pjm1;
			pjm1 = pj;
			pj = q[j];
			q[j] = a*(pj-2.0*pjm1+pjm2)-b1*qjm1-b2*qjm2;
			qjm2 = qjm1;
			qjm1 = q[j];
		}
	}
}


/*巴特沃斯-低通*/
void ISL_bflowpass (int npoles, float f3db, int n, float p[], float q[])
{
	int jpair, j;
	float r, scale, theta, a, b1, b2, pj, pjm1, pjm2, qjm1, qjm2;

	r = 2.0 * tan(PI * fabs(f3db));
	if (npoles % 2 != 0) {
		scale = r + 2.0;
		a = r / scale;
		b1 = (r - 2.0) / scale;
		pj = 0.0;
		qjm1 = 0.0;
		for (j = 0; j < n; j++) {
			pjm1 = pj;
			pj = p[j];
			q[j] = a * (pj + pjm1) - b1 * qjm1;
			qjm1 = q[j];
		}
	} else {
		for (j = 0; j < n; j++)
			q[j] = p[j];
	}
	for (jpair = 0; jpair < npoles / 2; jpair++) {
		theta = PI * (2 * jpair + 1) / (2 * npoles);
		scale = 4.0 + 4.0 * r * sin(theta) + r * r;
		a = r * r / scale;
		b1 = (2.0 * r * r - 8.0) / scale;
		b2 = (4.0 - 4.0 * r * sin(theta) + r * r) / scale;
		pjm1 = 0.0;
		pj = 0.0;
		qjm2 = 0.0;
		qjm1 = 0.0;
		for (j = 0; j < n; j++) {
			pjm2 = pjm1;
			pjm1 = pj;
			pj = q[j];
			q[j] = a * (pj + 2.0 * pjm1 + pjm2) - b1 * qjm1 - b2 * qjm2;
			qjm2 = qjm1;
			qjm1 = q[j];
		}
	}
}


/* 离散随机线性系统的卡尔曼滤波  */
int ISL_lman(int n, int m, int k, double *f, double *q, double *r, double *h, double *y,
		double *x, double *p, double *g)
{
	int i, j, kk, ii, l, jj, js = 0;
	double *e, *a, *b;

	e = new double[m * m];
	l = m;
	if (l < n)
		l = n;

	a = new double[l * l];
	b = new double[l * l];

	for (i = 0; i <= n - 1; i++)
		for (j = 0; j <= n - 1; j++) {
			ii = i * l + j;
			a[ii] = 0.0;
			for (kk = 0; kk <= n - 1; kk++)
				a[ii] = a[ii] + p[i * n + kk] * f[j * n + kk];
		}
	for (i = 0; i <= n - 1; i++)
		for (j = 0; j <= n - 1; j++) {
			ii = i * n + j;
			p[ii] = q[ii];
			for (kk = 0; kk <= n - 1; kk++)
				p[ii] = p[ii] + f[i * n + kk] * a[kk * l + j];
		}
	for (ii = 2; ii <= k; ii++) {
		for (i = 0; i <= n - 1; i++)
			for (j = 0; j <= m - 1; j++) {
				jj = i * l + j;
				a[jj] = 0.0;
				for (kk = 0; kk <= n - 1; kk++)
					a[jj] = a[jj] + p[i * n + kk] * h[j * n + kk];
			}
		for (i = 0; i <= m - 1; i++)
			for (j = 0; j <= m - 1; j++) {
				jj = i * m + j;
				e[jj] = r[jj];
				for (kk = 0; kk <= n - 1; kk++)
					e[jj] = e[jj] + h[i * n + kk] * a[kk * l + j];
			}
		js = ISL_rinv(e, m);
		if (js == 0) {
			if(e){ delete []e; e = NULL; }
			if(a){ delete []a; a = NULL; }
			if(b){ delete []b; b = NULL; }
			return (js);
		}
		for (i = 0; i <= n - 1; i++)
			for (j = 0; j <= m - 1; j++) {
				jj = i * m + j;
				g[jj] = 0.0;
				for (kk = 0; kk <= m - 1; kk++)
					g[jj] = g[jj] + a[i * l + kk] * e[j * m + kk];
			}
		for (i = 0; i <= n - 1; i++) {
			jj = (ii - 1) * n + i;
			x[jj] = 0.0;
			for (j = 0; j <= n - 1; j++)
				x[jj] = x[jj] + f[i * n + j] * x[(ii - 2) * n + j];
		}
		for (i = 0; i <= m - 1; i++) {
			jj = i * l;
			b[jj] = y[(ii - 1) * m + i];
			for (j = 0; j <= n - 1; j++)
				b[jj] = b[jj] - h[i * n + j] * x[(ii - 1) * n + j];
		}
		for (i = 0; i <= n - 1; i++) {
			jj = (ii - 1) * n + i;
			for (j = 0; j <= m - 1; j++)
				x[jj] = x[jj] + g[i * m + j] * b[j * l];
		}
		if (ii < k) {
			for (i = 0; i <= n - 1; i++)
				for (j = 0; j <= n - 1; j++) {
					jj = i * l + j;
					a[jj] = 0.0;
					for (kk = 0; kk <= m - 1; kk++)
						a[jj] = a[jj] - g[i * m + kk] * h[kk * n + j];
					if (i == j)
						a[jj] = 1.0 + a[jj];
				}
			for (i = 0; i <= n - 1; i++)
				for (j = 0; j <= n - 1; j++) {
					jj = i * l + j;
					b[jj] = 0.0;
					for (kk = 0; kk <= n - 1; kk++)
						b[jj] = b[jj] + a[i * l + kk] * p[kk * n + j];
				}
			for (i = 0; i <= n - 1; i++)
				for (j = 0; j <= n - 1; j++) {
					jj = i * l + j;
					a[jj] = 0.0;
					for (kk = 0; kk <= n - 1; kk++)
						a[jj] = a[jj] + b[i * l + kk] * f[j * n + kk];
				}
			for (i = 0; i <= n - 1; i++)
				for (j = 0; j <= n - 1; j++) {
					jj = i * n + j;
					p[jj] = q[jj];
					for (kk = 0; kk <= n - 1; kk++)
						p[jj] = p[jj] + f[i * n + kk] * a[j * l + kk];
				}
		}
	}
	if(e){ delete []e; e = NULL; }
	if(a){ delete []a; a = NULL; }
	if(b){ delete []b; b = NULL; }

	return (js);
}

/**********************************************************************************************************************
 *
 *  功能：alpha-beta-gamma 滤波
 *
 *  说明：
 *
 *  参数：
 *		Type				Name				In/Out		Description
 *		----				----				------		-----------
 *		int					n					In			量测数据的点数
 *		double *			x[n]				In			N个等间隔点上的量测值
 *		double				t					In			采样周期
 *		double				a					In			滤波器结构参数-alpha
 *		double				b					In			滤波器结构参数-beta
 *		double				c					In			滤波器结构参数-gamma
 *		double *			y[n]				Out			返回n个等间隔点上的滤波估值
 *
 *  返回：无
 *
**********************************************************************************************************************/
void ISL_kabg(int n, double *x, double t, double a, double b, double c, double *y)
{
	int i;
	double s1, ss, v1, vv, a1, aa;
	aa = 0.0;
	vv = 0.0;
	ss = 0.0;
	for (i = 0; i <= n - 1; i++) {
		s1 = ss + t * vv + t * t * aa / 2.0;
		v1 = vv + t * aa;
		a1 = aa;
		ss = s1 + a * (x[i] - s1);
		y[i] = ss;
		vv = v1 + b * (x[i] - s1);
		aa = a1 + 2.0 * c * (x[i] - s1) / (t * t);
	}
	return;
}



/*计算中值滤波,取一定的时窗对取出的数据进行滤波*/
void ISL_windowMidValue(float *inData, int num, int mfl, float *outData)
{
	if (mfl % 2 == 0)
		mfl++;

	int halfWin;
	float medium;
	float *dataTimeWin = new float[mfl];

	halfWin = mfl / 2;
	for (int i = 0; i < num; i++) {
		if (i < halfWin) {
			for (int j = i - halfWin; j <= i + halfWin; j++) {
				if (j < 0)
					dataTimeWin[j - i + halfWin] = inData[2 * i - j];
				else
					dataTimeWin[j - i + halfWin] = inData[j];
			}
		} else if (i >= num - halfWin) {
			for (int j = i - halfWin; j <= i + halfWin; j++) {
				if (j >= num)
					dataTimeWin[j - i + halfWin] = inData[2 * i - j];
				else
					dataTimeWin[j - i + halfWin] = inData[j];
			}
		} else {
			for (int j = i - halfWin; j <= i + halfWin; j++)
				dataTimeWin[j - i + halfWin] = inData[j];
		}

		/////中值滤波/////
		ISL_midValue(dataTimeWin, mfl, medium);
		outData[i] = medium;
	}

	delete[] dataTimeWin;
	dataTimeWin = NULL;
}

// ===================================================================================================
//				Advance filter
// ===================================================================================================


/** @brief	Dip Filter */
int ISL_dipFilter(	int nt, int ntr, float ** inData,
					float dt, float dtr,
					int nslopes, float *slopes, float *amps,
					float bias,
					int &ntfft, int &ntrfft, float **& outData )
{
	if(nslopes != 4 || slopes == NULL || amps == NULL)
		return -1;

	float sfft; 	/* scale factor for FFT */
	int nw; 		/* number of frequencies */
	float dw; 		/* frequency sampling interval */
	float fw; 		/* first frequency */
	int nk; 		/* number of wavenumbers */
	float dk; 		/* wavenumber sampling interval */
	float w, k; 	/* frequency and wavenumber */
	int it, ix, iw, ik; /* sample indices */
	float slope, amp; 	/* slope and amplitude for particular w,k */
	complex **cpfft; 	/* complex FFT workspace */
	float phase; 			/* phase shift for bias */
	complex cshift; 		/* complex phase shifter for bias */

	/* determine lengths and scale factors for prime-factor FFTs */
	ntfft = su_fft::ISL_npfar(nt);
	ntrfft = su_fft::ISL_npfa(ntr);
	sfft = 1.0/(ntfft*ntrfft);

	/* determine frequency and wavenumber sampling */
	nw = ntfft / 2 + 1;
	dw = 2.0 * PI / (ntfft * dt);
	fw = 0.000001 * dw; /* non-zero to avoid divide by zero w */
	nk = ntrfft;
	dk = 2.0 * PI / (ntrfft * dtr);

	/* allocate real and complex workspace for FFTs */
	cpfft = alloc2complex(nw, nk);
	if(outData){
		free2float(outData);
		outData = NULL;
	}
	outData = alloc2float(ntfft, ntrfft);

	/*zero all arrays*/
	memset((void *) outData[0], 0, ntfft * ntrfft * FSIZE);

	for (ix = 0; ix < ntrfft; ix++)
		for (it = 0; it < ntfft; it++)
			outData[ix][it] = 0.0;
	for (ix = 0; ix < ntr; ix++)
		for (it = 0; it < nt; it++)
			outData[ix][it] = inData[ix][it];

	/* Fourier transform t to w */
	su_fft::ISL_pfa2rc(1, 1, ntfft, ntr, outData[0], cpfft[0]);

	/* do linear moveout bias via phase shift */
	for (ix = 0; ix < ntr; ix++) {
		for (iw = 0, w = 0.0; iw < nw; iw++, w += dw) {
			phase = -ix * dtr * w * bias;
			cshift = ISL_cmplx(cos(phase), sin(phase));
			cpfft[ix][iw] = ISL_cmul(cpfft[ix][iw], cshift);
		}
	}

	/* Fourier transform x to k */
	su_fft::ISL_pfa2cc(-1, 2, nw, ntrfft, cpfft[0]);

	/* loop over wavenumbers */
	for (ik = 0; ik < nk; ik++) {

		/* determine wavenumber */
		k = (ik <= nk / 2) ? ik * dk : (ik - nk) * dk;

		/* loop over frequencies */
		for (iw = 0, w = fw; iw < nw; iw++, w += dw) {

			/* determine biased slope */
			slope = k / w + bias;

			/* linearly interpolate to find amplitude */
			Interpolation::ISL_intlin(nslopes, slopes, amps, amps[0], amps[nslopes - 1], 1, &slope, &amp);

			/* include fft scaling */
			amp *= sfft;

			/* filter real and imaginary parts */
			cpfft[ik][iw].r *= amp;
			cpfft[ik][iw].i *= amp;
		}
	}

	/* Fourier transform k to x */
	su_fft::ISL_pfa2cc(1, 2, nw, ntrfft, cpfft[0]);

	/* undo linear moveout bias via phase shift */
	for (ix = 0; ix < ntr; ix++) {
		for (iw = 0, w = 0.0; iw < nw; iw++, w += dw) {
			phase = ix * dtr * w * bias;
			cshift = ISL_cmplx(cos(phase), sin(phase));
			cpfft[ix][iw] = ISL_cmul(cpfft[ix][iw], cshift);
		}
	}

	/* Fourier transform w to t */
	su_fft::ISL_pfa2cr(-1, 1, ntfft, ntr, cpfft[0], outData[0]);

	/* clear memory */
	free2complex(cpfft);
	cpfft = NULL;

	return 0;
}

// =========================================================================================
// =========================================================================================
void polygonalFilter( float *f, float *amps, int npoly, int nfft, float dt, float *filter )
/*************************************************************************
polygonalFilter -- polygonal filter with sin tapering
**************************************************************************
Input:
f		array[npoly] of frequencies defining the filter
amps		array[npoly] of amplitude values
npoly		size of input f and amps arrays
dt		time sampling interval
nfft		number of points in the fft
pow		power flag, use pow=1 for 1D usage, =2 for 2D

Output:
filter		array[nfft] filter values
**************************************************************************
Notes: Filter is to be applied in the frequency domain
Also, this is not exactly the same filter that is used in "sufilter".
That filter is sin^2. This done because the filter is squared in
its application to make the amplitudes come out right.
**************************************************************************
Author:  CWP: John Stockwell   1992
*************************************************************************/
#define PIBY2   1.57079632679490
{
	int *intfr; /* .... integerizations of f		*/
	int icount, ifs; /* loop counting variables              */
	int taper = 0; /* flag counter				*/
	int nf; /* number of frequencies (incl Nyq)     */
	int nfm1; /* nf-1                                 */
	float onfft; /* reciprocal of nfft                   */
	float df; /* frequency spacing (from dt)          */

	intfr = alloc1int(npoly);

	nf = nfft / 2 + 1;
	nfm1 = nf - 1;
	onfft = 1.0 / nfft;

	/* Compute array of integerized frequencies that define the filter*/
	df = onfft / dt;
	for (ifs = 0; ifs < npoly; ++ifs) {
		intfr[ifs] = NINT(f[ifs] / df);
		if (intfr[ifs] > nfm1)
			intfr[ifs] = nfm1;
	}

	/* Build filter, with scale, and taper specified by amps[] values*/
	/* Do low frequency end first*/
	for (icount = 0; icount < intfr[0]; ++icount)
		filter[icount] = amps[0] * onfft;

	/* now do the middle frequencies */
	for (ifs = 0; ifs < npoly - 1; ++ifs) {
		if (amps[ifs] < amps[ifs + 1]) {
			++taper;
			for (icount = intfr[ifs]; icount <= intfr[ifs + 1]; ++icount) {
				float c = PIBY2 / (intfr[ifs + 1] - intfr[ifs] + 2);
				float s = sin(c * (icount - intfr[ifs] + 1));
				float adiff = amps[ifs + 1] - amps[ifs];
				filter[icount] = (amps[ifs] + adiff * s) * onfft;
			}
		} else if (amps[ifs] > amps[ifs + 1]) {
			++taper;
			for (icount = intfr[ifs]; icount <= intfr[ifs + 1]; ++icount) {
				float c = PIBY2 / (intfr[ifs + 1] - intfr[ifs] + 2);
				float s = sin(c * (intfr[ifs + 1] - icount + 1));
				float adiff = amps[ifs] - amps[ifs + 1];
				filter[icount] = (amps[ifs + 1] + adiff * s) * onfft;
			}
		} else if (!(taper)) {
			for (icount = intfr[ifs]; icount <= intfr[ifs + 1]; ++icount)
				filter[icount] = amps[ifs] * onfft;
		} else {
			for (icount = intfr[ifs] + 1; icount <= intfr[ifs + 1]; ++icount)
				filter[icount] = amps[ifs] * onfft;
		}
	}

	/* finally do the high frequency end */
	for (icount = intfr[npoly - 1] + 1; icount < nf; ++icount) {
		filter[icount] = amps[npoly - 1] * onfft;
	}
}


void polartorect (	int na, float da, float fa, int nr, float dr, float fr, float **q,
					int nx, float dx, float fx, int ny, float dy, float fy, float **p )
/*****************************************************************************
Convert a function of q(a,r) to p(x,y), where x = r*cos(a) and y = r*sin(a)
******************************************************************************
Input:
na		number of a samples
da		a sampling interval
fa		first a sample
nr		number of r samples
dr		r sampling interval
fr		first r sample
nx		number of x samples
dx		x sampling interval
fx		first x sample
ny		number of y samples
dy		y sampling interval
fy		first y sample
q		array[nr][na] containing samples of q(a,r)

Output:
p		array[ny][nx] containing samples of p(x,y)
******************************************************************************
Notes:
The polar angle a is measured in radians.

Linear extrapolation is used to determine the value of q(a,r) for
a and r coordinates not in the range corresponding to na, da, ....
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/15/90
******************************************************************************/
{
	int ix, iy, ia, ir;
	float x, y, a = 0.0, r, ai, ri, sa, sr;

	/* for all y */
	for (iy = 0, y = fy; iy < ny; ++iy, y += dy) {

		/* for all x */
		for (ix = 0, x = fx; ix < nx; ++ix, x += dx) {

			/* determine a and r */
			if (x != 0.0)
				a = atan2((double) y, (double) x);
			else if (y > 0.0)
				a = PI / 2.0;
			else if (y < 0.0)
				a = -PI / 2.0;
			else if (y == 0.0)
				a = 0.0;

			r = sqrt(x * x + y * y);

			/* determine sample indices */
			ai = (a - fa) / da;
			ia = ai;
			if (ia < 0)
				ai = ia = 0;
			if (ia > na - 2) {
				ai = na - 1;
				ia = na - 2;
			}
			ri = (r - fr) / dr;
			ir = ri;
			if (ir < 0)
				ri = ir = 0;
			if (ir > nr - 2) {
				ri = nr - 1;
				ir = nr - 2;
			}

			/* bilinear interpolation */
			sa = ai - ia;
			sr = ri - ir;
			p[iy][ix] = (1.0 - sr)
					* ((1.0 - sa) * q[ir][ia] + sa * q[ir][ia + 1])
					+ sr
							* ((1.0 - sa) * q[ir + 1][ia]
									+ sa * q[ir + 1][ia + 1]);
		}
	}
}

// =======================================================================
int ISL_kFilter(	int nt, int ntr, float ** inData,
					float dt, float dtr,
					int knum, float *k, float *amps,
					int nphi,
					int &ntfft, int &ntrfft, float ** & outData )
{

	if(knum != 4 || k == NULL || amps == NULL)
		return -1;

	float dx; 			/* sampling intervals */
	int nK1, nK2, nK; 	/* transform (output) dimensions	*/
	int nxfft; 			/* dimensions after padding for FFT	*/
	int iphi;			/* phi values				*/

    int ik; 			/* k counter				*/
	int ik1;			/* k1 counter				*/
	int ik2;			/* k2 counter				*/

	complex ** ct; 		/* complex FFT workspace		*/
	float ** kphi;		/* k,phi array				*/
	float * kfilt; 		/* k wavenumber filter			*/
	float ** karray;	/* array of filter values 		*/
	float ** kfilter;	/* wavenumber filter 			*/


	dx = NINT(sqrt(dt * dt + dtr * dtr));

	/* Determine lengths for prime-factor FFTs */
	ntfft = su_fft::ISL_npfaro(nt, LOOKFAC * nt);
	ntrfft = su_fft::ISL_npfa(ntr);

	if (ntfft >= PFA_MAX || ntrfft >= PFA_MAX)
		return -1;

	/* Determine number of wavenumbers in K1 and K2 */
	nK1 = ntfft / 2 + 1;
	nK2 = ntrfft / 2 + 1;
	nK = NINT(sqrt(nK1 * nK1 + nK2 * nK2));
	nxfft = 2 * (nK - 1); /* faked for polygonal filter in k */

	// ====================================================
	/* Allocate space */
	if(outData){
		free2float(outData);
		outData = NULL;
	}
	outData = alloc2float(ntfft, ntrfft);
	ct = alloc2complex(nK1, ntrfft);
	karray = alloc2float(ntfft, ntrfft);
	kfilter = alloc2float(ntfft, ntrfft);
	kfilt = alloc1float(nK);
	kphi = alloc2float(1, nK);

	/* Zero all arrays */
	memset((void *) outData[0], 0, ntfft*ntrfft*FSIZE);
	memset((void *) kfilter[0], 0, ntfft * ntrfft * FSIZE);
	memset((void *) ct[0], 0, nK1 * ntrfft * sizeof(complex));
	memset((void *) kfilt, 0, nK * FSIZE);
	memset((void *) kphi[0], 0, nK * FSIZE);

	/* Build Filter */
	polygonalFilter(k, amps, knum, nxfft, dx, kfilt);

	/* Build k,phi filter */
	for (iphi = 0; iphi < 1; ++iphi)
		for (ik = 0; ik < nK; ++ik)
			kphi[ik][iphi] = kfilt[ik] * kfilt[ik];

	/* convert polar to rectangular coordinates */
	polartorect(nphi, 0, 0, nK, 1, 0, kphi, ntfft, 1, -nK1, ntrfft, 1, -nK2, karray);

	/* positive k1, positive k2 */
	for (ik2 = 0; ik2 < nK2; ++ik2)
		for (ik1 = 0; ik1 < nK1; ++ik1)
			kfilter[ik2][ik1] = karray[nK2 - 1 - ik2][nK1 - 1 - ik1];

	/* positive k1, negative k2 */
	for (ik2 = nK2; ik2 < ntrfft; ++ik2)
		for (ik1 = 0; ik1 < nK1; ++ik1)
			kfilter[ik2][ik1] = karray[ik2 - nK2][nK1 - 1 - ik1];


	/* set Output Data value */
	for(int i=0; i<ntr; i++){
		for(int j=0; j<nt; j++){
			outData[i][j] = inData[i][j];
		}
	}

	/* Fourier transform dimension 1 */
	su_fft::ISL_pfa2rc(-1, 1, ntfft, ntr, outData[0], ct[0]);

	/* Fourier transform dimension 2 */
	su_fft::ISL_pfa2cc(-1, 2, nK1, ntrfft, ct[0]);

	/* Apply filter */
	for (ik2 = 0; ik2 < ntrfft; ++ik2)
		for (ik1 = 0; ik1 < nK1; ++ik1)
			ct[ik2][ik1] = ISL_crmul(ct[ik2][ik1], kfilter[ik2][ik1]);

	/* Inverse Fourier transform dimension 2 */
	su_fft::ISL_pfa2cc(1, 2, nK1, ntrfft, ct[0]);

	/* Inverse Fourier transform dimension 1 */
	su_fft::ISL_pfa2cr(1, 1, ntfft, ntr, ct[0], outData[0]);

	/* clear memory */
	free2complex(ct);
	ct = NULL;

	free2float(kphi);
	kphi = NULL;

	free1float(kfilt);
	kfilt = NULL;

	free2float(karray);
	karray = NULL;

	free2float(kfilter);
	kfilter = NULL;

	return 0;
}


/** @brief	Wiener predictive error filtering */
void ISL_pefFilter(	int nt, float dt, float * inData,
					float * & outData )
{
	int i, ilag, imix; 	/* counters	*/
	float pnoise = PNOISE;

	// ===============================
	float minlag, maxlag;
	int iminlag, imaxlag;

	minlag = 0.04;
	maxlag = 0.16;

	iminlag = NINT(minlag / dt);
	imaxlag = NINT(maxlag / dt);

	// ===============================
	float mincorr;		/* start time of correlation window	*/
	int imincorr;		/* .. in samples	*/
	float maxcorr;		/* end time of correlation window	*/
	int imaxcorr;		/* .. in samples	*/

	mincorr = 0.0;
	imincorr = NINT(mincorr/dt);

	maxcorr = 1.0;
	imaxcorr = NINT(maxcorr/dt);

	float mix = 1.0;
	int nmix = 1;

	/* compute filter sizes and correlation number */
	int nlag = imaxlag - iminlag + 1;
	int ncorr = imaxcorr - imincorr + 1;
	int lcorr = imaxlag + 1;

	/* Compute byte sizes in wiener/spiker and autocorr */
	size_t lagbytes = FSIZE * nlag;
	size_t maxlagbytes = FSIZE * lcorr;
	size_t mixbytes = maxlagbytes * nmix;

	/* Allocate memory */
	float * wiener = alloc1float(nlag);
	float * spiker = alloc1float(nlag);
	float * autocorr = alloc1float(lcorr);
	float * temp = alloc1float(lcorr);
	float ** mixacorr = alloc2float(lcorr, nmix);

	/* Set pointer to "cross" correlation */
	float * crosscorr = autocorr + iminlag;

	/* Zero out mixing array */
	memset((void *) mixacorr[0], 0, mixbytes);

	/* zero out filter vectors */
	memset((void *) wiener, 0, lagbytes);
	memset((void *) spiker, 0, lagbytes);
	memset((void *) autocorr, 0, maxlagbytes);
	memset((void *) temp, 0, maxlagbytes);

	// ================================================  开始计算  ========================================================
	/* Form autocorrelation vector */
	ISL_xcor(	ncorr, imincorr, inData,
				ncorr, imincorr, inData,
				lcorr, 0, autocorr);

	/* Whiten */
	autocorr[0] *= 1.0 + pnoise;

	/* Read autocorr into first column of mixacorr[][] */
	memcpy((void *) mixacorr[0], (const void *) autocorr, maxlagbytes);

	/* Loop over values of the autocorrelation array */
	for (ilag = 0; ilag < lcorr; ++ilag) {
		/* Weighted moving average (mix) */
		temp[ilag] += mixacorr[0][ilag] * mix;

		/* put mixed data back in seismic trace */
		autocorr[ilag] = temp[ilag];
	}

	/* to make space for autocorr from next trace */
	for (imix = nmix - 1; 0 < imix; --imix)
		for (ilag = 0; ilag < lcorr; ++ilag)
			mixacorr[imix][ilag] = mixacorr[imix - 1][ilag];

	/* Get inverse filter by Wiener-Levinson */
	ISL_stoepf(nlag, autocorr, crosscorr, wiener, spiker);

	/* Convolve pefilter with trace - don't do zero multiplies */
	for (i = 0; i < nt; ++i) {
		register int j;
		register int n = MIN(i, imaxlag);
		register float sum = inData[i];

		for (j = iminlag; j <= n; ++j)
			sum -= wiener[j - iminlag] * outData[i - j];

		outData[i] = sum;
	}

	free1float(wiener);
	free1float(spiker);
	free1float(autocorr);
	free1float(temp);
	free2float(mixacorr);
}


// ===================================================================================================
//				applies a zero-phase, sine-squared tapered filter
// ===================================================================================================
void ISL_taperFilter(	int nt, float dt, float * inData,
						float fnum, float *f, float *amps,
						float * & outData	)
{
    float *rt = NULL;		/* real trace                       */
    complex *ct = NULL;		/* complex transformed trace        */
    float *filter = NULL;	/* filter array                     */
    float nyq;              /* nyquist frequency                */
    int nfft;               /* number of points for fft trace   */
    int nf;                 /* number of frequencies (incl Nyq) */

    /* Set up FFT parameters */
	nfft = su_fft::ISL_npfaro(nt, LOOKFAC * nt);
	nf = nfft/2 + 1;

    nyq = 0.5/dt;

    /* Allocate fft arrays */
	rt = alloc1float( nfft );
	ct = alloc1complex( nf );
	filter = alloc1float( nf );

	/* Build the polygonal filter filter[]*/
	polygonalFilter(f, amps, fnum, nfft, dt, filter);

    /* Load trace into rt (zero-padded) */
	memcpy((void *) rt, (const void *) inData, nt * sizeof(float));
	memset((void *) (rt + nt), 0, (nfft - nt) * sizeof(float));

    /* FFT, filter, inverse FFT */
	su_fft::ISL_pfarc(1, nfft, rt, ct);
	for (int i = 0; i < nf; ++i)
		ct[i] = ISL_crmul(ct[i], filter[i]);
	su_fft::ISL_pfacr(-1, nfft, ct, rt);

    /* Load traces back in, recall filter had nfft factor */
	for (int i = 0; i < nt; ++i) {
		outData[i] = rt[i];
	}

	free1float(rt);
	free1float(filter);
	free1complex(ct);
}

// ===================================================================================================
//				TV Band
// ===================================================================================================
void bandpass(float *data, int nt, int nfft, int nfreq, float *filterj, float *ftracej)
{
	float *rt = (float*) malloc(sizeof(float) * nfft);
	complex *ct = (complex*) malloc(sizeof(complex) * nfreq);

	/* Load trace into rt (zero-padded) */
	memcpy((void *) rt, (const void *) data, nt * FSIZE);
	memset((void *) (rt + nt), 0, (nfft - nt) * FSIZE);

	/* FFT, filter, inverse FFT */
	su_fft::ISL_pfarc(1, nfft, rt, ct);
	for (int i = 0; i < nfreq; ++i)
		ct[i] = ISL_crmul(ct[i], filterj[i]);
	su_fft::ISL_pfacr(-1, nfft, ct, rt);

	/* Load traces back in, recall filter had nfft factor */
	for (int i = 0; i < nt; ++i)
		ftracej[i] = rt[i]; /* ftracej = rt ?? */

	/* free temp arrays */
	free(rt);
	free(ct);
}


void makefilter(float *f, int nfft, int nfreq, float dt, float *filter)
{
	float onfft = 1.0 / nfft;
	float df = onfft / dt;
	int nfreqm1 = nfreq - 1;
	int if1 = NINT(f[0] / df);
	int if2 = NINT(f[1] / df);
	int if3 = MIN(NINT(f[2]/df), nfreqm1);
	int if4 = MIN(NINT(f[3]/df), nfreqm1);

	// ====================================================
	float c = PIBY2 / (if2 - if1 + 2);
	for (int i = if1; i <= if2; ++i) {
		float s = sin(c * (i - if1 + 1));
		filter[i] = s * s * onfft;
	}

	// ====================================================
	float c1 = PIBY2 / (if4 - if3 + 2);
	for (int i = if3; i <= if4; ++i) {
		float s = sin(c1 * (if4 - i + 1));
		filter[i] = s * s * onfft;
	}

	// ====================================================
	for (int i = if2 + 1; i < if3; ++i)
		filter[i] = onfft;
	for (int i = 0; i < if1; ++i)
		filter[i] = 0.0;
	for (int i = if4 + 1; i < nfreq; ++i)
			filter[i] = 0.0;
}

void ISL_tvBandFilter(	int nt, float dt, float * inData,
						int nfilter, float *tf, float **f,
						float * & outData	)
{
	float *temp_tf = NULL;	/* times at which filters are centered  */
	float **temp_f = NULL;	/* 4 corner frequency array for each times */
	int *itf = NULL;		/* ... as integers                      */
	float ** ftrace = NULL;	/* filtered sub-traces                  */
	float ** filter = NULL;         /* filter arrays               	*/

	int jmin;				/* index of first filter itf value	*/
	int jmax;				/* index of last filter itf value	*/
	float tmin = 0;			/* first time on traces				*/
    int nfft;     	        /* fft sizes in each time gate		*/
    int nfreq;     	        /* number of frequencies			*/

	temp_tf = alloc1float(nfilter + 4);
	itf = alloc1int(nfilter + 4);
	int iCount = 0;
	for(int i = 0; i<nfilter + 4; i++){
		temp_tf[i] = 0.;
		itf[i] = 0;
		if(i>=2 && i<(nfilter + 2)){
			temp_tf[i] = tf[iCount];
			iCount++;
		}
	}

	jmin = 2;
	jmax = nfilter + 1;
	for (int i = jmin; i <= jmax; ++i){
	  	itf[i] = NINT((temp_tf[i] - tmin)/dt);
	}

    /* Make filters with scale for inverse transform */
	nfft = su_fft::ISL_npfaro(nt, LOOKFAC * nt);
	nfreq = nfft / 2 + 1;

	// ================================================
	filter = alloc2float(nfreq, (nfilter+4));
	temp_f = alloc2float(nfreq, (nfilter+4));
	iCount = 0;
	for(int i = 0; i<nfilter + 4; i++){
		temp_f[i][0] = temp_f[i][1] = temp_f[i][2] = temp_f[i][3] = 0;
		if(i>=2 && i<(nfilter + 2)){
			temp_f[i][0] = f[iCount][0];
			temp_f[i][1] = f[iCount][1];
			temp_f[i][2] = f[iCount][2];
			temp_f[i][3] = f[iCount][3];
			iCount++;
		}
	}

	for (int j = jmin; j <= jmax; ++j) {
		makefilter(temp_f[j], nfft, nfreq, dt, filter[j]);
	}

	/* User may not have given a filter for tmin and/or tmax--	*/
	/* Extend array so can always assume these filters are present.	*/
	/* Note don't really use any of the extra storage in **filter!	*/
	if (itf[jmin] > 0) {
		filter[jmin - 1] = filter[jmin];
		itf[jmin - 1] = 0;
		--jmin;
	}
	if (itf[jmax] < nt - 1) {
		filter[jmax + 1] = filter[jmax];
		itf[jmax + 1] = nt - 1;
		++jmax;
	}

	/* Extend array so can always consider time points to be interior */
	itf[jmin-1] = 0;      /* now jmin - 1 is a valid index */
	itf[jmax+1] = nt - 1; /* now jmax + 1 is a valid index */

	/* Main loop over traces */
	ftrace = alloc2float(nfreq, (nfilter+4));

	/* Construct filtered sub-traces */
	for (int j = jmin; j <= jmax; ++j) {
		bandpass(inData, nt, nfft, nfreq, filter[j], ftrace[j]);
	}

    /* Compose filtered trace from sub-traces */
	int ii;
	for (int j = jmin; j < jmax; ++j) {
		float fitfj;
		for (fitfj = ii = itf[j]; ii <= itf[j + 1]; ++ii) {
			float a = (ii - fitfj) / (itf[j + 1] - fitfj);
			outData[ii] = (1 - a) * ftrace[j][ii] + a * ftrace[j + 1][ii];
		}
	}

	free1float(temp_tf);
	free1int(itf);
	free2float(ftrace);
	free2float(filter);
	free2float(temp_f);
}

} /* End Of namespace ISLIB */

