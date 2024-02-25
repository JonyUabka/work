//***************************************************************************
// **
// **				ISeis Lib [Base Algorithm] - 平滑
// **
// **				Writer By	Mouri Song
// **							Baihong Liu
// **							Qiang Yang
// **
// **				(Soft Center IGP)
// **
// **				DATA : 2014-02-10
// **
// ****************************************************************************/

//#include "ISL_Smooth.h"
//#include "ISL_Interpolation.h"

//namespace ISLib
//{

//int Smooth::ISL_smooth_weights(int k1, int k2, int pos, int mid, float *ss,
//		float s0)
//{
//	if (k1 >= k2 || ss == NULL)
//		return 0;
//	//if(s0 <= 0.0 || s0 >= 1.0) s0 = 0.5;

//	float sx = (1.0 - s0) / 2.0;
//	int nn = k2 - k1 + 1;
//	if (s0 <= 0.0 || s0 >= 1.0)
//	{
//		for (int i = 0; i < nn; i++)
//		{
//			ss[i] = 1.0 / (float) nn;
//		}
//		return nn;
//	}

//	float tt;
//	if (k1 == (pos - mid) && k2 == (pos + mid))
//	{
//		tt = 0.0;
//		for (int i = 0; i < mid; i++)
//		{
//			ss[i] = (float) (i + 1) / (float) mid;
//			tt = tt + ss[i];
//		}
//		for (int i = 0; i < mid; i++)
//		{
//			ss[i] = sx * ss[i] / tt;
//			ss[nn - 1 - i] = ss[i];
//		}
//		ss[mid] = s0;
//		return nn;
//	}

//	int idx = pos - k1;
//	tt = 0.0;
//	for (int i = 0; i < idx; i++)
//	{
//		ss[i] = (float) (i + 1) / (float) idx;
//		tt = tt + ss[i];
//	}
//	for (int i = 0; i < idx; i++)
//	{
//		ss[i] = sx * ss[i] / tt;
//	}
//	tt = 0.0;
//	for (int i = nn - 1; i > idx; i--)
//	{
//		ss[i] = (float) (nn - i) / (float) (nn - idx);
//		tt = tt + ss[i];
//	}
//	for (int i = nn - 1; i > idx; i--)
//	{
//		ss[i] = sx * ss[i] / tt;
//	}
//	ss[idx] = sx;

//	return nn;
//}

///* 一维高斯平滑滤波器 */
//void Smooth::ISL_gaussian1d_smoothing(int ns, int nsr, float *data)
//{
//	int is; /* loop counter */
//	float sum = 0.0;
//	float fcut;
//	float r;
//	float fcutr = 1.0 / nsr;
//	static int n;
//	static int mean;
//	static float fcutl;
//	static float s[401]; /* smoothing filter array */
//	float *temp; /* temporary array */

//	/* allocate space */
//	temp = alloc1float(ns);

//	/* save input fcut */
//	fcut = fcutr;

//	/* don't smooth if nsr equal to zero */
//	if (nsr == 0 || ns <= 1)
//		return;

//	/* if halfwidth more than 100 samples, truncate */
//	if (nsr > 100)
//		fcut = 1.0 / 100;

//	/* initialize smoothing function if not the same as the last one used */
//	if (fcut != fcutl)
//	{
//		fcutl = fcut;

//		/* set span of 3, at width of 1.5*exp(-PI*1.5**2)=1/1174 */
//		n = 3.0 / fcut + 0.5;
//		n = 2 * n / 2 + 1; /* make it odd for symmetry */

//		/* mean is the index of the zero in the smoothing wavelet */
//		mean = n / 2;

//		/* s(n) is the smoothing gaussian */
//		for (is = 1; is <= n; is++)
//		{
//			r = is - mean - 1;
//			r = -r * r * fcut * fcut * PI;
//			s[is - 1] = exp(r);
//		}

//		/* normalize to unit area, will preserve DC frequency at full
//		 amplitude. Frequency at fcut will be half amplitude */
//		for (is = 0; is < n; is++)
//			sum += s[is];
//		for (is = 0; is < n; is++)
//			s[is] /= sum;
//	}
//	/* convolve by gaussian into buffer */
//	if (1.01 / fcutr > (float) ns)
//	{

//		/* replace drastic smoothing by averaging */
//		sum = 0.0;
//		for (is = 0; is < ns; is++)
//			sum += data[is];
//		sum /= ns;
//		for (is = 0; is < ns; is++)
//			data[is] = sum;

//	}
//	else
//	{

//		/* convolve with gaussian */
//		ISL_conv(n, -mean, s, ns, -mean, data, ns, -mean, temp);

//		/* copy filtered data back to output array */
//		for (is = 0; is < ns; is++)
//			data[is] = temp[is];
//	}

//	/* free allocated space */
//	free1float(temp);
//}

///* 高斯平滑直方图 */
//void Smooth::ISL_smooth_histogram(int nintlh, float *pdf)
//{
//	int i; /* loop counter */
//	int ng = 11; /* number of samples in smoothing filter */
//	int mg = 6; /* index of sample with zero lag */
//	float rv = 1.0; /* ?? */
//	float r; /* auxiliary variable */
//	float sum = 0.0; /* auxiliary variable */
//	float *filter = NULL; /* gaussian filter */
//	float *spdf = NULL; /* smoothed histogram */

//	/* allocate working space */
//	filter = alloc1float(ng);
//	spdf = alloc1float(nintlh);

//	/* compute gaussian filter */
//	for (i = 1; i <= ng; i++)
//	{

//		r = i - mg;
//		r = -r * r / rv / rv * PI;
//		filter[i - 1] = exp(r);
//	}

//	/* apply filter */
//	ISL_conv(ng, -mg + 1, filter, nintlh, 0, pdf, nintlh, 0, spdf);

//	/* normalize output histogram */
//	for (i = 0; i < nintlh; i++)
//		sum += spdf[i];
//	for (i = 0; i < nintlh; i++)
//		pdf[i] = spdf[i] / sum;

//	/* free allocated space */
//	free1float(filter);
//	free1float(spdf);
//}

///* 一维滚动窗口平均平滑法(生成滤波器，滚动波形) */
//void Smooth::ISL_rwa_smoothing_filter(int flag, int nl, int nr, float *filter)
//{
//	int i; /* loop counter */
//	int np = nl + nr + 1; /* number of filter points */

//	if (flag == 1)
//	{
//		float scale = 1.0 / np; /* scale for rectangular window */

//		for (i = 0; i < np; i++)
//			filter[i] = scale;

//	}
//	else if (flag == 2)
//	{
//		float scale = 0.0; /* scale for triangular window */

//		for (i = -nl; i < 0; i++)
//			filter[i + nl] = 1 + (float) i / nl;
//		for (i = 1; i < nr; i++)
//			filter[i + nl] = 1 - (float) i / nr;
//		filter[nl] = 1.0;

//		/* normalize */
//		for (i = 0; i < np; i++)
//			scale += filter[i];
//		for (i = 0; i < np; i++)
//			filter[i] /= scale;
//	}
//	else
//	{
//		printf("error in rwa_smoothing filter, flag should 1 or 2\n");
//		return;
//	}
//}

///* 9点平滑 */
//void Smooth::ISL_smooth_9_3(int n, float in[], float out[])
//{
//	int i;
//	if (n < 9)
//	{
//		for (i = 0; i <= n - 1; i++)
//			out[i] = in[i];
//	}
//	else
//	{
//		out[0] = 69.0 * in[0] + 4.0 * in[1] - 6.0 * in[2] + 4.0 * in[3] - in[4];
//		out[0] = out[0] / 70.0;
//		out[1] = 2.0 * in[0] + 27.0 * in[1] + 12.0 * in[2] - 8.0 * in[3];
//		out[1] = (out[1] + 2.0 * in[4]) / 35.0;
//		out[2] = -3.0 * in[0] + 12.0 * in[1] + 17.0 * in[2];
//		out[2] = (out[2] + 12.0 * in[3] - 3.0 * in[4]) / 35.0;
//		out[3] = -2.0 * in[0] + 3.0 * in[1] + 6.0 * in[2] + 7.0 * in[3];
//		out[3] = (out[3] + 6.0 * in[4] + 3.0 * in[5] - 2.0 * in[6]) / 21.0;

//		for (i = 4; i <= n - 5; i++)
//		{
//			out[i] = -21.0 * in[i - 4] + 14.0 * in[i - 3] + 39.0 * in[i - 2]
//					+ 54.0 * in[i - 1] + 59.0 * in[i];
//			out[i] = (out[i] + 54.0 * in[i + 1] + 39.0 * in[i + 2] + 14.0
//					* in[i + 3] - 21.0 * in[i + 4]) / 231.0;
//		}

//		out[n - 4] = -2.0 * in[n - 7] + 3.0 * in[n - 6] + 6.0 * in[n - 5] + 7.0
//				* in[n - 4];
//		out[n - 4] = (out[n - 4] + 6.0 * in[n - 3] + 3.0 * in[n - 2] - 2.0
//				* in[n - 1]) / 21.0;
//		out[n - 3] = -3.0 * in[n - 5] + 12.0 * in[n - 4] + 17.0 * in[n - 3];
//		out[n - 3] = (out[n - 3] + 12.0 * in[n - 2] - 3.0 * in[n - 1]) / 35.0;
//		out[n - 2] = 2.0 * in[n - 5] - 8.0 * in[n - 4] + 12.0 * in[n - 3];
//		out[n - 2] = (out[n - 2] + 27.0 * in[n - 2] + 2.0 * in[n - 1]) / 35.0;
//		out[n - 1] = -in[n - 5] + 4.0 * in[n - 4] - 6.0 * in[n - 3];
//		out[n - 1] = (out[n - 1] + 4.0 * in[n - 2] + 69.0 * in[n - 1]) / 70.0;
//	}
//	return;
//}

//void Smooth::ISL_smooth_5_3( int n, float in[], float out[] )
//{
//	int i;
//	if (n < 5) {
//		for (i = 0; i <= n - 1; i++)
//			out[i] = in[i];

//	} else {
//		out[0] = 69.0 * in[0] + 4.0 * in[1] - 6.0 * in[2] + 4.0 * in[3] - in[4];
//		out[0] = out[0] / 70.0;
//		out[1] = 2.0 * in[0] + 27.0 * in[1] + 12.0 * in[2] - 8.0 * in[3];
//		out[1] = (out[1] + 2.0 * in[4]) / 35.0;

//		for (i = 2; i <= n - 3; i++) {
//			out[i] = -3.0 * in[i - 2] + 12.0 * in[i - 1] + 17.0 * in[i];
//			out[i] = (out[i] + 12.0 * in[i + 1] - 3.0 * in[i + 2]) / 35.0;
//		}
//		out[n - 2] = 2.0 * in[n - 5] - 8.0 * in[n - 4] + 12.0 * in[n - 3];
//		out[n - 2] = (out[n - 2] + 27.0 * in[n - 2] + 2.0 * in[n - 1]) / 35.0;
//		out[n - 1] = -in[n - 5] + 4.0 * in[n - 4] - 6.0 * in[n - 3];
//		out[n - 1] = (out[n - 1] + 4.0 * in[n - 2] + 69.0 * in[n - 1]) / 70.0;
//	}
//	return;
//}

///* 多点平滑 */
//void Smooth::ISL_smooth_MuliPoints(float *vlc, int ns, int spn, int cpn)
//{
//	int i, j, k, n, m1, m2;
//	float *ss = NULL, Sum;

//	n = cpn;
//	if (n < 1)
//		n = 1;

//	m1 = (int) (spn * 0.5);
//	if (m1 < 1)
//		m1 = 1;
//	m2 = m1 * 2 + 1;

//	ss = new float[ns + m2];

//	for (i = 0; i < n; i++)
//	{
//		for (j = 0; j < m1; j++)
//			ss[j] = vlc[0];
//		for (j = 0; j < ns; j++)
//			ss[j + m1] = vlc[j];
//		for (j = 0; j < m1; j++)
//			ss[j + ns + m1] = vlc[ns - 1];
//		for (j = 0; j < ns; j++)
//		{
//			Sum = 0.;
//			for (k = 0; k < m2; k++)
//				Sum = Sum + ss[j + k];
//			vlc[j] = Sum / m2;
//		}
//	}
//	delete[] ss;
//}

///* 高斯离散点一维平滑 */
//void Smooth::ISL_smooth_disc(int nin, float xin[], float yin[], int intr_num,
//		float intr_idx[], int c)
//{
//	float * nout = new float[intr_num];

//	Interpolation::ISL_intlin(nin, xin, yin, 0, 0, intr_num, intr_idx, nout, 1);

//	float max, min, top;
//	ISL_findMaxValue(nout, intr_num, max, min, top);

//	float * nouto = new float[intr_num];
//	//	smooth_9_3(intr_num, nout, nouto);
//	ISL_gaussian1d_smoothing(intr_num, c, nout);

//	for (int i = 0; i < nin; i++)
//	{
//		for (int j = 0; j < intr_num; j++)
//		{
//			if (ISL_floatEqual(xin[i], intr_idx[j]) == 1)
//			{
//				yin[i] = nout[j] * max * 0.5;
//				continue;
//			}
//		}
//	}
//	if (nout)
//	{
//		delete[] nout;
//		nout = NULL;
//	}
//	if (nouto)
//	{
//		delete[] nouto;
//		nouto = NULL;
//	}
//}

///* 加权平均平滑法 */
//int Smooth::ISL_smooth_WA(float *vals, int num, int count, float s0)
//{
//	if (vals == NULL)
//		return 0;
//	if (count > 65)
//		count = 65;
//	//if(s0 <= 0.0 || s0 >= 1.0) s0 = 0.5;

//	int nn;
//	nn = count / 2;
//	if (nn <= 0)
//		nn = 1;
//	nn = 2 * nn + 1;

//	if (num <= nn)
//		return 0;

//	float *outs = new float[num];
//	if (outs == NULL)
//		return 0;

//	int j, jj, k, k1, k2;
//	float ss[65];
//	for (int i = 0; i < num; i++)
//	{
//		outs[i] = 0.0;
//		k1 = i - nn / 2;
//		if (k1 < 0)
//			k1 = 0;
//		k2 = i + nn / 2;
//		if (k2 > (num - 1))
//			k2 = num - 1;
//		outs[i] = 0.0;
//		jj = Smooth::ISL_smooth_weights(k1, k2, i, nn / 2, ss, s0);

//		if (jj <= 0)
//		{
//			outs[i] = vals[i];
//			continue;
//		}
//		k = k1;

//		for (j = 0; j < jj; j++)
//		{
//			if (k > k2)
//				break;
//			outs[i] = outs[i] + ss[j] * vals[k];
//			k++;
//		}
//	}

//	for (int i = 0; i < num; i++)
//		vals[i] = outs[i];

//	delete[] outs;
//	outs = NULL;

//	return nn;
//}

///* 帽子平滑法 */
//void Smooth::ISL_smooth_CAP(int num, float *yin, float *yout, int step,
//		int flag)
//{
//	float *xout = new float[num];
//	float *rangeArray = new float[step];

//	int iPeak = num / step;
//	float *peakArray = new float[iPeak];

//	float *xin = new float[iPeak];
//	int iCount = 0;
//	for (int i = 0; i < iPeak; i++)
//	{

//		for (int j = 0; j < step; j++)
//		{
//			rangeArray[j] = yin[iCount];
//			iCount++;
//			xout[i * step + j] = iCount;
//			yout[i * step + j] = 0;
//		}

//		float max = 0, min = 0, top = 0;
//		ISL_findMaxValue(rangeArray, step, max, min, top);
//		peakArray[i] = max;

//		int idx = ISL_findValue(rangeArray, step, max);
//		if (idx != -1)
//		{
//			xin[i] = iCount - step + idx;
//		}
//		else
//		{
//			xin[i] = iCount - step;
//		}
//	}

//	float max = 0, min = 0, top = 0;
//	ISL_findMaxValue(peakArray, iPeak, max, min, top);

//	Interpolation::ISL_intlin(iPeak, xin, peakArray, min, min, num, xout, yout);

//	if (flag != 2)
//	{
//		float smoothArray[num];
//		Smooth::ISL_smooth_9_3(num, yout, smoothArray);
//		for (int i = 0; i < num; i++)
//		{
//			yout[i] = smoothArray[i];
//		}
//	}
//	else
//	{
//		Smooth::ISL_smooth_MuliPoints(yout, num, 2, 2);
//	}

//	if(xout){ delete []xout; xout = NULL; }
//	if(xin){ delete []xin; xin = NULL; }
//	if(rangeArray){ delete []rangeArray; rangeArray = NULL; }
//	if(peakArray){ delete []peakArray; peakArray = NULL; }
//}

//}/*End of ISLib
